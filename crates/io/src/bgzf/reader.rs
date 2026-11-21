use std::collections::VecDeque;
use std::fmt;
use std::fs::File;
use std::io::{BufReader, ErrorKind, Read, Seek, SeekFrom};
use std::thread::{self, JoinHandle};

use crossbeam_channel::{Receiver, Sender, bounded, never};
use libdeflater::Decompressor;
use rayon::prelude::*;

use super::block::{Block, FIXED_HEADER_LEN, Frame, block_len, extra_len, has_magic};
use crate::{Error, Workers};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Readahead {
	batch_bytes: usize,
	batches: usize,
}

impl Readahead {
	pub const DEFAULT: Self = Self {
		batch_bytes: 16 << 20,
		batches: 4,
	};
	pub const COMPACT: Self = Self {
		batch_bytes: 1 << 20,
		batches: 1,
	};
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockSpan {
	pub first: u64,
	pub last: u64,
}

pub enum Input {
	Stream(Box<dyn Read + Send>),
	Spans { file: File, spans: Vec<BlockSpan> },
}

impl fmt::Debug for Input {
	fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
		match self {
			Self::Stream(_) => formatter.write_str("Stream"),
			Self::Spans { spans, .. } => formatter
				.debug_struct("Spans")
				.field("spans", spans)
				.finish_non_exhaustive(),
		}
	}
}

#[derive(Debug)]
pub struct BlockReader {
	batches: Receiver<Result<Vec<Block>, Error>>,
	ready: VecDeque<Block>,
	threads: Vec<JoinHandle<()>>,
	workers: Workers,
}

impl BlockReader {
	pub fn new(input: Input, workers: &Workers) -> Self {
		Self::with_readahead(input, workers, Readahead::DEFAULT)
	}

	pub fn with_readahead(input: Input, workers: &Workers, readahead: Readahead) -> Self {
		let (frames_out, frames_in) = bounded(readahead.batches);
		let (blocks_out, blocks_in) = bounded(readahead.batches);
		let scanner = thread::spawn(move || scan(input, &frames_out, readahead));
		let inflating = workers.clone();
		let inflater = thread::spawn(move || inflate(&frames_in, &blocks_out, &inflating));
		Self {
			batches: blocks_in,
			ready: VecDeque::new(),
			threads: vec![scanner, inflater],
			workers: workers.clone(),
		}
	}

	pub fn next_batch(&mut self) -> Option<Result<Vec<Block>, Error>> {
		if !self.ready.is_empty() {
			return Some(Ok(self.ready.drain(..).collect()));
		}
		self.batches.recv().ok()
	}

	pub fn next_block(&mut self) -> Option<Result<Block, Error>> {
		loop {
			if let Some(block) = self.ready.pop_front() {
				return Some(Ok(block));
			}
			match self.batches.recv().ok()? {
				Ok(batch) => self.ready.extend(batch),
				Err(error) => return Some(Err(error)),
			}
		}
	}
}

impl Drop for BlockReader {
	fn drop(&mut self) {
		self.batches = never();
		for thread in self.threads.drain(..) {
			let _ = self.workers.join(thread);
		}
	}
}

impl Iterator for BlockReader {
	type Item = Result<Block, Error>;

	fn next(&mut self) -> Option<Self::Item> {
		self.next_block()
	}
}

fn scan(input: Input, out: &Sender<Result<Vec<Frame>, Error>>, readahead: Readahead) {
	let buffer = readahead.batch_bytes.min(8 << 20);
	let result = match input {
		Input::Stream(stream) => scan_stream(
			&mut BufReader::with_capacity(buffer, stream),
			0,
			u64::MAX,
			out,
			readahead,
		),
		Input::Spans { file, spans } => scan_spans(
			&mut BufReader::with_capacity(buffer, file),
			&spans,
			out,
			readahead,
		),
	};
	if let Err(error) = result {
		let _ = out.send(Err(error));
	}
}

fn scan_spans(
	file: &mut BufReader<File>,
	spans: &[BlockSpan],
	out: &Sender<Result<Vec<Frame>, Error>>,
	readahead: Readahead,
) -> Result<(), Error> {
	for span in spans {
		file.seek(SeekFrom::Start(span.first))?;
		scan_stream(file, span.first, span.last, out, readahead)?;
	}
	Ok(())
}

fn scan_stream(
	stream: &mut impl Read,
	mut offset: u64,
	last: u64,
	out: &Sender<Result<Vec<Frame>, Error>>,
	readahead: Readahead,
) -> Result<(), Error> {
	let mut batch = Vec::new();
	let mut batch_bytes = 0;
	while offset <= last {
		let Some(frame) = read_frame(stream, offset)? else {
			break;
		};
		offset += frame.bytes.len() as u64;
		batch_bytes += frame.bytes.len();
		batch.push(frame);
		if batch_bytes >= readahead.batch_bytes {
			out.send(Ok(std::mem::take(&mut batch)))
				.map_err(|_| Error::Disconnected)?;
			batch_bytes = 0;
		}
	}
	if !batch.is_empty() {
		out.send(Ok(batch)).map_err(|_| Error::Disconnected)?;
	}
	Ok(())
}

fn read_frame(stream: &mut impl Read, offset: u64) -> Result<Option<Frame>, Error> {
	let mut header = [0; FIXED_HEADER_LEN];
	if !fill(stream, &mut header, offset)? {
		return Ok(None);
	}
	if !has_magic(&header) {
		return Err(Error::NotBgzf { offset });
	}
	let extra = extra_len(&header);
	let mut bytes = Vec::with_capacity(FIXED_HEADER_LEN + extra);
	bytes.extend_from_slice(&header);
	bytes.resize(FIXED_HEADER_LEN + extra, 0);
	read_exact(stream, &mut bytes[FIXED_HEADER_LEN..], offset)?;
	let len = block_len(&bytes, offset)?;
	let known = bytes.len();
	if len < known {
		return Err(Error::Corrupt { offset });
	}
	bytes.resize(len, 0);
	read_exact(stream, &mut bytes[known..], offset)?;
	Ok(Some(Frame { offset, bytes }))
}

fn fill(stream: &mut impl Read, buffer: &mut [u8], offset: u64) -> Result<bool, Error> {
	let mut filled = 0;
	while filled < buffer.len() {
		match stream.read(&mut buffer[filled..]) {
			Ok(0) if filled == 0 => return Ok(false),
			Ok(0) => return Err(Error::Truncated { offset }),
			Ok(read) => filled += read,
			Err(error) if error.kind() == ErrorKind::Interrupted => {}
			Err(error) => return Err(error.into()),
		}
	}
	Ok(true)
}

fn read_exact(stream: &mut impl Read, buffer: &mut [u8], offset: u64) -> Result<(), Error> {
	stream
		.read_exact(buffer)
		.map_err(|error| match error.kind() {
			ErrorKind::UnexpectedEof => Error::Truncated { offset },
			_ => error.into(),
		})
}

fn inflate(
	frames: &Receiver<Result<Vec<Frame>, Error>>,
	out: &Sender<Result<Vec<Block>, Error>>,
	workers: &Workers,
) {
	for batch in frames {
		let blocks = batch.and_then(|frames| {
			workers.install(|| {
				frames
					.into_par_iter()
					.map_init(Decompressor::new, |decompressor, frame| {
						frame.inflate(decompressor)
					})
					.collect::<Result<Vec<_>, _>>()
			})
		});
		let failed = blocks.is_err();
		if out.send(blocks).is_err() || failed {
			return;
		}
	}
}
