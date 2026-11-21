use std::fs::File;
use std::io::Read;
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::thread;

use gyrase_io::{Block, BlockReader, BlockSpan, Input, VirtualOffset, Workers};

use super::batch::{Located, RecordBatch};
use super::frame::Framer;
use super::query::{Chunk, IntervalSet, Query};
use crate::error::{Error, HeaderError};
use crate::header::{ByteSource, Header};

const HEADER_WINDOWS: [u64; 3] = [1 << 20, 1 << 24, u64::MAX];
const BATCHES_READY: usize = 2;
const BATCHES_RECYCLED: usize = 4;

#[derive(Debug)]
enum Filter {
	Everything,
	Intervals(IntervalSet),
}

#[derive(Debug)]
struct BatchReader {
	blocks: BlockReader,
	leftover: Vec<Block>,
	framer: Framer,
	filter: Filter,
	reference_count: u32,
}

#[derive(Debug)]
struct Pipeline {
	ready: Receiver<Result<RecordBatch, Error>>,
	recycle: SyncSender<RecordBatch>,
}

#[derive(Debug)]
pub struct BamReader {
	header: Header,
	batches: Pipeline,
	current: RecordBatch,
	cursor: usize,
}

impl BamReader {
	pub fn open(source: Box<dyn Read + Send>, workers: &Workers) -> Result<Self, Error> {
		let mut blocks = BlockReader::new(Input::Stream(source), workers);
		let mut cursor = HeaderCursor::new(&mut blocks);
		let header = Header::read_bam(&mut cursor)?;
		let (leftover, start) = cursor.into_rest();
		let chunk = Chunk {
			start,
			end: VirtualOffset::from_raw(u64::MAX),
		};
		let batches = BatchReader::new(
			blocks,
			leftover,
			vec![chunk],
			Filter::Everything,
			reference_count(&header),
		);
		Ok(Self::assemble(header, batches.spawn()))
	}

	pub fn query(
		file: File,
		header: Header,
		query: Query,
		workers: &Workers,
	) -> Result<Self, Error> {
		let (chunks, intervals) = query.into_parts();
		let spans = chunks
			.iter()
			.map(|chunk| BlockSpan {
				first: chunk.start.block(),
				last: chunk.end.block(),
			})
			.collect();
		let blocks = BlockReader::new(Input::Spans { file, spans }, workers);
		let batches = BatchReader::new(
			blocks,
			Vec::new(),
			chunks,
			Filter::Intervals(intervals),
			reference_count(&header),
		);
		Ok(Self::assemble(header, batches.spawn()))
	}

	pub fn read_header(file: &File, workers: &Workers) -> Result<Header, Error> {
		for last in HEADER_WINDOWS {
			let input = Input::Spans {
				file: file.try_clone()?,
				spans: vec![BlockSpan { first: 0, last }],
			};
			let mut blocks = BlockReader::new(input, workers);
			match Header::read_bam(&mut HeaderCursor::new(&mut blocks)) {
				Err(Error::Header(HeaderError::Truncated)) if last != u64::MAX => {}
				result => return result,
			}
		}
		Err(HeaderError::Truncated.into())
	}

	pub fn header(&self) -> &Header {
		&self.header
	}

	pub fn read_batch(&mut self, batch: &mut RecordBatch) -> Result<bool, Error> {
		self.batches.read_batch(batch)
	}

	pub fn next_record(&mut self) -> Result<Option<Located<'_>>, Error> {
		while self.cursor >= self.current.len() {
			if !self.batches.read_batch(&mut self.current)? {
				return Ok(None);
			}
			self.cursor = 0;
		}
		let index = self.cursor;
		self.cursor += 1;
		Ok(Some(self.current.at(index)?))
	}

	fn assemble(header: Header, batches: Pipeline) -> Self {
		Self {
			header,
			batches,
			current: RecordBatch::new(),
			cursor: 0,
		}
	}
}

impl Pipeline {
	fn read_batch(&mut self, batch: &mut RecordBatch) -> Result<bool, Error> {
		match self.ready.recv() {
			Ok(message) => {
				let used = std::mem::replace(batch, message?);
				let _ = self.recycle.try_send(used);
				Ok(true)
			}
			Err(_) => Ok(false),
		}
	}
}

impl BatchReader {
	fn new(
		blocks: BlockReader,
		leftover: Vec<Block>,
		chunks: Vec<Chunk>,
		filter: Filter,
		reference_count: u32,
	) -> Self {
		Self {
			blocks,
			leftover,
			framer: Framer::new(chunks),
			filter,
			reference_count,
		}
	}

	fn spawn(self) -> Pipeline {
		let (ready_out, ready) = sync_channel(BATCHES_READY);
		let (recycle, recycled) = sync_channel(BATCHES_RECYCLED);
		thread::spawn(move || self.run(&ready_out, &recycled));
		Pipeline { ready, recycle }
	}

	fn run(
		mut self,
		ready: &SyncSender<Result<RecordBatch, Error>>,
		recycled: &Receiver<RecordBatch>,
	) {
		loop {
			let mut batch = recycled.try_recv().unwrap_or_default();
			let message = match self.read_batch(&mut batch) {
				Ok(true) => Ok(batch),
				Ok(false) => return,
				Err(error) => Err(error),
			};
			let failed = message.is_err();
			if ready.send(message).is_err() || failed {
				return;
			}
		}
	}

	fn next_blocks(&mut self) -> Result<Option<Vec<Block>>, Error> {
		if !self.leftover.is_empty() {
			return Ok(Some(std::mem::take(&mut self.leftover)));
		}
		self.blocks.next_batch().transpose().map_err(Error::from)
	}

	fn read_batch(&mut self, batch: &mut RecordBatch) -> Result<bool, Error> {
		while let Some(blocks) = self.next_blocks()? {
			batch.reset(blocks, self.reference_count);
			self.framer.frame(batch)?;
			if let Filter::Intervals(intervals) = &self.filter {
				batch.retain(|record| intervals.overlaps(record));
			}
			if !batch.is_empty() {
				return Ok(true);
			}
		}
		self.framer.finish()?;
		Ok(false)
	}
}

fn reference_count(header: &Header) -> u32 {
	u32::try_from(header.references().len()).unwrap_or(u32::MAX)
}

struct HeaderCursor<'a> {
	reader: &'a mut BlockReader,
	loaded: Vec<Block>,
	block: usize,
	within: usize,
}

impl<'a> HeaderCursor<'a> {
	fn new(reader: &'a mut BlockReader) -> Self {
		Self {
			reader,
			loaded: Vec::new(),
			block: 0,
			within: 0,
		}
	}

	fn into_rest(mut self) -> (Vec<Block>, VirtualOffset) {
		let Some(current) = self.loaded.get(self.block) else {
			return (Vec::new(), VirtualOffset::default());
		};
		if self.within >= current.data.len() {
			let next = VirtualOffset::new(current.offset + u64::from(current.compressed_len), 0);
			return (self.loaded.split_off(self.block + 1), next);
		}
		let start = VirtualOffset::new(
			current.offset,
			u16::try_from(self.within).unwrap_or_default(),
		);
		(self.loaded.split_off(self.block), start)
	}
}

impl ByteSource for HeaderCursor<'_> {
	fn take(&mut self, len: usize, out: &mut Vec<u8>) -> Result<(), Error> {
		out.clear();
		while out.len() < len {
			let Some(current) = self.loaded.get(self.block) else {
				let batch = self.reader.next_batch().ok_or(HeaderError::Truncated)??;
				self.loaded.extend(batch);
				continue;
			};
			let available = current.data.get(self.within..).unwrap_or_default();
			if available.is_empty() {
				if self.block + 1 >= self.loaded.len() {
					let batch = self.reader.next_batch().ok_or(HeaderError::Truncated)??;
					self.loaded.extend(batch);
				}
				self.block += 1;
				self.within = 0;
				continue;
			}
			let take = available.len().min(len - out.len());
			out.extend_from_slice(available.get(..take).unwrap_or_default());
			self.within += take;
		}
		Ok(())
	}
}
