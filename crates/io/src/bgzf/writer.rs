use std::io::{self, Write};
use std::thread::{self, JoinHandle};

use crossbeam_channel::{Receiver, Sender, bounded, unbounded};
use libdeflater::{Compressor, Crc};
use rayon::prelude::*;

use super::block::{EMPTY_BLOCK, FOOTER_LEN, MAX_BLOCK_LEN, MAX_PAYLOAD_LEN};
use super::{BlockPosition, Level, VirtualOffset};
use crate::{Error, Workers};

const BLOCKS_PER_BATCH: usize = 256;
const BATCHES_IN_FLIGHT: usize = 4;
const HEADER: [u8; 16] = [
	0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
];
const HEADER_LEN: usize = 18;

#[derive(Debug)]
pub struct BgzfWriter<W: Write + Send + 'static> {
	payload: Vec<u8>,
	batch: Vec<Vec<u8>>,
	ordinal: u64,
	commands: Sender<Command>,
	completed: Receiver<Vec<u64>>,
	offsets: Vec<u64>,
	compressor: Option<JoinHandle<Result<Finished<W>, Error>>>,
	workers: Workers,
}

#[derive(Debug)]
enum Command {
	Blocks(Vec<Vec<u8>>),
	Finish,
}

#[derive(Debug)]
pub struct Finished<W> {
	pub inner: W,
	pub block_offsets: Vec<u64>,
}

impl<W> Finished<W> {
	pub fn resolve(&self, position: BlockPosition) -> Option<VirtualOffset> {
		let index = usize::try_from(position.ordinal).ok()?;
		Some(VirtualOffset::new(
			*self.block_offsets.get(index)?,
			position.within,
		))
	}
}

impl<W: Write + Send + 'static> BgzfWriter<W> {
	pub fn new(inner: W, level: Level, workers: &Workers) -> Self {
		let (commands, received) = bounded(BATCHES_IN_FLIGHT);
		let (reported, completed) = unbounded();
		let compressing = workers.clone();
		let compressor =
			thread::spawn(move || compress(inner, level, &received, &reported, &compressing));
		Self {
			payload: Vec::with_capacity(MAX_PAYLOAD_LEN),
			batch: Vec::with_capacity(BLOCKS_PER_BATCH),
			ordinal: 0,
			commands,
			completed,
			offsets: Vec::new(),
			compressor: Some(compressor),
			workers: workers.clone(),
		}
	}

	pub fn resolve(&mut self, position: BlockPosition) -> Option<VirtualOffset> {
		while let Ok(batch) = self.completed.try_recv() {
			self.offsets.extend(batch);
		}
		let index = usize::try_from(position.ordinal).ok()?;
		Some(VirtualOffset::new(
			*self.offsets.get(index)?,
			position.within,
		))
	}

	pub fn position(&self) -> BlockPosition {
		let within = u16::try_from(self.payload.len()).unwrap_or(u16::MAX);
		BlockPosition {
			ordinal: self.ordinal,
			within,
		}
	}

	pub fn end_block(&mut self) -> Result<(), Error> {
		if self.payload.is_empty() {
			return Ok(());
		}
		self.batch.push(std::mem::replace(
			&mut self.payload,
			Vec::with_capacity(MAX_PAYLOAD_LEN),
		));
		self.ordinal += 1;
		if self.batch.len() == BLOCKS_PER_BATCH {
			self.send_batch()?;
		}
		Ok(())
	}

	pub fn finish(mut self) -> Result<Finished<W>, Error> {
		self.end_block()?;
		self.send_batch()?;
		self.send(Command::Finish)?;
		self.join()
	}

	fn send_batch(&mut self) -> Result<(), Error> {
		if self.batch.is_empty() {
			return Ok(());
		}
		let batch = std::mem::replace(&mut self.batch, Vec::with_capacity(BLOCKS_PER_BATCH));
		self.send(Command::Blocks(batch))
	}

	fn send(&mut self, command: Command) -> Result<(), Error> {
		match self.commands.send(command) {
			Ok(()) => Ok(()),
			Err(_) => Err(self.join().err().unwrap_or(Error::Disconnected)),
		}
	}

	fn join(&mut self) -> Result<Finished<W>, Error> {
		let compressor = self.compressor.take().ok_or(Error::Disconnected)?;
		self.workers
			.join(compressor)
			.map_err(|_| Error::Disconnected)?
	}
}

impl<W: Write + Send + 'static> Drop for BgzfWriter<W> {
	fn drop(&mut self) {
		let (abandoned, _) = bounded(0);
		drop(std::mem::replace(&mut self.commands, abandoned));
		let _ = self.join();
	}
}

impl<W: Write + Send + 'static> Write for BgzfWriter<W> {
	fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
		let len = bytes.len().min(MAX_PAYLOAD_LEN - self.payload.len());
		self.payload.extend_from_slice(&bytes[..len]);
		if self.payload.len() == MAX_PAYLOAD_LEN {
			self.end_block()?;
		}
		Ok(len)
	}

	fn flush(&mut self) -> io::Result<()> {
		Ok(())
	}
}

fn compress<W: Write>(
	mut inner: W,
	level: Level,
	commands: &Receiver<Command>,
	reported: &Sender<Vec<u64>>,
	workers: &Workers,
) -> Result<Finished<W>, Error> {
	let mut block_offsets = Vec::new();
	let mut offset = 0u64;
	loop {
		let batch = match commands.recv() {
			Ok(Command::Blocks(batch)) => batch,
			Ok(Command::Finish) => break,
			Err(_) => return Err(Error::Disconnected),
		};
		let blocks = workers.install(|| {
			batch
				.par_iter()
				.map_init(
					|| Compressor::new(level.deflate()),
					|compressor, payload| encode(compressor, payload),
				)
				.collect::<Result<Vec<_>, _>>()
		})?;
		let first = block_offsets.len();
		for block in blocks {
			block_offsets.push(offset);
			inner.write_all(&block)?;
			offset += block.len() as u64;
		}
		let _ = reported.send(block_offsets.get(first..).unwrap_or_default().to_vec());
	}
	block_offsets.push(offset);
	let _ = reported.send(vec![offset]);
	inner.write_all(&EMPTY_BLOCK)?;
	inner.flush()?;
	Ok(Finished {
		inner,
		block_offsets,
	})
}

fn encode(compressor: &mut Compressor, payload: &[u8]) -> Result<Vec<u8>, Error> {
	let mut block = vec![0; MAX_BLOCK_LEN];
	let budget = MAX_BLOCK_LEN - HEADER_LEN - FOOTER_LEN;
	let compressed =
		match compressor.deflate_compress(payload, &mut block[HEADER_LEN..HEADER_LEN + budget]) {
			Ok(len) => len,
			Err(_) => Compressor::new(Level::STORED.deflate())
				.deflate_compress(payload, &mut block[HEADER_LEN..HEADER_LEN + budget])
				.map_err(|_| Error::Compression)?,
		};
	let len = HEADER_LEN + compressed + FOOTER_LEN;
	let size = u16::try_from(len - 1).map_err(|_| Error::Compression)?;
	let payload_len = u32::try_from(payload.len()).map_err(|_| Error::Compression)?;
	let mut checksum = Crc::new();
	checksum.update(payload);
	block[..16].copy_from_slice(&HEADER);
	block[16..HEADER_LEN].copy_from_slice(&size.to_le_bytes());
	let footer = HEADER_LEN + compressed;
	block[footer..footer + 4].copy_from_slice(&checksum.sum().to_le_bytes());
	block[footer + 4..len].copy_from_slice(&payload_len.to_le_bytes());
	block.truncate(len);
	Ok(block)
}
