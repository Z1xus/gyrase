use std::hint::black_box;

use gyrase_io::{Block, Chunk, VirtualOffset};

use super::batch::{Entry, Location, RecordBatch};
use crate::bytes::u32_le;
use crate::error::{BamError, Error};

const PREFIX_LEN: usize = 4;
const MIN_RECORD_LEN: usize = 32;
const MAX_RECORD_LEN: usize = 1 << 26;
const READ_AHEAD: usize = 4096;

#[derive(Debug, Clone, Copy)]
enum Walk {
	Waiting { chunk: usize },
	Inside { chunk: usize },
}

#[derive(Debug, Clone, Copy)]
enum Carry {
	Empty,
	Partial { start: VirtualOffset },
}

#[derive(Debug)]
pub(super) struct Framer {
	chunks: Vec<Chunk>,
	walk: Walk,
	carry: Carry,
	scratch: Vec<u8>,
}

struct Output<'a> {
	stitched: &'a mut Vec<u8>,
	entries: &'a mut Vec<Entry>,
}

impl Framer {
	pub(super) fn new(chunks: Vec<Chunk>) -> Self {
		Self {
			chunks,
			walk: Walk::Waiting { chunk: 0 },
			carry: Carry::Empty,
			scratch: Vec::new(),
		}
	}

	pub(super) fn frame(&mut self, batch: &mut RecordBatch) -> Result<(), Error> {
		let (blocks, stitched, entries) = batch.parts_mut();
		let mut output = Output { stitched, entries };
		for (index, block) in blocks.iter().enumerate() {
			let index = u32::try_from(index).map_err(|_| BamError::BatchTooLarge)?;
			self.frame_block(index, block, &mut output)?;
		}
		Ok(())
	}

	pub(super) fn finish(&self) -> Result<(), BamError> {
		match self.carry {
			Carry::Empty => Ok(()),
			Carry::Partial { .. } => Err(BamError::TruncatedRecord),
		}
	}

	fn enter(&mut self, block: &Block) -> Result<(usize, Chunk, usize), BamError> {
		match self.walk {
			Walk::Waiting { chunk } => {
				let current = *self.chunks.get(chunk).ok_or(BamError::UnexpectedBlock)?;
				let within = usize::from(current.begin.within());
				if block.offset != current.begin.block() || within > block.data.len() {
					return Err(BamError::MisalignedChunk);
				}
				self.walk = Walk::Inside { chunk };
				Ok((chunk, current, within))
			}
			Walk::Inside { chunk } => {
				let current = *self.chunks.get(chunk).ok_or(BamError::UnexpectedBlock)?;
				Ok((chunk, current, 0))
			}
		}
	}

	fn frame_block(
		&mut self,
		index: u32,
		block: &Block,
		output: &mut Output<'_>,
	) -> Result<(), Error> {
		let (chunk_index, chunk, first) = self.enter(block)?;
		let data = block.data.as_slice();
		let mut position = first;
		if let Carry::Partial { start } = self.carry {
			position = self.extend_carry(data, position, start, block, output)?;
		}
		while position < data.len() && matches!(self.carry, Carry::Empty) {
			let start = VirtualOffset::new(block.offset, offset_within(position)?);
			if start >= chunk.end {
				break;
			}
			let remaining = data.get(position..).unwrap_or_default();
			black_box(remaining.get(READ_AHEAD).copied());
			match record_len(remaining)? {
				Some(len) if remaining.len() >= PREFIX_LEN + len => {
					let end_position = position + PREFIX_LEN + len;
					output.entries.push(Entry {
						location: Location::Block {
							index,
							offset: to_u32(position + PREFIX_LEN)?,
							len: to_u32(len)?,
						},
						start,
						end: end_offset(block, end_position),
					});
					position = end_position;
				}
				_ => {
					self.scratch.clear();
					self.scratch.extend_from_slice(remaining);
					self.carry = Carry::Partial { start };
					position = data.len();
				}
			}
		}
		if block.offset + u64::from(block.compressed_len) > chunk.end.block() {
			self.finish()?;
			self.walk = Walk::Waiting {
				chunk: chunk_index + 1,
			};
		}
		Ok(())
	}

	fn extend_carry(
		&mut self,
		data: &[u8],
		position: usize,
		start: VirtualOffset,
		block: &Block,
		output: &mut Output<'_>,
	) -> Result<usize, Error> {
		let mut position = position;
		if self.scratch.len() < PREFIX_LEN {
			position += self.append(data, position, PREFIX_LEN - self.scratch.len());
			if self.scratch.len() < PREFIX_LEN {
				return Ok(position);
			}
		}
		let len = record_len(&self.scratch)?.ok_or(BamError::TruncatedRecord)?;
		let total = PREFIX_LEN + len;
		position += self.append(data, position, total - self.scratch.len());
		if self.scratch.len() == total {
			let body = self.scratch.get(PREFIX_LEN..).unwrap_or_default();
			output.entries.push(Entry {
				location: Location::Stitched {
					offset: to_u32(output.stitched.len())?,
					len: to_u32(len)?,
				},
				start,
				end: end_offset(block, position),
			});
			output.stitched.extend_from_slice(body);
			self.carry = Carry::Empty;
		}
		Ok(position)
	}

	fn append(&mut self, data: &[u8], position: usize, wanted: usize) -> usize {
		let available = data.get(position..).unwrap_or_default();
		let take = available.len().min(wanted);
		self.scratch
			.extend_from_slice(available.get(..take).unwrap_or_default());
		take
	}
}

fn record_len(bytes: &[u8]) -> Result<Option<usize>, BamError> {
	let Some(raw) = u32_le(bytes, 0) else {
		return Ok(None);
	};
	usize::try_from(raw)
		.ok()
		.filter(|len| (MIN_RECORD_LEN..=MAX_RECORD_LEN).contains(len))
		.map(Some)
		.ok_or(BamError::RecordSize {
			len: u64::from(raw),
		})
}

fn offset_within(position: usize) -> Result<u16, BamError> {
	u16::try_from(position).map_err(|_| BamError::UnexpectedBlock)
}

fn to_u32(value: usize) -> Result<u32, BamError> {
	u32::try_from(value).map_err(|_| BamError::BatchTooLarge)
}

fn end_offset(block: &Block, position: usize) -> VirtualOffset {
	match u16::try_from(position) {
		Ok(within) if position < block.data.len() => VirtualOffset::new(block.offset, within),
		_ => VirtualOffset::new(block.offset + u64::from(block.compressed_len), 0),
	}
}
