use std::collections::VecDeque;
use std::fs::File;
use std::io::{self, BufRead, Read};

use gyrase_index::Chunk;
use gyrase_io::{Block, BlockReader, BlockSpan, Input, VirtualOffset, Workers};

use crate::Error;

#[derive(Debug)]
enum Stage {
	Between,
	Inside { end: VirtualOffset },
}

#[derive(Debug)]
pub struct ChunkReader {
	blocks: BlockReader,
	chunks: VecDeque<Chunk>,
	stage: Stage,
	current: Block,
	position: usize,
	limit: usize,
}

impl ChunkReader {
	pub fn new(file: File, chunks: &[Chunk], workers: &Workers) -> Self {
		let spans = chunks
			.iter()
			.map(|chunk| BlockSpan {
				first: chunk.begin.block(),
				last: chunk.end.block(),
			})
			.collect();
		Self {
			blocks: BlockReader::new(Input::Spans { file, spans }, workers),
			chunks: chunks.iter().copied().collect(),
			stage: Stage::Between,
			current: Block {
				offset: 0,
				compressed_len: 0,
				data: Vec::new(),
			},
			position: 0,
			limit: 0,
		}
	}

	fn next_block(&mut self) -> Result<(), Error> {
		self.current = self
			.blocks
			.next_block()
			.ok_or(Error::Chunk("missing block"))??;
		Ok(())
	}

	fn refill(&mut self) -> Result<bool, Error> {
		loop {
			match self.stage {
				Stage::Between => {
					let Some(chunk) = self.chunks.pop_front() else {
						return Ok(false);
					};
					self.next_block()?;
					if self.current.offset != chunk.begin.block() {
						return Err(Error::Chunk("chunk does not start at a block"));
					}
					self.position = usize::from(chunk.begin.within());
					self.stage = Stage::Inside { end: chunk.end };
				}
				Stage::Inside { end } if self.current.offset == end.block() => {
					self.stage = Stage::Between;
					continue;
				}
				Stage::Inside { .. } => {
					let Some(block) = self.blocks.next_block() else {
						self.stage = Stage::Between;
						continue;
					};
					self.current = block?;
					self.position = 0;
				}
			}
			self.limit = match self.stage {
				Stage::Inside { end } if self.current.offset == end.block() => {
					usize::from(end.within()).min(self.current.data.len())
				}
				_ => self.current.data.len(),
			};
			if self.position < self.limit {
				return Ok(true);
			}
		}
	}
}

impl BufRead for ChunkReader {
	fn fill_buf(&mut self) -> io::Result<&[u8]> {
		while self.position >= self.limit {
			if !self.refill().map_err(io::Error::other)? {
				return Ok(&[]);
			}
		}
		Ok(self
			.current
			.data
			.get(self.position..self.limit)
			.unwrap_or_default())
	}

	fn consume(&mut self, amount: usize) {
		self.position = (self.position + amount).min(self.limit);
	}
}

impl Read for ChunkReader {
	fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
		let available = self.fill_buf()?;
		let len = available.len().min(buffer.len());
		buffer
			.get_mut(..len)
			.unwrap_or_default()
			.copy_from_slice(available.get(..len).unwrap_or_default());
		self.consume(len);
		Ok(len)
	}
}
