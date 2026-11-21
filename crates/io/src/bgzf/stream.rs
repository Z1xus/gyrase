use std::io::{self, BufRead, Read};

use super::{Block, BlockReader, Input, Readahead, VirtualOffset};
use crate::Workers;

#[derive(Debug)]
pub struct BgzfReader {
	blocks: BlockReader,
	current: Block,
	position: usize,
}

impl BgzfReader {
	pub fn new(input: Input, workers: &Workers) -> Self {
		Self::with_readahead(input, workers, Readahead::DEFAULT)
	}

	pub fn with_readahead(input: Input, workers: &Workers, readahead: Readahead) -> Self {
		Self {
			blocks: BlockReader::with_readahead(input, workers, readahead),
			current: Block {
				offset: 0,
				compressed_len: 0,
				data: Vec::new(),
			},
			position: 0,
		}
	}

	pub fn virtual_offset(&self) -> VirtualOffset {
		match u16::try_from(self.position) {
			Ok(within) if self.position < self.current.data.len() => {
				VirtualOffset::new(self.current.offset, within)
			}
			_ => VirtualOffset::new(
				self.current.offset + u64::from(self.current.compressed_len),
				0,
			),
		}
	}

	fn advance(&mut self) -> io::Result<bool> {
		while self.position >= self.current.data.len() {
			match self.blocks.next_block() {
				Some(block) => {
					self.current = block?;
					self.position = 0;
				}
				None => return Ok(false),
			}
		}
		Ok(true)
	}
}

impl Read for BgzfReader {
	fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
		let available = self.fill_buf()?;
		let len = available.len().min(buffer.len());
		buffer[..len].copy_from_slice(&available[..len]);
		self.consume(len);
		Ok(len)
	}
}

impl BufRead for BgzfReader {
	fn fill_buf(&mut self) -> io::Result<&[u8]> {
		if self.advance()? {
			Ok(&self.current.data[self.position..])
		} else {
			Ok(&[])
		}
	}

	fn consume(&mut self, amount: usize) {
		self.position = (self.position + amount).min(self.current.data.len());
	}
}
