use super::stream::Stream;
use crate::cram::Error;
use crate::cram::cursor::Cursor;

const TOP: u32 = 1 << 24;
const START_BYTES: usize = 5;

#[derive(Debug)]
pub(super) struct RangeDecoder<'a> {
	range: u32,
	code: u32,
	stream: Stream<'a>,
}

impl<'a> RangeDecoder<'a> {
	pub(super) fn new(data: &'a [u8]) -> Result<Self, Error> {
		let mut cursor = Cursor::new(data);
		let mut code = 0u32;
		for byte in cursor.take(START_BYTES)? {
			code = (code << 8) | u32::from(*byte);
		}
		Ok(Self {
			range: u32::MAX,
			code,
			stream: Stream::new(cursor.rest()),
		})
	}

	pub(super) fn frequency(&mut self, total: u32) -> u32 {
		if total != 0 && self.range >= total {
			self.range /= total;
			self.code.checked_div(self.range).unwrap_or(0)
		} else {
			0
		}
	}

	pub(super) fn decode(&mut self, cumulative: u32, frequency: u32) {
		self.code = self.code.wrapping_sub(cumulative.wrapping_mul(self.range));
		self.range = self.range.wrapping_mul(frequency);
		while self.range < TOP {
			self.code = (self.code << 8) | self.stream.byte();
			self.range <<= 8;
			if self.range == 0 {
				break;
			}
		}
	}

	pub(super) fn is_overrun(&self) -> bool {
		self.stream.is_overrun()
	}
}
