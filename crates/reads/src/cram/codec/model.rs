use super::range_coder::RangeDecoder;
use crate::cram::{Codec, Error};

const MAX_FREQUENCY: u32 = (1 << 16) - 17;
const STEP: u16 = 16;

#[derive(Debug, Clone, Copy)]
struct Entry {
	frequency: u16,
	symbol: u16,
}

#[derive(Debug, Clone)]
pub(super) struct Model {
	total: u32,
	entries: Box<[Entry]>,
}

impl Model {
	pub(super) fn new(symbols: usize) -> Self {
		let entries = (0..symbols)
			.map(|symbol| Entry {
				frequency: 1,
				symbol: u16::try_from(symbol).unwrap_or(u16::MAX),
			})
			.collect();
		Self {
			total: u32::try_from(symbols).unwrap_or(u32::MAX),
			entries,
		}
	}

	pub(super) fn decode(
		&mut self,
		decoder: &mut RangeDecoder<'_>,
		codec: Codec,
	) -> Result<u16, Error> {
		let target = decoder.frequency(self.total);
		if target > MAX_FREQUENCY {
			return Err(Error::Corrupt(codec));
		}
		let mut cumulative = 0u32;
		let mut found = None;
		for (index, entry) in self.entries.iter().enumerate() {
			cumulative += u32::from(entry.frequency);
			if cumulative > target {
				found = Some((index, *entry));
				break;
			}
		}
		let Some((index, entry)) = found else {
			return Err(Error::Corrupt(codec));
		};
		decoder.decode(
			cumulative - u32::from(entry.frequency),
			u32::from(entry.frequency),
		);
		self.update(index);
		Ok(entry.symbol)
	}

	fn update(&mut self, index: usize) {
		if let Some(entry) = self.entries.get_mut(index) {
			entry.frequency += STEP;
		}
		self.total += u32::from(STEP);
		if self.total > MAX_FREQUENCY {
			self.normalise();
		}
		if let Some(previous) = index.checked_sub(1)
			&& let Some(pair) = self.entries.get_mut(previous..=index)
			&& let [before, after] = pair
			&& after.frequency > before.frequency
		{
			std::mem::swap(before, after);
		}
	}

	fn normalise(&mut self) {
		self.total = 0;
		for entry in &mut self.entries {
			entry.frequency -= entry.frequency >> 1;
			self.total += u32::from(entry.frequency);
		}
	}
}
