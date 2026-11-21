use super::boxed_array;
use crate::cram::Error;

pub(super) const TOTAL_BITS: u32 = 12;
pub(super) const TOTAL: u32 = 1 << TOTAL_BITS;
pub(super) const SLOTS: usize = 1 << TOTAL_BITS;
pub(super) const SLOT_MASK: u32 = TOTAL - 1;

#[derive(Debug)]
pub(super) struct Order0Table {
	symbols: Box<[u8; SLOTS]>,
	entries: Box<[u32; SLOTS]>,
}

impl Order0Table {
	pub(super) fn new() -> Result<Self, Error> {
		Ok(Self {
			symbols: boxed_array(0)?,
			entries: boxed_array(0)?,
		})
	}

	pub(super) fn fill(&mut self, symbol: u8, start: u32, frequency: u32) -> Option<u32> {
		let end = start.checked_add(frequency).filter(|end| *end <= TOTAL)?;
		let range = start as usize..end as usize;
		let symbols = self.symbols.get_mut(range.clone())?;
		let entries = self.entries.get_mut(range)?;
		for ((slot_symbol, slot_entry), offset) in symbols.iter_mut().zip(entries).zip(0u32..) {
			*slot_symbol = symbol;
			*slot_entry = frequency | (offset << 16);
		}
		Some(end)
	}

	pub(super) fn fill_last_slot(&mut self) {
		let [.., previous, last] = &mut *self.symbols;
		*last = *previous;
		let [.., previous, last] = &mut *self.entries;
		*last = *previous + (1 << 16);
	}

	pub(super) fn advance(&self, state: u32) -> (u8, u32) {
		let slot = (state & SLOT_MASK) as usize;
		let entry = self.entries[slot];
		let next = (entry & 0xffff)
			.wrapping_mul(state >> TOTAL_BITS)
			.wrapping_add(entry >> 16);
		(self.symbols[slot], next)
	}

	pub(super) fn peek(&self, state: u32) -> u8 {
		self.symbols[(state & SLOT_MASK) as usize]
	}
}
