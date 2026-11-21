use std::mem::size_of;
use std::sync::atomic::AtomicBool;

use gyrase_io::Workers;
use rayon::prelude::*;

use crate::bytes::{i32_le, u16_le, u32_le};
use crate::error::{Error, SortError};
use crate::record::Flags;

const BYTES_PER_SLOT: usize = 128;
const LENGTH_PREFIX: usize = 4;
const CONTIG: usize = 0;
const POSITION: usize = 4;
const FLAGS: usize = 14;
const UNPLACED: u64 = u64::MAX;

#[derive(Debug, Clone, Copy)]
struct Slot {
	key: u64,
	offset: usize,
}

#[derive(Debug)]
pub(super) struct Buffer {
	arena: Vec<u8>,
	slots: Vec<Slot>,
	arena_limit: usize,
	slot_limit: usize,
}

pub(super) fn key(core: &[u8]) -> u64 {
	let (Some(contig), Some(position)) = (i32_le(core, CONTIG), i32_le(core, POSITION)) else {
		return UNPLACED;
	};
	let (Ok(contig), Ok(position)) = (u64::try_from(contig), u64::try_from(position)) else {
		return UNPLACED;
	};
	let reverse = u64::from(u16_le(core, FLAGS).unwrap_or(0) & Flags::REVERSE.bits() != 0);
	contig << 32 | position << 1 | reverse
}

impl Buffer {
	pub(super) fn new(budget: usize) -> Result<Self, SortError> {
		let slot_limit = budget / BYTES_PER_SLOT;
		let arena_limit = budget - slot_limit * size_of::<Slot>();
		let mut buffer = Self {
			arena: Vec::new(),
			slots: Vec::new(),
			arena_limit,
			slot_limit,
		};
		buffer
			.arena
			.try_reserve_exact(arena_limit)
			.map_err(|_| SortError::Reserve { bytes: arena_limit })?;
		buffer
			.slots
			.try_reserve_exact(slot_limit)
			.map_err(|_| SortError::Reserve {
				bytes: slot_limit * size_of::<Slot>(),
			})?;
		Ok(buffer)
	}

	pub(super) fn len(&self) -> usize {
		self.slots.len()
	}

	pub(super) fn fits(&self, record_len: usize) -> bool {
		self.slots.len() < self.slot_limit
			&& self.arena.len() + LENGTH_PREFIX + record_len <= self.arena_limit
	}

	pub(super) fn push(&mut self, core: &[u8], rest: &[u8], len: u32) {
		self.slots.push(Slot {
			key: key(core),
			offset: self.arena.len(),
		});
		self.arena.extend_from_slice(&len.to_le_bytes());
		self.arena.extend_from_slice(core);
		self.arena.extend_from_slice(rest);
	}

	pub(super) fn sort(&mut self, workers: &Workers, cancel: &AtomicBool) -> Result<(), Error> {
		super::check(cancel)?;
		let slots = &mut self.slots;
		workers.install(|| slots.par_sort_unstable_by_key(|slot| (slot.key, slot.offset)));
		super::check(cancel)
	}

	pub(super) fn clear(&mut self) {
		self.arena.clear();
		self.slots.clear();
	}

	pub(super) fn key_at(&self, index: usize) -> Option<u64> {
		Some(self.slots.get(index)?.key)
	}

	pub(super) fn record_at(&self, index: usize) -> Option<&[u8]> {
		let start = self.slots.get(index)?.offset + LENGTH_PREFIX;
		let len = u32_le(&self.arena, start - LENGTH_PREFIX)? as usize;
		self.arena.get(start..start + len)
	}
}
