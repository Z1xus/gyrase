use std::mem::size_of;
use std::sync::atomic::AtomicBool;

use gyrase_io::Workers;
use rayon::prelude::*;

use crate::bytes::{i32_le, u16_le, u32_le};
use crate::error::{Error, SortError};
use crate::record::{Flags, RecordRef};

const BYTES_PER_SLOT: usize = 128;
const LENGTH_PREFIX: usize = 4;
const ARENA: usize = 4 << 20;
const CONTIG: usize = 0;
const POSITION: usize = 4;
const FLAGS: usize = 14;
const UNPLACED: u64 = u64::MAX;

#[derive(Debug, Clone, Copy)]
struct Slot {
	key: u64,
	arena: u32,
	offset: u32,
}

#[derive(Debug, Default)]
pub struct Encoded {
	bytes: Vec<u8>,
	keys: Vec<(u64, u32)>,
	longest: u64,
}

#[derive(Debug)]
pub(super) struct Buffer {
	arenas: Vec<Vec<u8>>,
	slots: Vec<Slot>,
	bytes: usize,
	arena_limit: usize,
	slot_limit: usize,
	longest: u64,
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

fn arena_index(index: usize) -> u32 {
	u32::try_from(index).unwrap_or(u32::MAX)
}

fn span(record: RecordRef<'_>) -> u64 {
	let start = u64::from(record.position().unwrap_or(0));
	record
		.reference_end()
		.map_or(0, |end| end.saturating_sub(start))
}

impl Encoded {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn take(&mut self) -> Self {
		let (bytes, keys) = (self.bytes.len(), self.keys.len());
		Self {
			bytes: std::mem::replace(&mut self.bytes, Vec::with_capacity(bytes)),
			keys: std::mem::replace(&mut self.keys, Vec::with_capacity(keys)),
			longest: std::mem::take(&mut self.longest),
		}
	}

	pub fn bytes(&self) -> usize {
		self.bytes.len()
	}

	pub fn records(&self) -> usize {
		self.keys.len()
	}

	pub fn push(&mut self, record: RecordRef<'_>) -> Result<(), SortError> {
		let (core, rest) = record.parts();
		let len = u32::try_from(core.len() + rest.len()).map_err(|_| SortError::RecordTooLarge)?;
		let offset = u32::try_from(self.bytes.len()).map_err(|_| SortError::RecordTooLarge)?;
		self.keys.push((key(core), offset));
		self.bytes.extend_from_slice(&len.to_le_bytes());
		self.bytes.extend_from_slice(core);
		self.bytes.extend_from_slice(rest);
		self.longest = self.longest.max(span(record));
		Ok(())
	}
}

impl Buffer {
	pub(super) fn new(budget: usize) -> Result<Self, SortError> {
		let slot_limit = budget / BYTES_PER_SLOT;
		let arena_limit = budget - slot_limit * size_of::<Slot>();
		let mut buffer = Self {
			arenas: Vec::new(),
			slots: Vec::new(),
			bytes: 0,
			arena_limit,
			slot_limit,
			longest: 0,
		};
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

	pub(super) fn longest(&self) -> u64 {
		self.longest
	}

	pub(super) fn fits(&self, bytes: usize, records: usize) -> bool {
		self.slots.len() + records <= self.slot_limit
			&& self.bytes + records * LENGTH_PREFIX + bytes <= self.arena_limit
	}

	pub(super) fn push(&mut self, record: RecordRef<'_>, len: u32) {
		let (core, rest) = record.parts();
		let (arena, tail) = self.tail(LENGTH_PREFIX + len as usize);
		let offset = u32::try_from(tail.len()).unwrap_or(u32::MAX);
		tail.extend_from_slice(&len.to_le_bytes());
		tail.extend_from_slice(core);
		tail.extend_from_slice(rest);
		self.slots.push(Slot {
			key: key(core),
			arena,
			offset,
		});
		self.bytes += LENGTH_PREFIX + len as usize;
		self.longest = self.longest.max(span(record));
	}

	pub(super) fn append(&mut self, chunk: &Encoded) {
		let (arena, tail) = self.tail(chunk.bytes.len());
		let base = u32::try_from(tail.len()).unwrap_or(u32::MAX);
		tail.extend_from_slice(&chunk.bytes);
		self.bytes += chunk.bytes.len();
		self.longest = self.longest.max(chunk.longest);
		self.slots
			.extend(chunk.keys.iter().map(|&(key, offset)| Slot {
				key,
				arena,
				offset: base + offset,
			}));
	}

	fn tail(&mut self, needed: usize) -> (u32, &mut Vec<u8>) {
		if self
			.arenas
			.last()
			.is_none_or(|tail| tail.capacity() - tail.len() < needed)
		{
			self.arenas.push(Vec::with_capacity(ARENA.max(needed)));
		}
		let index = self.arenas.len() - 1;
		(arena_index(index), &mut self.arenas[index])
	}

	pub(super) fn sort(&mut self, workers: &Workers, cancel: &AtomicBool) -> Result<(), Error> {
		super::check(cancel)?;
		let slots = &mut self.slots;
		workers
			.install(|| slots.par_sort_unstable_by_key(|slot| (slot.key, slot.arena, slot.offset)));
		super::check(cancel)
	}

	pub(super) fn clear(&mut self) {
		self.arenas.clear();
		self.slots.clear();
		self.bytes = 0;
	}

	pub(super) fn lower_bound(&self, key: u64) -> usize {
		self.slots.partition_point(|slot| slot.key < key)
	}

	pub(super) fn key_at(&self, index: usize) -> Option<u64> {
		Some(self.slots.get(index)?.key)
	}

	pub(super) fn record_at(&self, index: usize) -> Option<&[u8]> {
		let slot = self.slots.get(index)?;
		let arena = self.arenas.get(slot.arena as usize)?;
		let len = u32_le(arena, slot.offset as usize)? as usize;
		let start = slot.offset as usize + LENGTH_PREFIX;
		arena.get(start..start + len)
	}

	pub(super) fn record_mut(&mut self, index: usize) -> Option<&mut [u8]> {
		let slot = *self.slots.get(index)?;
		let arena = self.arenas.get_mut(slot.arena as usize)?;
		let len = u32_le(arena, slot.offset as usize)? as usize;
		let start = slot.offset as usize + LENGTH_PREFIX;
		arena.get_mut(start..start + len)
	}
}
