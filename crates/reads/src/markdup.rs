mod key;

use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::{HashMap, VecDeque};
use std::hash::{BuildHasherDefault, Hasher};

use self::key::{Extent, PairKey, SingleKey, extent, mate_extent, pair_key, single_key};
use crate::error::{Error, MarkError};
use crate::header::{Header, Line};
use crate::record::{Flags, RecordRef, Value};

const FLAGS_AT: usize = 14;
const MIX: u64 = 0x517c_c1b7_2722_0a95;
const MIN_QUALITY: u8 = 15;
const EXCLUDED: Flags = Flags::new(
	Flags::SECONDARY.bits()
		| Flags::SUPPLEMENTARY.bits()
		| Flags::UNMAPPED.bits()
		| Flags::QC_FAIL.bits(),
);

#[derive(Debug, Default, Clone, Copy)]
struct Mixer(u64);

impl Hasher for Mixer {
	fn finish(&self) -> u64 {
		self.0 ^ self.0 >> 32
	}

	fn write(&mut self, bytes: &[u8]) {
		for &byte in bytes {
			self.write_u64(u64::from(byte));
		}
	}

	fn write_u8(&mut self, value: u8) {
		self.write_u64(u64::from(value));
	}

	fn write_u32(&mut self, value: u32) {
		self.write_u64(u64::from(value));
	}

	fn write_u64(&mut self, value: u64) {
		self.0 = (self.0.rotate_left(5) ^ value).wrapping_mul(MIX);
	}

	fn write_i64(&mut self, value: i64) {
		self.write_u64(value.cast_unsigned());
	}
}

type Fast<K, V> = HashMap<K, V, BuildHasherDefault<Mixer>>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MarkSettings {
	pub max_length: u32,
	pub mate_window: u32,
}

impl Default for MarkSettings {
	fn default() -> Self {
		Self {
			max_length: 300,
			mate_window: 500_000,
		}
	}
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MarkStats {
	pub reads: u64,
	pub examined: u64,
	pub excluded: u64,
	pub pairs: u64,
	pub singles: u64,
	pub pair_duplicates: u64,
	pub single_duplicates: u64,
	pub unmatched: u64,
	pub flagged: u64,
}

#[derive(Debug, Clone, Copy)]
struct Mate {
	extent: Extent,
	score: i64,
}

#[derive(Debug, Clone, Copy)]
enum MateState {
	NotNeeded,
	Waiting,
	Known(Mate),
	Missing,
}

#[derive(Debug)]
struct Item {
	bytes: Vec<u8>,
	contig: i32,
	position: i32,
	pos: i64,
	library: u32,
	examined: bool,
	has_mate: bool,
	score: i64,
	extent: Extent,
	mate: MateState,
	pair_key: Option<PairKey>,
	single_key: Option<SingleKey>,
}

#[derive(Debug)]
struct Waiting {
	seq: u64,
	mate: Mate,
}

#[derive(Debug)]
pub(crate) struct Marker {
	settings: MarkSettings,
	items: VecDeque<Item>,
	base: u64,
	decided: usize,
	groups: HashMap<Vec<u8>, u32>,
	pairs: Fast<PairKey, u64>,
	singles: Fast<SingleKey, u64>,
	waiting: Fast<Vec<u8>, Waiting>,
	latest: (i32, i32),
	previous_placed: Option<(u32, u32)>,
	stats: MarkStats,
}

impl Marker {
	pub(crate) fn new(settings: MarkSettings, header: &Header) -> Self {
		Self {
			settings,
			items: VecDeque::new(),
			base: 0,
			decided: 0,
			groups: libraries(header),
			pairs: Fast::default(),
			singles: Fast::default(),
			waiting: Fast::default(),
			latest: (0, 0),
			previous_placed: None,
			stats: MarkStats::default(),
		}
	}

	pub(crate) fn push(
		&mut self,
		bytes: &[u8],
		emit: &mut impl FnMut(&[u8]) -> Result<(), Error>,
	) -> Result<(), Error> {
		let record = RecordRef::new(bytes)?;
		if let (Some(contig), Some(position)) = (record.contig(), record.position()) {
			let placed = (contig.get(), position);
			if self
				.previous_placed
				.is_some_and(|previous| placed < previous)
			{
				return Err(MarkError::NotSorted.into());
			}
			self.previous_placed = Some(placed);
		}
		let flags = record.flags();
		let contig = record
			.contig()
			.map_or(-1, |contig| contig.get().cast_signed());
		let position = record.position().map_or(-1, u32::cast_signed);
		let mate_contig = record
			.mate_contig()
			.map_or(-1, |contig| contig.get().cast_signed());
		let mate_position = record.mate_position().map_or(-1, u32::cast_signed);
		let examined = flags.bits() & EXCLUDED.bits() == 0;
		let has_mate = flags.is_paired()
			&& !flags.is_mate_unmapped()
			&& !(mate_contig == -1 && mate_position == -1);
		let library = record
			.tags()
			.get(*b"RG")
			.ok()
			.flatten()
			.and_then(|value| self.groups.get(value.as_string()?).copied())
			.unwrap_or(0);
		let mut item = Item {
			bytes: bytes.to_vec(),
			contig,
			position,
			pos: i64::from(position),
			library,
			examined,
			has_mate,
			score: 0,
			extent: Extent::default(),
			mate: MateState::NotNeeded,
			pair_key: None,
			single_key: None,
		};
		if flags.is_duplicate() {
			set_duplicate(&mut item.bytes, false);
		}
		self.stats.reads += 1;
		if examined {
			item.score = score(record);
			item.extent = extent(record);
			if has_mate {
				item.mate = self.resolve(
					record,
					&item,
					i64::from(mate_position),
					mate_contig != contig,
				);
			}
		}
		self.latest = (contig, position);
		self.items.push_back(item);
		self.drain(false, emit)
	}

	pub(crate) fn finish(
		mut self,
		emit: &mut impl FnMut(&[u8]) -> Result<(), Error>,
	) -> Result<MarkStats, Error> {
		self.drain(true, emit)?;
		for item in self.items.drain(..) {
			self.stats.flagged += u64::from(is_duplicate(&item.bytes));
			emit(&item.bytes)?;
		}
		Ok(self.stats)
	}

	fn resolve(
		&mut self,
		record: RecordRef<'_>,
		item: &Item,
		mate_position: i64,
		other_contig: bool,
	) -> MateState {
		let tags = record.tags();
		let text = tags.get(*b"MC").ok().flatten();
		let mate_score = tags.get(*b"ms").ok().flatten();
		if let (Some(text), Some(mate_score)) = (
			text.as_ref().and_then(|value| value.as_string()),
			mate_score.as_ref().and_then(Value::as_integer),
		) {
			return MateState::Known(Mate {
				extent: mate_extent(mate_position, text),
				score: mate_score,
			});
		}
		if other_contig {
			self.stats.unmatched += 1;
			return MateState::Missing;
		}
		let own = Mate {
			extent: item.extent,
			score: item.score,
		};
		let seq = self.base + self.items.len() as u64;
		if let Some(waiting) = self.waiting.remove(record.name()) {
			if let Some(other) = item_mut(&mut self.items, self.base, waiting.seq) {
				other.mate = MateState::Known(own);
			}
			MateState::Known(waiting.mate)
		} else {
			self.waiting
				.insert(record.name().to_vec(), Waiting { seq, mate: own });
			MateState::Waiting
		}
	}

	fn drain(
		&mut self,
		last: bool,
		emit: &mut impl FnMut(&[u8]) -> Result<(), Error>,
	) -> Result<(), Error> {
		while let Some(item) = self.items.get(self.decided) {
			if matches!(item.mate, MateState::Waiting) {
				let window = i64::from(self.settings.mate_window);
				let expired = last
					|| self.latest.0 != item.contig
					|| i64::from(self.latest.1) > i64::from(item.position) + window;
				if !expired {
					break;
				}
				let name = RecordRef::new(&item.bytes)?.name().to_vec();
				self.waiting.remove(&name);
				self.stats.unmatched += 1;
				if let Some(item) = self.items.get_mut(self.decided) {
					item.mate = MateState::Missing;
				}
			}
			self.decide(self.decided)?;
			let (contig, position) = self
				.items
				.get(self.decided)
				.map_or((0, 0), |item| (item.contig, item.position));
			self.decided += 1;
			self.evict((contig, position), emit)?;
		}
		Ok(())
	}

	fn evict(
		&mut self,
		current: (i32, i32),
		emit: &mut impl FnMut(&[u8]) -> Result<(), Error>,
	) -> Result<(), Error> {
		let max_length = i64::from(self.settings.max_length);
		while self.decided > 0 {
			let Some(front) = self.items.front() else {
				break;
			};
			let in_window = front.pos + max_length > i64::from(current.1)
				&& front.contig == current.0
				&& (current.0 != -1 || current.1 != -1);
			if in_window {
				break;
			}
			let Some(item) = self.items.pop_front() else {
				break;
			};
			if let Some(key) = item.pair_key {
				self.pairs.remove(&key);
			}
			if let Some(key) = item.single_key {
				self.singles.remove(&key);
			}
			self.decided -= 1;
			self.base += 1;
			self.stats.flagged += u64::from(is_duplicate(&item.bytes));
			emit(&item.bytes)?;
		}
		Ok(())
	}

	fn decide(&mut self, index: usize) -> Result<(), Error> {
		let Some(item) = self.items.get(index) else {
			return Ok(());
		};
		if !item.examined {
			self.stats.excluded += 1;
			return Ok(());
		}
		self.stats.examined += 1;
		let seq = self.base + index as u64;
		let record = RecordRef::new(&item.bytes)?;
		let flags = record.flags();
		let library = item.library;
		let contig = i64::from(item.contig) + 1;
		let single = single_key(library, contig, flags, item.extent);
		let mate_contig = record
			.mate_contig()
			.map_or(-1, |contig| i64::from(contig.get().cast_signed()))
			+ 1;
		let paired = match item.mate {
			MateState::Known(mate) => Some((
				pair_key(
					library,
					contig,
					mate_contig,
					flags,
					item.extent,
					mate.extent,
				),
				item.score + mate.score,
			)),
			_ => None,
		};
		let own_score = item.score;
		if let Some(item) = self.items.get_mut(index) {
			item.pos = single.coord();
		}
		if let Some((pair, pair_score)) = paired {
			self.stats.pairs += 1;
			self.decide_pair_single(seq, single);
			self.decide_pair(seq, pair, pair_score);
		} else {
			self.stats.singles += 1;
			self.decide_single(seq, single, own_score);
		}
		Ok(())
	}

	fn decide_pair_single(&mut self, seq: u64, key: SingleKey) {
		match self.singles.entry(key) {
			Entry::Vacant(slot) => {
				slot.insert(seq);
				if let Some(item) = item_mut(&mut self.items, self.base, seq) {
					item.single_key = Some(key);
				}
			}
			Entry::Occupied(mut slot) => {
				let other = *slot.get();
				if item_ref(&self.items, self.base, other).is_some_and(|item| !item.has_mate) {
					slot.insert(seq);
					mark(&mut self.items, self.base, other);
					self.stats.single_duplicates += 1;
				}
			}
		}
	}

	fn decide_pair(&mut self, seq: u64, key: PairKey, score: i64) {
		match self.pairs.entry(key) {
			Entry::Vacant(slot) => {
				slot.insert(seq);
				if let Some(item) = item_mut(&mut self.items, self.base, seq) {
					item.pair_key = Some(key);
				}
			}
			Entry::Occupied(mut slot) => {
				let other = *slot.get();
				let Some(old) = item_ref(&self.items, self.base, other) else {
					return;
				};
				let old_score = old.score
					+ match old.mate {
						MateState::Known(mate) => mate.score,
						_ => 0,
					};
				let new_wins = match score.cmp(&old_score) {
					Ordering::Greater => true,
					Ordering::Less => false,
					Ordering::Equal => {
						let name = item_ref(&self.items, self.base, seq)
							.and_then(|item| RecordRef::new(&item.bytes).ok())
							.map(|record| record.name());
						let old = RecordRef::new(&old.bytes).ok().map(|record| record.name());
						name.zip(old).is_some_and(|(name, old)| name < old)
					}
				};
				let duplicate = if new_wins {
					slot.insert(seq);
					other
				} else {
					seq
				};
				mark(&mut self.items, self.base, duplicate);
				self.stats.pair_duplicates += 1;
			}
		}
	}

	fn decide_single(&mut self, seq: u64, key: SingleKey, score: i64) {
		match self.singles.entry(key) {
			Entry::Vacant(slot) => {
				slot.insert(seq);
				if let Some(item) = item_mut(&mut self.items, self.base, seq) {
					item.single_key = Some(key);
				}
			}
			Entry::Occupied(mut slot) => {
				let other = *slot.get();
				let duplicate = match item_ref(&self.items, self.base, other) {
					Some(old) if !old.has_mate && score > old.score => {
						slot.insert(seq);
						other
					}
					_ => seq,
				};
				mark(&mut self.items, self.base, duplicate);
				self.stats.single_duplicates += 1;
			}
		}
	}
}

fn item_ref(items: &VecDeque<Item>, base: u64, seq: u64) -> Option<&Item> {
	items.get(usize::try_from(seq.checked_sub(base)?).ok()?)
}

fn item_mut(items: &mut VecDeque<Item>, base: u64, seq: u64) -> Option<&mut Item> {
	items.get_mut(usize::try_from(seq.checked_sub(base)?).ok()?)
}

fn mark(items: &mut VecDeque<Item>, base: u64, seq: u64) {
	if let Some(item) = item_mut(items, base, seq) {
		set_duplicate(&mut item.bytes, true);
	}
}

fn score(record: RecordRef<'_>) -> i64 {
	record
		.qualities()
		.as_phred()
		.iter()
		.filter(|&&quality| quality >= MIN_QUALITY)
		.map(|&quality| i64::from(quality))
		.sum()
}

fn libraries(header: &Header) -> HashMap<Vec<u8>, u32> {
	let mut names: HashMap<&str, u32> = HashMap::new();
	let mut groups = HashMap::new();
	for line in header.lines() {
		if let Line::ReadGroup(fields) = line {
			let library = match fields.get(*b"LB").filter(|name| !name.is_empty()) {
				Some(name) => {
					let next = u32::try_from(names.len() + 1).unwrap_or(u32::MAX);
					*names.entry(name).or_insert(next)
				}
				None => 0,
			};
			groups.insert(
				fields.get(*b"ID").unwrap_or_default().as_bytes().to_vec(),
				library,
			);
		}
	}
	groups
}

pub(crate) fn is_duplicate(bytes: &[u8]) -> bool {
	bytes
		.get(FLAGS_AT..FLAGS_AT + 2)
		.is_some_and(|flags| Flags::new(u16::from_le_bytes([flags[0], flags[1]])).is_duplicate())
}

pub(crate) fn set_duplicate(bytes: &mut [u8], on: bool) {
	if let Some(flags) = bytes.get_mut(FLAGS_AT..FLAGS_AT + 2) {
		let current = u16::from_le_bytes([flags[0], flags[1]]);
		let updated = if on {
			current | Flags::DUPLICATE.bits()
		} else {
			current & !Flags::DUPLICATE.bits()
		};
		flags.copy_from_slice(&updated.to_le_bytes());
	}
}
