use std::collections::VecDeque;

use super::kinds::Kinds;
use crate::record::Record;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Duplicates {
	#[default]
	Keep,
	Exact,
	Snps,
	Indels,
	Both,
	Any,
}

#[derive(Debug, Default)]
pub(crate) struct Window {
	queue: VecDeque<Record>,
	spare: Vec<Record>,
	span: i64,
}

impl Window {
	pub(crate) fn new(span: u32) -> Self {
		Self {
			queue: VecDeque::new(),
			spare: Vec::new(),
			span: i64::from(span),
		}
	}

	pub(crate) fn len(&self) -> usize {
		self.queue.len()
	}

	pub(crate) fn insert(&mut self, record: Record) {
		let mut at = self.queue.len();
		while let Some(previous) = at.checked_sub(1).and_then(|index| self.queue.get(index)) {
			if previous.contig != record.contig || previous.start <= record.start {
				break;
			}
			at -= 1;
		}
		self.queue.insert(at, record);
	}

	pub(crate) fn spans_positions(&self) -> bool {
		match (self.queue.front(), self.queue.back()) {
			(Some(front), Some(back)) => front.contig != back.contig || front.start != back.start,
			_ => false,
		}
	}

	pub(crate) fn due(&self) -> usize {
		let (Some(front), Some(back)) = (self.queue.front(), self.queue.back()) else {
			return 0;
		};
		self.queue
			.iter()
			.take_while(|record| {
				if front.contig == back.contig {
					record.start < back.start
						&& i64::from(back.start) - i64::from(record.start) >= self.span
				} else {
					record.contig == front.contig
				}
			})
			.count()
	}

	pub(crate) fn pop(&mut self) -> Option<Record> {
		self.queue.pop_front()
	}

	pub(crate) fn recycle(&mut self, mut record: Record) {
		record.clear();
		self.spare.push(record);
	}

	pub(crate) fn spares(&mut self) -> &mut Vec<Record> {
		&mut self.spare
	}
}

#[derive(Debug, Default)]
struct Seen {
	alleles: usize,
	reference: Vec<u8>,
	keys: Vec<u8>,
	ends: Vec<usize>,
}

impl Seen {
	fn has_key(&self, key: &[u8]) -> bool {
		let mut start = 0;
		self.ends.iter().any(|&end| {
			let found = self.keys.get(start..end) == Some(key);
			start = end;
			found
		})
	}
}

#[derive(Debug, Default)]
pub(crate) struct Dedup {
	rule: Duplicates,
	svlen_key: Option<u32>,
	group: Option<(u32, u32)>,
	kinds: Kinds,
	seen: Vec<Seen>,
	used: usize,
	key: Vec<u8>,
}

fn push_key(out: &mut Vec<u8>, allele: &[u8]) {
	if allele.first() == Some(&b'<') {
		out.extend_from_slice(allele);
		return;
	}
	let mut in_mate = false;
	for &byte in allele {
		if byte == b'[' || byte == b']' {
			in_mate = !in_mate;
			out.push(byte);
		} else if in_mate {
			out.push(byte);
		} else {
			out.push(byte.to_ascii_uppercase());
		}
	}
}

impl Dedup {
	pub(crate) fn new(rule: Duplicates, svlen_key: Option<u32>) -> Self {
		Self {
			rule,
			svlen_key,
			..Self::default()
		}
	}

	pub(crate) fn drops(&mut self, record: &Record) -> bool {
		if self.rule == Duplicates::Keep {
			return false;
		}
		let kinds = Kinds::of(record);
		let site = (record.contig, record.start);
		if self.group == Some(site) {
			if self.repeats(record, kinds) {
				return true;
			}
		} else {
			self.group = Some(site);
			self.kinds = Kinds::default();
			self.used = 0;
		}
		self.kinds = self.kinds.union(kinds);
		if self.rule == Duplicates::Exact {
			self.remember(record);
		}
		false
	}

	fn repeats(&mut self, record: &Record, kinds: Kinds) -> bool {
		let substitutions = kinds.has_substitution() && self.kinds.has_substitution();
		let indels = kinds.has_indel() && self.kinds.has_indel();
		match self.rule {
			Duplicates::Keep => false,
			Duplicates::Any => true,
			Duplicates::Snps => substitutions,
			Duplicates::Indels => indels,
			Duplicates::Both => substitutions || indels,
			Duplicates::Exact => self.exact(record),
		}
	}

	fn svlen<'a>(&self, record: &'a Record) -> &'a [i32] {
		self.svlen_key
			.and_then(|key| record.info_field(key))
			.map_or(&[][..], |field| record.store.ints(&field.cells))
	}

	fn exact(&mut self, record: &Record) -> bool {
		let reference = record.allele(0).unwrap_or_default();
		let alleles = record.allele_count();
		let svlen = self.svlen(record);
		let Some(used) = self.seen.get(..self.used) else {
			return false;
		};
		used.iter().any(|seen| {
			seen.alleles == alleles
				&& seen.reference.eq_ignore_ascii_case(reference)
				&& (1..alleles).all(|index| {
					alt_key(&mut self.key, record, index, svlen);
					if alleles == 2 {
						seen.keys == self.key
					} else {
						seen.has_key(&self.key)
					}
				})
		})
	}

	fn remember(&mut self, record: &Record) {
		if self.used == self.seen.len() {
			self.seen.push(Seen::default());
		}
		let svlen = self.svlen(record);
		let Some(seen) = self.seen.get_mut(self.used) else {
			return;
		};
		seen.alleles = record.allele_count();
		seen.reference.clear();
		seen.reference
			.extend_from_slice(record.allele(0).unwrap_or_default());
		seen.keys.clear();
		seen.ends.clear();
		for index in 1..seen.alleles {
			alt_key(&mut self.key, record, index, svlen);
			seen.keys.extend_from_slice(&self.key);
			seen.ends.push(seen.keys.len());
		}
		self.used += 1;
	}
}

fn alt_key(out: &mut Vec<u8>, record: &Record, index: usize, svlen: &[i32]) {
	out.clear();
	let Some(alt) = record.allele(index) else {
		return;
	};
	push_key(out, alt);
	if alt.first() == Some(&b'<')
		&& record.allele_count() == 2
		&& let Some(value) = svlen.first()
	{
		out.extend_from_slice(format!(".{value}").as_bytes());
	}
}
