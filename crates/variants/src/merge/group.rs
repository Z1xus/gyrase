use std::ops::Range;

use super::Collapse;
use super::input::Taken;
use super::pairing;
use crate::Error;
use crate::filter::variant_type;
use crate::record::Record;

const REF: u16 = 1;
const SNPS: u16 = 2 | 4;
const INDELS: u16 = 8;
const INSERTION: u16 = 128;
const DELETION: u16 = 256;

#[derive(Debug)]
pub(super) struct Member {
	pub(super) input: usize,
	pub(super) slot: usize,
	pub(super) map: Vec<usize>,
}

#[derive(Debug, Default)]
pub(super) struct Group {
	pub(super) members: Vec<Member>,
	pub(super) alleles: Vec<Vec<u8>>,
}

impl Group {
	pub(super) fn map(&self, input: usize) -> Option<&[usize]> {
		self.members
			.iter()
			.find(|member| member.input == input)
			.map(|member| member.map.as_slice())
	}
}

#[derive(Debug)]
struct Candidate {
	input: usize,
	slot: usize,
	types: u16,
	included: bool,
	done: bool,
	map: Vec<usize>,
}

#[derive(Debug, Default)]
pub(super) struct Groups {
	emitted: Vec<Group>,
	candidates: Vec<Candidate>,
	runs: Vec<Range<usize>>,
	reader_types: Vec<u16>,
	alleles: Vec<Vec<u8>>,
	counts: Vec<u32>,
	kinds: Vec<u16>,
	used: Vec<bool>,
}

fn kind(reference: &[u8], alt: &[u8]) -> u16 {
	match variant_type(reference, alt) {
		"snp" => 2,
		"mnp" => 4,
		"indel" if alt.len() > reference.len() => INDELS | INSERTION,
		"indel" => INDELS | DELETION,
		"other" => 16,
		"bnd" => 32,
		"overlap" => 64,
		_ => 0,
	}
}

fn types(record: &Record) -> u16 {
	let reference = record.allele(0).unwrap_or_default();
	let bits = record
		.alleles()
		.skip(1)
		.fold(0, |all, alt| all | kind(reference, alt));
	if bits == 0 { REF } else { bits }
}

fn allele_type(reference: &[u8], alt: &[u8]) -> u16 {
	let bits = kind(reference, alt);
	if bits == 0 { REF } else { bits }
}

fn same_allele(left: &[u8], right: &[u8]) -> bool {
	if left.first() == Some(&b'<')
		|| right.first() == Some(&b'<')
		|| left.contains(&b'[')
		|| left.contains(&b']')
	{
		left == right
	} else {
		left.eq_ignore_ascii_case(right)
	}
}

fn symbolic(allele: &[u8]) -> bool {
	matches!(allele.first(), Some(b'<' | b'*'))
}

struct Offset {
	tail: Vec<u8>,
	difference: isize,
}

impl Offset {
	fn new(left: &[u8], right: &[u8]) -> Option<Self> {
		let common = left
			.iter()
			.zip(right)
			.take_while(|(a, b)| a.eq_ignore_ascii_case(b))
			.count();
		if common < left.len() && common < right.len() {
			return None;
		}
		let (tail, difference) = if left.len() > right.len() {
			(
				left.get(common..),
				isize::try_from(left.len() - common).ok()?,
			)
		} else {
			(
				right.get(common..),
				-isize::try_from(right.len() - common).ok()?,
			)
		};
		Some(Self {
			tail: tail.unwrap_or_default().to_ascii_uppercase(),
			difference,
		})
	}

	fn find(&self, known: &[Vec<u8>], other: &[u8]) -> Option<usize> {
		known.iter().position(|allele| {
			let common = allele
				.iter()
				.zip(other)
				.take_while(|(a, b)| a.eq_ignore_ascii_case(b))
				.count();
			if common < allele.len() && common < other.len() {
				return false;
			}
			if self.difference == 0 {
				return common == allele.len() && common == other.len();
			}
			let rest = if common < allele.len() {
				if self.difference < 0 {
					return false;
				}
				allele.get(common..)
			} else {
				if self.difference > 0 {
					return false;
				}
				other.get(common..)
			};
			rest.is_some_and(|rest| rest.eq_ignore_ascii_case(&self.tail))
		})
	}
}

fn merge_alleles(
	incoming: &[&[u8]],
	map: &mut Vec<usize>,
	known: &mut Vec<Vec<u8>>,
) -> Result<(), Error> {
	map.clear();
	map.push(0);
	let (Some(&incoming_ref), Some(known_ref)) = (incoming.first(), known.first()) else {
		return Err(Error::Incompatible("record without alleles"));
	};
	let incoming_length = incoming_ref.len();
	let known_length = known_ref.len();
	let shared = incoming_length.min(known_length);
	let prefixes = incoming_ref
		.get(..shared)
		.zip(known_ref.get(..shared))
		.is_some_and(|(a, b)| a.eq_ignore_ascii_case(b));
	if !prefixes {
		return Err(Error::Incompatible("reference alleles differ"));
	}
	if incoming_length > known_length {
		let tail = incoming_ref.get(known_length..).unwrap_or_default();
		for allele in known.iter_mut().filter(|allele| !symbolic(allele)) {
			allele.extend_from_slice(tail);
		}
	}
	for &alt in incoming.iter().skip(1) {
		let mut expanded;
		let alt = if known_length > incoming_length && !symbolic(alt) {
			expanded = alt.to_vec();
			expanded.extend_from_slice(
				known
					.first()
					.and_then(|reference| reference.get(incoming_length..))
					.unwrap_or_default(),
			);
			expanded.as_slice()
		} else {
			alt
		};
		let found = known
			.iter()
			.skip(1)
			.position(|allele| same_allele(alt, allele));
		if let Some(position) = found {
			map.push(position + 1);
		} else {
			map.push(known.len());
			known.push(alt.to_vec());
		}
	}
	Ok(())
}

fn trailing_common(alleles: &[Vec<u8>]) -> usize {
	let Some(reference) = alleles.first() else {
		return 0;
	};
	if reference.len() < 2 {
		return 0;
	}
	let mut length = 1;
	while length < reference.len() {
		let tail = reference
			.len()
			.checked_sub(length)
			.and_then(|at| reference.get(at));
		let mut finished = false;
		for alt in alleles.iter().skip(1) {
			finished |= length >= alt.len();
			let other = alt.len().checked_sub(length).and_then(|at| alt.get(at));
			if other != tail {
				finished = true;
				break;
			}
		}
		if finished {
			break;
		}
		length += 1;
	}
	length - 1
}

impl Collapse {
	fn snps(self) -> bool {
		matches!(self, Self::Snps | Self::Both)
	}

	fn indels(self) -> bool {
		matches!(self, Self::Indels | Self::Both)
	}

	fn multiallelic(self) -> bool {
		!matches!(self, Self::None | Self::Exact)
	}
}

fn compatible(
	collapse: Collapse,
	selected: u16,
	rec_types: u16,
	record: &Record,
	known: &[Vec<u8>],
) -> bool {
	if selected == 0 {
		return true;
	}
	if collapse == Collapse::Exact && known.len() != record.allele_count() {
		return false;
	}
	if collapse == Collapse::All {
		return true;
	}
	if selected & !REF == 0 || rec_types & !REF == 0 {
		return true;
	}
	if collapse.multiallelic() {
		if collapse.snps() && rec_types & SNPS != 0 && selected & SNPS != 0 {
			return true;
		}
		if collapse.indels() && rec_types & INDELS != 0 && selected & INDELS != 0 {
			return true;
		}
	}
	let common = selected & rec_types;
	if common == 0 || (common != selected && common != rec_types) {
		return false;
	}
	if collapse == Collapse::Exact && (common != selected || common != rec_types) {
		return false;
	}
	let (Some(first), Some(reference)) = (known.first(), record.allele(0)) else {
		return false;
	};
	let Some(reference) = Offset::new(first, reference) else {
		return false;
	};
	let own_reference = record.allele(0).unwrap_or_default();
	let mut alts = record.alleles().skip(1);
	let rest = known.get(1..).unwrap_or_default();
	if collapse == Collapse::Exact {
		alts.all(|alt| reference.find(rest, alt).is_some())
	} else {
		alts.any(|alt| {
			allele_type(own_reference, alt) != REF && reference.find(rest, alt).is_some()
		})
	}
}

impl Groups {
	pub(super) fn build(
		&mut self,
		taken: &[Taken],
		ready: &[usize],
		collapse: Collapse,
	) -> Result<(), Error> {
		self.emitted.clear();
		self.candidates.clear();
		self.runs.clear();
		let firsts = pairing::first_row(taken, ready);
		for (position, &input) in ready.iter().enumerate() {
			let source = taken.get(input).ok_or(Error::Limits("inputs"))?;
			let start = self.candidates.len();
			let first = firsts.get(position).copied().flatten();
			let order = first
				.into_iter()
				.chain((0..source.count()).filter(|&slot| Some(slot) != first));
			for slot in order {
				let record = source.record(slot).ok_or(Error::Limits("records"))?;
				self.candidates.push(Candidate {
					input,
					slot,
					types: types(record),
					included: false,
					done: false,
					map: Vec::new(),
				});
			}
			self.runs.push(start..self.candidates.len());
		}
		while self.round(taken, collapse)? {}
		Ok(())
	}

	fn record<'a>(taken: &'a [Taken], candidate: &Candidate) -> Result<&'a Record, Error> {
		taken
			.get(candidate.input)
			.and_then(|source| source.record(candidate.slot))
			.ok_or(Error::Limits("records"))
	}

	fn round(&mut self, taken: &[Taken], collapse: Collapse) -> Result<bool, Error> {
		let mut remaining = 0u16;
		self.reader_types.clear();
		for run in &self.runs {
			let mut mask = 0;
			for candidate in self.candidates.get_mut(run.clone()).unwrap_or_default() {
				candidate.included = false;
				if !candidate.done {
					mask |= candidate.types;
				}
			}
			remaining |= mask;
			self.reader_types.push(mask);
		}
		if remaining == 0 {
			return Ok(false);
		}
		self.alleles.clear();
		self.counts.clear();
		self.kinds.clear();
		let mut selected = 0u16;
		for index in 0..self.candidates.len() {
			let Some(candidate) = self.candidates.get(index) else {
				break;
			};
			if candidate.done {
				continue;
			}
			let types = candidate.types;
			let record = Self::record(taken, candidate)?;
			if !compatible(collapse, selected, types, record, &self.alleles) {
				continue;
			}
			selected |= types;
			self.update_alleles(index, record)?;
		}
		self.stage(collapse, remaining)
	}

	fn update_alleles(&mut self, index: usize, record: &Record) -> Result<(), Error> {
		let reference = record.allele(0).unwrap_or_default();
		let incoming: Vec<&[u8]> = record.alleles().collect();
		let candidate = self
			.candidates
			.get_mut(index)
			.ok_or(Error::Limits("records"))?;
		candidate.included = true;
		if self.alleles.is_empty() {
			self.alleles = incoming.iter().map(|allele| allele.to_vec()).collect();
			candidate.map = (0..incoming.len()).collect();
			self.counts = vec![1; incoming.len()];
			self.kinds = incoming
				.iter()
				.enumerate()
				.map(|(at, allele)| {
					if at == 0 {
						REF
					} else {
						allele_type(reference, allele)
					}
				})
				.collect();
			return Ok(());
		}
		merge_alleles(&incoming, &mut candidate.map, &mut self.alleles)?;
		self.counts.resize(self.alleles.len(), 0);
		self.kinds.resize(self.alleles.len(), REF);
		for (at, allele) in incoming.iter().enumerate().skip(1) {
			let target = *candidate.map.get(at).ok_or(Error::Limits("alleles"))?;
			if let Some(kind) = self.kinds.get_mut(target) {
				*kind = allele_type(reference, allele);
			}
			if let Some(count) = self.counts.get_mut(target) {
				*count += 1;
			}
		}
		if let Some(count) = self.counts.first_mut() {
			*count += 1;
		}
		Ok(())
	}

	fn stage(&mut self, collapse: Collapse, remaining: u16) -> Result<bool, Error> {
		let mut frequent: Option<usize> = None;
		for at in 1..self.alleles.len() {
			if self.kinds.get(at).is_some_and(|kind| kind & REF != 0) {
				continue;
			}
			if frequent.is_none_or(|best| self.counts.get(best) < self.counts.get(at)) {
				frequent = Some(at);
			}
		}
		let wanted = frequent
			.and_then(|at| self.kinds.get(at).copied())
			.unwrap_or(REF);
		let mut chosen: Vec<usize> = Vec::new();
		for (run, &reader) in self.runs.iter().zip(&self.reader_types) {
			let usable = || {
				run.clone().filter(|&index| {
					self.candidates
						.get(index)
						.is_some_and(|candidate| candidate.included && !candidate.done)
				})
			};
			let mut found = usable().find(|&index| {
				let Some(candidate) = self.candidates.get(index) else {
					return false;
				};
				if wanted & reader != 0 && wanted & candidate.types == 0 {
					return false;
				}
				(wanted & reader == 0 && candidate.types & REF != 0)
					|| wanted == REF
					|| frequent.is_some_and(|at| candidate.map.contains(&at))
			});
			if found.is_none() && collapse.multiallelic() {
				found = usable().find(|&index| {
					self.candidates.get(index).is_some_and(|candidate| {
						Self::fallback(collapse, remaining, candidate.types)
					})
				});
			}
			chosen.extend(found);
		}
		if chosen.is_empty() {
			return Err(Error::Incompatible("records could not be merged"));
		}
		self.emit(&chosen)?;
		Ok(true)
	}

	fn fallback(collapse: Collapse, remaining: u16, line: u16) -> bool {
		if collapse == Collapse::All {
			return true;
		}
		let snps = collapse.snps();
		let indels = collapse.indels();
		if (remaining & SNPS != 0 && line & SNPS != 0 && snps)
			|| (remaining & INDELS != 0 && line & INDELS != 0 && indels)
		{
			return true;
		}
		if line & REF != 0 {
			(remaining & SNPS != 0 && snps)
				|| (remaining & INDELS != 0 && indels)
				|| remaining & REF != 0
		} else if remaining & REF != 0 {
			(line & SNPS != 0 && snps) || (line & INDELS != 0 && indels)
		} else {
			false
		}
	}

	fn emit(&mut self, chosen: &[usize]) -> Result<(), Error> {
		self.used.clear();
		self.used.resize(self.alleles.len(), false);
		for &index in chosen {
			let candidate = self.candidates.get(index).ok_or(Error::Limits("records"))?;
			for &target in candidate.map.iter().skip(1) {
				if let Some(flag) = self.used.get_mut(target) {
					*flag = true;
				}
			}
		}
		let mut renumber = vec![0usize; self.alleles.len()];
		let mut alleles: Vec<Vec<u8>> = Vec::new();
		for (at, allele) in self.alleles.iter().enumerate() {
			if at == 0 || self.used.get(at).copied().unwrap_or(false) {
				if let Some(slot) = renumber.get_mut(at) {
					*slot = alleles.len();
				}
				alleles.push(allele.clone());
			}
		}
		let trim = trailing_common(&alleles);
		for allele in &mut alleles {
			allele.truncate(allele.len().saturating_sub(trim));
		}
		let mut members = Vec::with_capacity(chosen.len());
		for &index in chosen {
			let candidate = self
				.candidates
				.get_mut(index)
				.ok_or(Error::Limits("records"))?;
			candidate.done = true;
			let map = candidate
				.map
				.iter()
				.map(|&target| renumber.get(target).copied().unwrap_or(0))
				.collect();
			members.push(Member {
				input: candidate.input,
				slot: candidate.slot,
				map,
			});
		}
		self.emitted.push(Group { members, alleles });
		Ok(())
	}

	pub(super) fn iter(&self) -> impl Iterator<Item = &Group> {
		self.emitted.iter()
	}
}
