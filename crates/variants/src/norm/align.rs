use gyrase_genome::{ContigId, Interval, Reference};

use super::kinds::Kinds;
use crate::Error;
use crate::record::Record;
use crate::value::Cells;

const FIRST_FETCH: u64 = 100;
const LARGEST_FETCH: u64 = 65536;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Outcome {
	Unchanged,
	Changed,
	Mismatch,
	Duplicate,
}

pub(crate) struct Place<'a> {
	pub reference: &'a Reference,
	pub contig: ContigId,
	pub end_key: Option<u32>,
}

#[derive(Debug, Default)]
pub(crate) struct Aligner {
	alleles: Vec<Vec<u8>>,
	count: usize,
	fetched: Vec<u8>,
}

const fn iupac_first(byte: u8) -> Option<u8> {
	Some(match byte {
		b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n' => byte,
		b'R' | b'W' | b'M' | b'D' | b'H' | b'V' => b'A',
		b'Y' | b'S' | b'B' => b'C',
		b'K' => b'G',
		b'r' | b'w' | b'm' | b'd' | b'h' | b'v' => b'a',
		b'y' | b's' | b'b' => b'c',
		b'k' => b'g',
		_ => return None,
	})
}

fn has_non_acgtn(allele: &[u8]) -> bool {
	allele.iter().any(|&byte| iupac_first(byte) != Some(byte))
}

fn unambiguous(bases: &mut [u8]) {
	for base in bases {
		if let Some(first) = iupac_first(*base) {
			*base = first;
		}
	}
}

fn same_base(left: Option<&u8>, right: Option<&u8>) -> bool {
	matches!((left, right), (Some(left), Some(right)) if left.eq_ignore_ascii_case(right))
}

impl Aligner {
	fn alleles(&self) -> &[Vec<u8>] {
		self.alleles.get(..self.count).unwrap_or_default()
	}

	fn alleles_mut(&mut self) -> &mut [Vec<u8>] {
		self.alleles.get_mut(..self.count).unwrap_or_default()
	}

	fn load(&mut self, record: &Record) {
		self.count = record.allele_count();
		if self.alleles.len() < self.count {
			self.alleles.resize_with(self.count, Vec::new);
		}
		for (slot, allele) in self.alleles.iter_mut().zip(record.alleles()) {
			slot.clear();
			slot.extend_from_slice(allele);
		}
	}

	fn fetch(&mut self, place: &Place<'_>, start: u64, end: u64) -> Result<(), Error> {
		let interval = Interval {
			contig: place.contig,
			start,
			end,
		};
		place.reference.fetch_cased(&interval, &mut self.fetched)?;
		unambiguous(&mut self.fetched);
		Ok(())
	}

	fn contig_length(place: &Place<'_>) -> u64 {
		place
			.reference
			.contigs()
			.get(place.contig)
			.map_or(0, |contig| contig.length)
	}

	pub(crate) fn realign(
		&mut self,
		place: &Place<'_>,
		record: &mut Record,
	) -> Result<Outcome, Error> {
		let start = u64::from(record.start);
		let reference_allele = record.allele(0).unwrap_or_default();
		let reference_len = reference_allele.len() as u64;
		let stop = (start + reference_len).min(Self::contig_length(place));
		self.fetch(place, start.min(stop), stop)?;
		if has_non_acgtn(reference_allele) || !self.fetched.eq_ignore_ascii_case(reference_allele) {
			return Ok(Outcome::Mismatch);
		}
		if record.allele_count() == 1 {
			if reference_len > 1 {
				record.truncate_reference(1);
			}
			return Ok(Outcome::Unchanged);
		}
		let kinds = Kinds::of(record);
		if kinds.has_breakend() || kinds.has_overlap() {
			return Ok(Outcome::Unchanged);
		}
		if record
			.alleles()
			.skip(1)
			.any(|allele| allele.first() == Some(&b'<'))
		{
			return Ok(Outcome::Unchanged);
		}
		self.load(record);
		let reference = self.alleles().first().map_or(&[][..], Vec::as_slice);
		for allele in self.alleles().iter().skip(1) {
			if has_non_acgtn(allele) {
				return Ok(Outcome::Mismatch);
			}
			if allele.eq_ignore_ascii_case(reference) {
				return Ok(Outcome::Duplicate);
			}
		}
		let new_start = self.left_align(place, start)?;
		let unchanged = new_start == start
			&& self
				.alleles()
				.first()
				.is_some_and(|allele| allele.eq_ignore_ascii_case(reference_allele));
		if unchanged {
			return Ok(Outcome::Unchanged);
		}
		let new_reference_len = self.alleles().first().map_or(0, Vec::len);
		record.set_alleles(self.alleles().iter().map(Vec::as_slice));
		record.set_start(u32::try_from(new_start).map_err(|_| Error::Limits("position"))?);
		if let Some(key) = place.end_key {
			let moved = new_start != start || reference_len != new_reference_len as u64;
			let single = record
				.info_field(key)
				.is_some_and(|field| matches!(field.cells, Cells::Ints(_)) && field.width == 1);
			if moved && single {
				let end = i32::try_from(new_start + new_reference_len as u64)
					.map_err(|_| Error::Limits("position"))?;
				record.set_info_ints(key, &[end]);
			}
		}
		Ok(Outcome::Changed)
	}

	fn left_align(&mut self, place: &Place<'_>, start: u64) -> Result<u64, Error> {
		loop {
			let alleles = self.alleles();
			let last = alleles.first().and_then(|allele| allele.last());
			let agree = alleles.iter().all(|allele| same_base(allele.last(), last));
			let shortest = alleles.iter().map(Vec::len).min().unwrap_or(0);
			if !agree || (shortest <= 1 && start == 0) {
				break;
			}
			let mut emptied = false;
			for allele in self.alleles_mut() {
				allele.pop();
				emptied |= allele.is_empty();
			}
			if emptied {
				return self.shift_left(place, start);
			}
		}
		let mut trimmed = 0;
		loop {
			let alleles = self.alleles();
			let first = alleles.first().and_then(|allele| allele.get(trimmed));
			let agree = alleles
				.iter()
				.all(|allele| same_base(allele.get(trimmed), first));
			let shortest = alleles
				.iter()
				.map(|allele| allele.len().saturating_sub(trimmed))
				.min()
				.unwrap_or(0);
			if !agree || shortest <= 1 {
				break;
			}
			trimmed += 1;
		}
		if trimmed > 0 {
			for allele in self.alleles_mut() {
				allele.drain(..trimmed);
			}
		}
		Ok(start + trimmed as u64)
	}

	fn shift_left(&mut self, place: &Place<'_>, start: u64) -> Result<u64, Error> {
		let longest = self.alleles().iter().map(Vec::len).max().unwrap_or(0) as u64;
		let mut block = FIRST_FETCH;
		let mut shift = 0u64;
		while start.saturating_sub(shift) > 1 {
			let last = start - shift - 1;
			let first = (last + 1).saturating_sub(block).max(1);
			self.fetch(place, first, last + 1)?;
			let mut remaining = self.fetched.len();
			while remaining > 0 {
				let base = self.fetched.get(remaining - 1).map(u8::to_ascii_uppercase);
				let movable = self.alleles().iter().all(|allele| {
					let length = allele.len() as u64;
					length == 0
						|| (length - 1)
							.checked_sub(shift % length)
							.and_then(|index| usize::try_from(index).ok())
							.and_then(|index| allele.get(index))
							.map(u8::to_ascii_uppercase)
							== base
				});
				if !movable {
					break;
				}
				remaining -= 1;
				shift += 1;
			}
			if remaining > 0 {
				break;
			}
			if block < LARGEST_FETCH {
				block = if block > LARGEST_FETCH / 2 {
					LARGEST_FETCH
				} else {
					2 * block
				};
			}
		}
		let new_start = start - shift - 1;
		let reach = shift.min(longest);
		self.fetch(place, new_start, new_start + reach + 1)?;
		let prefix = std::mem::take(&mut self.fetched);
		for allele in self.alleles_mut() {
			let length = allele.len() as u64;
			let kept = usize::try_from(length - shift.min(length)).unwrap_or(0);
			let used = usize::try_from(shift.min(length)).unwrap_or(0);
			allele.truncate(kept);
			drop(allele.splice(
				0..0,
				prefix.get(..=used).unwrap_or_default().iter().copied(),
			));
		}
		self.fetched = prefix;
		Ok(new_start)
	}
}
