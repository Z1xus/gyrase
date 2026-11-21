use std::ops::Range;

use gyrase_genome::ContigId;

use crate::alignment::Strand;
use crate::config::Scoring;

const MAX_UNGAPPED_MISMATCHES: u32 = 8;
const DUPLICATE_OVERLAP_DIVISOR: u32 = 2;

#[derive(Debug, Clone)]
pub(crate) struct Candidate {
	pub(crate) owner: u32,
	pub(crate) strand: Strand,
	pub(crate) contig: ContigId,
	pub(crate) start: u32,
	pub(crate) stop: u32,
	pub(crate) query_start: u32,
	pub(crate) query_stop: u32,
	pub(crate) score: i32,
	pub(crate) operations: Range<u32>,
}

impl Candidate {
	pub(crate) fn reference_length(&self) -> u32 {
		self.stop - self.start
	}

	pub(crate) fn query_length(&self) -> u32 {
		self.query_stop - self.query_start
	}

	pub(crate) fn read_span(&self, read_length: u32) -> Range<u32> {
		match self.strand {
			Strand::Forward => self.query_start..self.query_stop,
			Strand::Reverse => read_length - self.query_stop..read_length - self.query_start,
		}
	}

	fn duplicates(&self, other: &Self) -> bool {
		if self.strand != other.strand || self.contig != other.contig {
			return false;
		}
		if self.start == other.start || self.stop == other.stop {
			return true;
		}
		let overlap = self
			.stop
			.min(other.stop)
			.saturating_sub(self.start.max(other.start));
		overlap * DUPLICATE_OVERLAP_DIVISOR >= self.reference_length().min(other.reference_length())
	}
}

pub(crate) fn settle(
	candidates: &mut Vec<Candidate>,
	ranges: &mut Vec<Range<usize>>,
	owners: usize,
) {
	candidates.sort_unstable_by(|left, right| {
		left.owner
			.cmp(&right.owner)
			.then(right.score.cmp(&left.score))
			.then(left.strand.index().cmp(&right.strand.index()))
			.then(left.contig.cmp(&right.contig))
			.then(left.start.cmp(&right.start))
			.then(left.stop.cmp(&right.stop))
	});
	let mut written = 0usize;
	let mut group_start = 0usize;
	let mut group_owner = None;
	for read in 0..candidates.len() {
		let Some(owner) = candidates.get(read).map(|candidate| candidate.owner) else {
			break;
		};
		if group_owner != Some(owner) {
			group_owner = Some(owner);
			group_start = written;
		}
		let duplicate = match (candidates.get(read), candidates.get(group_start..written)) {
			(Some(candidate), Some(kept)) => {
				kept.iter().any(|earlier| earlier.duplicates(candidate))
			}
			_ => true,
		};
		if !duplicate {
			candidates.swap(written, read);
			written += 1;
		}
	}
	candidates.truncate(written);
	ranges.clear();
	ranges.resize(owners, 0..0);
	let mut position = 0usize;
	while let Some(owner) = candidates.get(position).map(|candidate| candidate.owner) {
		let run = candidates
			.get(position..)
			.unwrap_or_default()
			.iter()
			.take_while(|candidate| candidate.owner == owner)
			.count();
		if let Some(slot) = ranges.get_mut(owner as usize) {
			*slot = position..position + run;
		}
		position += run;
	}
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Ungapped {
	pub(crate) score: i32,
}

pub(crate) fn ungapped(scoring: &Scoring, query: &[u8], reference: &[u8]) -> Option<Ungapped> {
	if query.len() != reference.len() {
		return None;
	}
	let matched = i32::from(scoring.matched);
	let mismatch = i32::from(scoring.mismatch);
	let clip = i32::from(scoring.clip);
	let mut mismatches = 0u32;
	let mut prefix = 0i32;
	for (&q, &r) in query.iter().zip(reference) {
		if q > 3 || r > 3 {
			return None;
		}
		if q == r {
			prefix += matched;
		} else {
			prefix -= mismatch;
			mismatches += 1;
			if mismatches > MAX_UNGAPPED_MISMATCHES {
				return None;
			}
		}
		if prefix <= -clip {
			return None;
		}
	}
	let mut suffix = 0i32;
	for (&q, &r) in query.iter().zip(reference).rev() {
		suffix += if q == r { matched } else { -mismatch };
		if suffix <= -clip {
			return None;
		}
	}
	Some(Ungapped { score: prefix })
}
