use gyrase_genome::{ContigId, Interval};
use gyrase_io::VirtualOffset;

use crate::record::RecordRef;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Chunk {
	pub start: VirtualOffset,
	pub end: VirtualOffset,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Query {
	chunks: Vec<Chunk>,
	intervals: IntervalSet,
}

impl Query {
	pub fn new(chunks: Vec<Chunk>, intervals: Vec<Interval>) -> Self {
		Self {
			chunks: merge_chunks(chunks),
			intervals: IntervalSet::new(intervals),
		}
	}

	pub(super) fn into_parts(self) -> (Vec<Chunk>, IntervalSet) {
		(self.chunks, self.intervals)
	}
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct IntervalSet {
	intervals: Vec<Interval>,
}

impl IntervalSet {
	fn new(mut intervals: Vec<Interval>) -> Self {
		intervals.retain(|interval| !interval.is_empty());
		intervals.sort_unstable();
		let mut merged: Vec<Interval> = Vec::with_capacity(intervals.len());
		for interval in intervals {
			match merged.last_mut() {
				Some(last) if last.contig == interval.contig && interval.start <= last.end => {
					last.end = last.end.max(interval.end);
				}
				_ => merged.push(interval),
			}
		}
		Self { intervals: merged }
	}

	pub(super) fn overlaps(&self, record: RecordRef<'_>) -> bool {
		let (Some(contig), Some(start), Some(end)) =
			(record.contig(), record.position(), record.reference_end())
		else {
			return false;
		};
		self.overlaps_span(contig, u64::from(start), end)
	}

	fn overlaps_span(&self, contig: ContigId, start: u64, end: u64) -> bool {
		let first = self
			.intervals
			.partition_point(|interval| (interval.contig, interval.end) <= (contig, start));
		self.intervals
			.get(first)
			.is_some_and(|interval| interval.contig == contig && interval.start < end)
	}
}

fn merge_chunks(mut chunks: Vec<Chunk>) -> Vec<Chunk> {
	chunks.retain(|chunk| chunk.start < chunk.end);
	chunks.sort_unstable();
	let mut merged: Vec<Chunk> = Vec::with_capacity(chunks.len());
	for chunk in chunks {
		match merged.last_mut() {
			Some(last) if chunk.start <= last.end => last.end = last.end.max(chunk.end),
			_ => merged.push(chunk),
		}
	}
	merged
}
