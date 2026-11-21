use std::cmp::Reverse;

use gyrase_genome::Interval;
use gyrase_io::{BlockSpan, Chunk};

use crate::Index;

pub const DEFAULT_MAX_GAP: u64 = 4 << 20;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Plan {
	chunks: Vec<Chunk>,
	spans: Vec<BlockSpan>,
}

impl Plan {
	pub fn chunks(&self) -> &[Chunk] {
		&self.chunks
	}

	pub fn spans(&self) -> &[BlockSpan] {
		&self.spans
	}
}

impl Index {
	pub fn plan(&self, intervals: &[Interval], max_gap: u64) -> Plan {
		let mut bins = Vec::new();
		let mut chunks = Vec::new();
		for interval in intervals {
			if let Some(reference) = self.references().get(interval.contig.index()) {
				reference.collect_chunks(
					self.geometry(),
					interval.start,
					interval.end,
					&mut bins,
					&mut chunks,
				);
			}
		}
		let chunks = merge(chunks);
		let spans = coalesce(&chunks, max_gap);
		Plan { chunks, spans }
	}
}

fn merge(mut chunks: Vec<Chunk>) -> Vec<Chunk> {
	chunks.sort_unstable_by_key(|chunk| (chunk.begin, Reverse(chunk.end)));
	let mut merged: Vec<Chunk> = Vec::with_capacity(chunks.len());
	for chunk in chunks {
		let Some(last) = merged.last_mut() else {
			merged.push(chunk);
			continue;
		};
		if last.end >= chunk.end {
			continue;
		}
		if last.end >= chunk.begin {
			last.end = chunk.begin;
		}
		if last.end.block() == chunk.begin.block() {
			last.end = chunk.end;
		} else {
			merged.push(chunk);
		}
	}
	merged
}

fn coalesce(chunks: &[Chunk], max_gap: u64) -> Vec<BlockSpan> {
	let mut spans: Vec<BlockSpan> = Vec::new();
	for chunk in chunks {
		let next = BlockSpan {
			first: chunk.begin.block(),
			last: chunk.end.block(),
		};
		match spans.last_mut() {
			Some(span) if next.first.saturating_sub(span.last) <= max_gap => {
				span.last = span.last.max(next.last);
			}
			_ => spans.push(next),
		}
	}
	spans
}
