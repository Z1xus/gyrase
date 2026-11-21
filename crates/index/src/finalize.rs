use std::collections::HashMap;

use gyrase_io::{Chunk, VirtualOffset};

use crate::order::bucket_order;
use crate::reference::fill_loff;
use crate::{Bin, Error, Geometry, RefIndex};

const MIN_MARKER_BLOCKS: u64 = 0x10000;
const UNSET: u64 = u64::MAX;

#[derive(Debug)]
pub(crate) struct Draft {
	pub(crate) bins: Vec<(u32, Vec<Chunk>)>,
	pub(crate) puts: Vec<u32>,
	pub(crate) linear: Vec<u64>,
}

pub(crate) fn finalize(geometry: Geometry, draft: Draft) -> Result<RefIndex, Error> {
	let linear = fill_gaps(draft.linear);
	let positions: HashMap<u32, usize> = draft
		.bins
		.iter()
		.enumerate()
		.map(|(position, (id, _))| (*id, position))
		.collect();
	let mut bins: Vec<Bin> = draft
		.bins
		.into_iter()
		.map(|(id, chunks)| Bin {
			id,
			loff: VirtualOffset::default(),
			chunks,
		})
		.collect();
	fill_loff(geometry, &mut bins, &linear);
	let mut slots: Vec<Option<Bin>> = bins.into_iter().map(Some).collect();
	merge_small_bins(geometry, &mut slots, &positions);
	tidy_chunks(geometry, &mut slots);
	let bins = bucket_order(draft.puts)
		.into_iter()
		.filter_map(|id| slots.get_mut(*positions.get(&id)?)?.take())
		.collect();
	RefIndex::new(bins, linear)
}

fn fill_gaps(mut linear: Vec<u64>) -> Vec<VirtualOffset> {
	let mut following = None;
	for entry in linear.iter_mut().rev() {
		if *entry == UNSET
			&& let Some(next) = following
		{
			*entry = next;
		}
		following = Some(*entry);
	}
	linear.into_iter().map(VirtualOffset::from_raw).collect()
}

fn is_small(chunks: &[Chunk]) -> bool {
	match (chunks.first(), chunks.last()) {
		(Some(first), Some(last)) => {
			last.end.block().wrapping_sub(first.begin.block()) < MIN_MARKER_BLOCKS
		}
		_ => false,
	}
}

fn merge_small_bins(
	geometry: Geometry,
	slots: &mut [Option<Bin>],
	positions: &HashMap<u32, usize>,
) {
	for level in (1..=geometry.depth()).rev() {
		for position in 0..slots.len() {
			let Some(bin) = slots.get_mut(position).and_then(Option::as_mut) else {
				continue;
			};
			if geometry.level(bin.id) != Some(level) {
				continue;
			}
			if level < geometry.depth() {
				bin.chunks.sort_unstable_by_key(|chunk| chunk.begin);
			}
			let small = is_small(&bin.chunks);
			let parent = bin.id.saturating_sub(1) >> 3;
			let Some(&parent_position) = positions.get(&parent) else {
				continue;
			};
			if !small || slots.get(parent_position).is_none_or(Option::is_none) {
				continue;
			}
			let Some(child) = slots.get_mut(position).and_then(Option::take) else {
				continue;
			};
			if let Some(parent) = slots.get_mut(parent_position).and_then(Option::as_mut) {
				parent.chunks.extend(child.chunks);
			}
		}
	}
}

fn tidy_chunks(geometry: Geometry, slots: &mut [Option<Bin>]) {
	for bin in slots.iter_mut().flatten() {
		if bin.id == 0 {
			bin.chunks.sort_unstable_by_key(|chunk| chunk.begin);
		}
		if bin.id < geometry.bin_count() {
			merge_same_block(&mut bin.chunks);
		}
	}
}

fn merge_same_block(chunks: &mut Vec<Chunk>) {
	let mut merged: Vec<Chunk> = Vec::with_capacity(chunks.len());
	for chunk in chunks.drain(..) {
		match merged.last_mut() {
			Some(last) if last.end.block() >= chunk.begin.block() => {
				last.end = last.end.max(chunk.end);
			}
			_ => merged.push(chunk),
		}
	}
	*chunks = merged;
}
