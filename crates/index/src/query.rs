use gyrase_io::VirtualOffset;

use crate::{Bin, Chunk, Geometry, RefIndex};

const CHILDREN_PER_BIN_BITS: u32 = 3;
const FIRST_CHILD_RESIDUE: u32 = 1;

impl RefIndex {
	pub(crate) fn collect_chunks(
		&self,
		geometry: Geometry,
		start: u64,
		end: u64,
		bins: &mut Vec<u32>,
		out: &mut Vec<Chunk>,
	) {
		if self.bins().is_empty() || start >= end || start >= geometry.max_position() {
			return;
		}
		let min_offset = self.min_offset(geometry, start);
		let max_offset = self.max_offset(geometry, end);
		geometry.bins(start, end, bins);
		for bin in bins.iter().filter_map(|&id| self.bin(id)) {
			out.extend(bin.chunks.iter().filter_map(|chunk| {
				let begin = chunk.begin.raw().max(min_offset);
				let end = chunk.end.raw().min(max_offset);
				(begin < end).then(|| Chunk {
					begin: VirtualOffset::from_raw(begin),
					end: VirtualOffset::from_raw(end),
				})
			}));
		}
	}

	fn anchor_bin(&self, geometry: Geometry, window: u64) -> Option<&Bin> {
		let mut bin = u64::from(geometry.first_leaf()) + window;
		loop {
			if let Some(found) = u32::try_from(bin).ok().and_then(|id| self.bin(id)) {
				return Some(found);
			}
			let parent = bin.saturating_sub(1) >> CHILDREN_PER_BIN_BITS;
			let first_sibling = (parent << CHILDREN_PER_BIN_BITS) + 1;
			bin = if bin > first_sibling { bin - 1 } else { parent };
			if bin == 0 {
				return self.bin(0);
			}
		}
	}

	fn min_offset(&self, geometry: Geometry, start: u64) -> u64 {
		let window = start >> geometry.min_shift();
		let anchor = self.anchor_bin(geometry, window);
		let anchor_loff = anchor.map_or(0, |bin| bin.loff.raw());
		let Some(entry) = usize::try_from(window)
			.ok()
			.and_then(|window| self.linear().get(window))
		else {
			return anchor_loff;
		};
		let mut min_offset = anchor_loff.max(entry.raw());
		let has_unmapped = self
			.metadata(geometry)
			.is_none_or(|metadata| metadata.unmapped != 0);
		if !has_unmapped {
			return min_offset;
		}
		let earlier = usize::try_from(window)
			.ok()
			.and_then(|window| self.linear().get(..window))
			.and_then(|before| before.iter().rev().find(|offset| offset.raw() < min_offset));
		if let Some(offset) = earlier {
			min_offset = offset.raw();
		}
		match anchor {
			Some(bin) if earlier.is_none() || min_offset < bin.loff.raw() => bin.loff.raw(),
			_ => min_offset,
		}
	}

	fn max_offset(&self, geometry: Geometry, end: u64) -> u64 {
		if end > geometry.max_position() {
			return u64::MAX;
		}
		let first_leaf = u64::from(geometry.first_leaf());
		let mut bin = first_leaf + (end.saturating_sub(1) >> geometry.min_shift()) + 1;
		if bin >= u64::from(geometry.bin_count()) {
			bin = 0;
		}
		loop {
			while bin % (1 << CHILDREN_PER_BIN_BITS) == u64::from(FIRST_CHILD_RESIDUE) {
				bin = bin.saturating_sub(1) >> CHILDREN_PER_BIN_BITS;
			}
			if bin == 0 {
				return u64::MAX;
			}
			if let Some(first) = u32::try_from(bin)
				.ok()
				.and_then(|id| self.bin(id))
				.and_then(|found| found.chunks.first())
			{
				return first.begin.raw();
			}
			bin += 1;
		}
	}
}
