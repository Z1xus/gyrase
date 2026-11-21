use crate::alignment::Strand;
use crate::num::i64_of;

const QUERY_BITS: u32 = 20;
const QUERY_MASK: u64 = (1 << QUERY_BITS) - 1;
const MAX_LOOKBACK: usize = 50;
const MAX_DIAGONAL_DIFFERENCE: u64 = 50;
const MAX_REFERENCE_GAP: u64 = 1024;
const DIAGONAL_PENALTY: f32 = 0.1;
const MIN_CHAIN_SCORE: f32 = 1.5;
const NO_PARENT: usize = usize::MAX;

#[derive(Debug, Clone, Copy)]
pub(crate) struct Chain {
	pub(crate) strand: Strand,
	pub(crate) score: f32,
	pub(crate) reference_start: u64,
	pub(crate) reference_end: u64,
	pub(crate) diagonal_low: i64,
	pub(crate) diagonal_high: i64,
	pub(crate) diagonal_middle: i64,
}

impl Chain {
	pub(crate) fn diagonal_spread(&self) -> usize {
		usize::try_from(self.diagonal_high - self.diagonal_low).unwrap_or(0)
	}
}

pub(crate) fn pack_anchor(reference: u64, query: u32) -> u64 {
	(reference << QUERY_BITS) | u64::from(query)
}

const fn unpack(anchor: u64) -> (u64, u64) {
	(anchor >> QUERY_BITS, anchor & QUERY_MASK)
}

#[derive(Debug, Default)]
pub(crate) struct Chainer {
	scores: Vec<f32>,
	parents: Vec<usize>,
	used: Vec<bool>,
	order: Vec<usize>,
	members: Vec<usize>,
	diagonals: Vec<i64>,
}

impl Chainer {
	pub(crate) fn chain(
		&mut self,
		anchors: &mut Vec<u64>,
		strand: Strand,
		seed_length: u32,
		chains: &mut Vec<Chain>,
	) {
		anchors.sort_unstable();
		anchors.dedup();
		let count = anchors.len();
		self.scores.clear();
		self.parents.clear();
		self.used.clear();
		self.used.resize(count, false);
		let length = f32::from(u16::try_from(seed_length).unwrap_or(u16::MAX));
		for (index, &anchor) in anchors.iter().enumerate() {
			let (reference, query) = unpack(anchor);
			let mut best = length;
			let mut parent = NO_PARENT;
			let first = index.saturating_sub(MAX_LOOKBACK);
			for (offset, &earlier) in anchors
				.get(first..index)
				.unwrap_or_default()
				.iter()
				.enumerate()
				.rev()
			{
				let (earlier_reference, earlier_query) = unpack(earlier);
				let reference_gap = reference - earlier_reference;
				if reference_gap > MAX_REFERENCE_GAP {
					break;
				}
				if earlier_query >= query || reference_gap == 0 {
					continue;
				}
				let query_gap = query - earlier_query;
				let difference = reference_gap.abs_diff(query_gap);
				if difference > MAX_DIAGONAL_DIFFERENCE {
					continue;
				}
				let gain = f32::from(
					u16::try_from(reference_gap.min(query_gap).min(u64::from(seed_length)))
						.unwrap_or(0),
				);
				let earlier_score = self.scores.get(first + offset).copied().unwrap_or(0.0);
				let candidate = earlier_score + gain
					- DIAGONAL_PENALTY * f32::from(u16::try_from(difference).unwrap_or(u16::MAX));
				if candidate > best {
					best = candidate;
					parent = first + offset;
				}
			}
			self.scores.push(best);
			self.parents.push(parent);
		}
		self.order.clear();
		self.order.extend(0..count);
		let scores = &self.scores;
		self.order.sort_unstable_by(|&left, &right| {
			let left_score = scores.get(left).copied().unwrap_or(0.0);
			let right_score = scores.get(right).copied().unwrap_or(0.0);
			right_score.total_cmp(&left_score).then(left.cmp(&right))
		});
		for position in 0..self.order.len() {
			let Some(&end) = self.order.get(position) else {
				break;
			};
			if self.used.get(end).copied().unwrap_or(true) {
				continue;
			}
			self.members.clear();
			let mut current = end;
			let mut inherited = 0.0;
			loop {
				if let Some(flag) = self.used.get_mut(current) {
					*flag = true;
				}
				self.members.push(current);
				let parent = self.parents.get(current).copied().unwrap_or(NO_PARENT);
				if parent == NO_PARENT {
					break;
				}
				if self.used.get(parent).copied().unwrap_or(true) {
					inherited = self.scores.get(parent).copied().unwrap_or(0.0);
					break;
				}
				current = parent;
			}
			let score = self.scores.get(end).copied().unwrap_or(0.0) - inherited;
			if score < MIN_CHAIN_SCORE * length {
				continue;
			}
			if let Some(chain) = summarize(
				anchors,
				&self.members,
				&mut self.diagonals,
				strand,
				score,
				seed_length,
			) {
				chains.push(chain);
			}
		}
	}
}

fn summarize(
	anchors: &[u64],
	members: &[usize],
	diagonals: &mut Vec<i64>,
	strand: Strand,
	score: f32,
	seed_length: u32,
) -> Option<Chain> {
	let first = unpack(*anchors.get(*members.last()?)?);
	let last = unpack(*anchors.get(*members.first()?)?);
	let mut diagonal_low = i64::MAX;
	let mut diagonal_high = i64::MIN;
	diagonals.clear();
	for &member in members {
		let (reference, query) = unpack(*anchors.get(member)?);
		let diagonal = i64_of(reference) - i64_of(query);
		diagonal_low = diagonal_low.min(diagonal);
		diagonal_high = diagonal_high.max(diagonal);
		diagonals.push(diagonal);
	}
	diagonals.sort_unstable();
	Some(Chain {
		strand,
		score,
		reference_start: first.0,
		reference_end: last.0 + u64::from(seed_length),
		diagonal_low,
		diagonal_high,
		diagonal_middle: *diagonals.get(diagonals.len() / 2)?,
	})
}
