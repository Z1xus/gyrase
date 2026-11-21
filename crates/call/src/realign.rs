use std::ops::Range;

use gyrase_pileup::{Pileup, Read};

use crate::genotype::MAX_ALLELES;
use crate::likelihood::{Emissions, Evidence, phred_to_ln};

const PAD: u64 = 8;
const FLANK: u64 = 40;
const SLACK: usize = 2;
const MAX_READS: usize = 200;
const DECISIVE: f64 = -std::f64::consts::LN_10;

struct Weights {
	matched: Vec<f64>,
	mismatched: Vec<f64>,
}

impl Weights {
	fn load(&mut self, emissions: &Emissions, qualities: &[u8]) {
		self.matched.clear();
		self.matched
			.extend(qualities.iter().map(|&quality| emissions.matched(quality)));
		self.mismatched.clear();
		self.mismatched.extend(
			qualities
				.iter()
				.map(|&quality| emissions.mismatched(quality)),
		);
	}

	fn placement(&self, bases: &[u8], window: &[u8]) -> f64 {
		window
			.iter()
			.zip(bases)
			.zip(self.matched.iter().zip(&self.mismatched))
			.map(|((expected, base), (matched, mismatched))| {
				if expected == base {
					*matched
				} else {
					*mismatched
				}
			})
			.sum()
	}

	fn best(&self, bases: &[u8], haplotype: &[u8], expected: usize) -> f64 {
		(expected.saturating_sub(SLACK)..=expected + SLACK)
			.filter_map(|offset| haplotype.get(offset..offset + bases.len()))
			.map(|window| self.placement(bases, window))
			.fold(f64::NEG_INFINITY, f64::max)
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Edit<'a> {
	pub start: u64,
	pub deleted: u64,
	pub inserted: &'a [u8],
}

impl Edit<'_> {
	fn shift(self) -> isize {
		let inserted = isize::try_from(self.inserted.len()).unwrap_or(isize::MAX);
		inserted - isize::try_from(self.deleted).unwrap_or(isize::MAX)
	}
}

pub(crate) fn realign(
	pileup: &Pileup<'_>,
	emissions: &Emissions,
	focus: Range<u64>,
	alternates: &[Edit<'_>],
	competitors: &[Edit<'_>],
	cap: impl Fn(&Read<'_>) -> Option<u8>,
) -> Vec<Evidence> {
	let reach = alternates
		.iter()
		.map(|edit| edit.start + edit.deleted)
		.fold(focus.end, u64::max);
	let from = focus.start.saturating_sub(FLANK).max(pileup.window.start);
	let to = reach + FLANK;
	let longest_deletion = alternates
		.iter()
		.chain(competitors)
		.map(|edit| edit.deleted)
		.max()
		.unwrap_or(0);
	let begin = from.saturating_sub(PAD).max(pileup.window.start);
	let end = to + PAD + longest_deletion;
	let reference = Haplotype::new(pileup, begin, end, Vec::new());
	let competitors: Vec<Edit<'_>> = competitors
		.iter()
		.copied()
		.filter(|edit| (begin..end).contains(&edit.start))
		.collect();
	let alternates: Vec<Vec<Haplotype<'_>>> = alternates
		.iter()
		.map(|&alternate| {
			std::iter::once(vec![alternate])
				.chain(
					competitors
						.iter()
						.filter(|competitor| apart(alternate, **competitor))
						.map(|&competitor| {
							let mut edits = vec![alternate, competitor];
							edits.sort_by_key(|edit| edit.start);
							edits
						}),
				)
				.map(|edits| Haplotype::new(pileup, begin, end, edits))
				.collect()
		})
		.collect();
	let competitors: Vec<Haplotype<'_>> = competitors
		.iter()
		.map(|&edit| Haplotype::new(pileup, begin, end, vec![edit]))
		.filter(|competitor| {
			competitor.sequence != reference.sequence
				&& alternates
					.iter()
					.flatten()
					.all(|alternate| alternate.sequence != competitor.sequence)
		})
		.collect();
	let mut weights = Weights {
		matched: Vec::new(),
		mismatched: Vec::new(),
	};
	informative(pileup, focus.start)
		.iter()
		.filter_map(|read| {
			let floor = phred_to_ln(cap(read)?);
			let (start, bases, qualities) = read.slice(from, to)?;
			let expected = usize::try_from(start - begin).ok()?;
			weights.load(emissions, qualities);
			let score =
				|haplotype: &Haplotype<'_>| haplotype.score(&weights, bases, start, expected);
			let mut likelihoods = [f64::NEG_INFINITY; MAX_ALLELES];
			likelihoods[0] = std::iter::once(&reference)
				.chain(&competitors)
				.map(score)
				.fold(f64::NEG_INFINITY, f64::max);
			for (slot, haplotypes) in likelihoods[1..].iter_mut().zip(&alternates) {
				*slot = haplotypes
					.iter()
					.map(score)
					.fold(f64::NEG_INFINITY, f64::max);
			}
			let used = &mut likelihoods[..=alternates.len()];
			let peak = used.iter().copied().fold(f64::NEG_INFINITY, f64::max);
			if !peak.is_finite() {
				return None;
			}
			for value in used.iter_mut() {
				*value = (*value - peak).max(floor);
			}
			Some(Evidence { likelihoods })
		})
		.collect()
}

pub(crate) fn depths(evidence: &[Evidence], alleles: usize) -> Vec<u32> {
	let mut depths = vec![0; alleles];
	for item in evidence {
		if let Some(winner) = decisive(&item.likelihoods[..alleles]) {
			depths[winner] += 1;
		}
	}
	depths
}

struct Haplotype<'a> {
	sequence: Vec<u8>,
	edits: Vec<Edit<'a>>,
}

impl<'a> Haplotype<'a> {
	fn new(pileup: &Pileup<'_>, begin: u64, end: u64, edits: Vec<Edit<'a>>) -> Self {
		let window = pileup.window;
		let offset = |position: u64| {
			usize::try_from(position.saturating_sub(window.start))
				.unwrap_or(usize::MAX)
				.min(window.bases.len())
		};
		let slice = |from: u64, to: u64| {
			window
				.bases
				.get(offset(from)..offset(to).max(offset(from)))
				.unwrap_or_default()
		};
		let mut sequence = Vec::new();
		let mut cursor = begin;
		for edit in &edits {
			sequence.extend_from_slice(slice(cursor, edit.start));
			sequence.extend_from_slice(edit.inserted);
			cursor = edit.start + edit.deleted;
		}
		sequence.extend_from_slice(slice(cursor, end));
		Self { sequence, edits }
	}

	fn score(&self, weights: &Weights, bases: &[u8], start: u64, expected: usize) -> f64 {
		let shift: isize = self
			.edits
			.iter()
			.filter(|edit| start >= edit.start)
			.map(|edit| edit.shift())
			.sum();
		let shifted = (shift != 0)
			.then(|| expected.checked_add_signed(shift))
			.flatten();
		std::iter::once(expected)
			.chain(shifted)
			.map(|center| weights.best(bases, &self.sequence, center))
			.fold(f64::NEG_INFINITY, f64::max)
	}
}

fn apart(first: Edit<'_>, second: Edit<'_>) -> bool {
	first.start + first.deleted < second.start || second.start + second.deleted < first.start
}

fn decisive(likelihoods: &[f64]) -> Option<usize> {
	let (winner, &peak) = likelihoods
		.iter()
		.enumerate()
		.max_by(|left, right| left.1.total_cmp(right.1))?;
	likelihoods
		.iter()
		.enumerate()
		.all(|(index, &value)| index == winner || value <= peak + DECISIVE)
		.then_some(winner)
}

fn informative<'a>(pileup: &Pileup<'a>, position: u64) -> Vec<Read<'a>> {
	let mut reads: Vec<Read<'a>> = pileup.reads_covering(position).collect();
	reads.sort_by_key(|read| (read.fragment, read.start, read.clipped_start));
	reads.dedup_by(|later, earlier| later.fragment.is_some() && later.fragment == earlier.fragment);
	let stride = reads.len().div_ceil(MAX_READS).max(1);
	reads.into_iter().step_by(stride).collect()
}
