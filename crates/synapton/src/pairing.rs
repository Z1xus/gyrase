use crate::alignment::Strand;
use crate::candidate::Candidate;
use crate::config::Scoring;
use crate::num::{float_len, float_of, round_i32, round_i64, u32_of};

const DEFAULT_MEAN: f64 = 400.0;
const DEFAULT_STD: f64 = 100.0;
const MIN_SAMPLES: usize = 100;
const OUTLIER_BOUND: f64 = 2.0;
const MAX_BOUND: f64 = 3.0;
const MAX_DEVIATIONS: f64 = 4.0;
const LIKELIHOOD_SCALE: f64 = 0.721;
const ROUNDING: f64 = 0.499;
const MAX_PAIR_CANDIDATES: usize = 8;

#[derive(Debug, Clone, Copy)]
pub(crate) struct InsertModel {
	mean: f64,
	std: f64,
	low: i64,
	high: i64,
}

impl Default for InsertModel {
	fn default() -> Self {
		Self::from_moments(DEFAULT_MEAN, DEFAULT_STD)
	}
}

impl InsertModel {
	fn from_moments(mean: f64, std: f64) -> Self {
		let low = round_i64((mean - MAX_DEVIATIONS * std + ROUNDING).floor());
		let high = round_i64((mean + MAX_DEVIATIONS * std + ROUNDING).floor());
		Self {
			mean,
			std,
			low: low.max(1),
			high,
		}
	}

	pub(crate) fn estimate(distances: &mut [i64]) -> Self {
		if distances.len() < MIN_SAMPLES {
			return Self::default();
		}
		distances.sort_unstable();
		let quantile = |fraction: f64| {
			let index = usize::try_from(round_i64((float_len(distances.len()) * fraction).floor()))
				.unwrap_or(0);
			float_of(distances.get(index).copied().unwrap_or(0))
		};
		let (p25, p75) = (quantile(0.25), quantile(0.75));
		let spread = p75 - p25;
		let (inner_low, inner_high) = (p25 - OUTLIER_BOUND * spread, p75 + OUTLIER_BOUND * spread);
		let kept: Vec<f64> = distances
			.iter()
			.map(|&distance| float_of(distance))
			.filter(|&distance| distance >= inner_low && distance <= inner_high)
			.collect();
		if kept.len() < MIN_SAMPLES {
			return Self::default();
		}
		let mean = kept.iter().sum::<f64>() / float_len(kept.len());
		let variance =
			kept.iter().map(|value| (value - mean).powi(2)).sum::<f64>() / float_len(kept.len());
		let std = variance.sqrt().max(1.0);
		let mut model = Self::from_moments(mean, std);
		let outer_low = round_i64(p25 - MAX_BOUND * spread);
		let outer_high = round_i64(p75 + MAX_BOUND * spread);
		model.low = model.low.min(outer_low).max(1);
		model.high = model.high.max(outer_high);
		model
	}

	pub(crate) fn contains(&self, distance: i64) -> bool {
		distance >= self.low && distance <= self.high
	}

	pub(crate) fn window(&self) -> (i64, i64) {
		(self.low, self.high)
	}

	fn bonus(&self, scoring: &Scoring, distance: i64) -> i32 {
		let deviation = (float_of(distance) - self.mean).abs() / self.std;
		let tail = (2.0 * erfc(deviation / std::f64::consts::SQRT_2)).max(f64::MIN_POSITIVE);
		round_i32(LIKELIHOOD_SCALE * tail.ln() * f64::from(scoring.matched) + ROUNDING)
	}
}

fn erfc(x: f64) -> f64 {
	let z = x.abs();
	let t = 1.0 / (1.0 + 0.5 * z);
	let polynomial = -z * z - 1.265_512_23
		+ t * (1.000_023_68
			+ t * (0.374_091_96
				+ t * (0.096_784_18
					+ t * (-0.186_288_06
						+ t * (0.278_868_07
							+ t * (-1.135_203_98
								+ t * (1.488_515_87 + t * (-0.822_152_23 + t * 0.170_872_77))))))));
	let value = t * polynomial.exp();
	if x >= 0.0 { value } else { 2.0 - value }
}

pub(crate) fn distance(first: &Candidate, second: &Candidate) -> Option<i64> {
	if first.contig != second.contig || first.strand == second.strand {
		return None;
	}
	let (forward, reverse) = match first.strand {
		Strand::Forward => (first, second),
		Strand::Reverse => (second, first),
	};
	let span = i64::from(reverse.stop) - i64::from(forward.start);
	(span > 0 && reverse.start >= forward.start.saturating_sub(1)).then_some(span)
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct PairChoice {
	pub(crate) first: usize,
	pub(crate) second: usize,
	pub(crate) score: i32,
	pub(crate) runner_up: i32,
	pub(crate) competing: u32,
}

pub(crate) fn choose_pair(
	scoring: &Scoring,
	model: &InsertModel,
	first: &[Candidate],
	second: &[Candidate],
) -> Option<PairChoice> {
	let unpaired = first.first()?.score + second.first()?.score - scoring.unpaired;
	let scored = || {
		first
			.iter()
			.take(MAX_PAIR_CANDIDATES)
			.enumerate()
			.flat_map(|(left_index, left)| {
				second
					.iter()
					.take(MAX_PAIR_CANDIDATES)
					.enumerate()
					.filter_map(move |(right_index, right)| {
						let span = distance(left, right).filter(|&span| model.contains(span))?;
						let score = left.score + right.score + model.bonus(scoring, span);
						Some((left_index, right_index, score))
					})
			})
	};
	let (left, right, score) = scored().fold(None, |best, found| match best {
		Some((_, _, top)) if found.2 <= top => best,
		_ => Some(found),
	})?;
	let near = i32::from(scoring.matched + scoring.mismatch);
	let others = || {
		scored()
			.filter(move |&(other_left, other_right, _)| (other_left, other_right) != (left, right))
			.map(|(_, _, other)| other)
	};
	(score > unpaired).then(|| PairChoice {
		first: left,
		second: right,
		score,
		runner_up: others().max().unwrap_or(i32::MIN).max(unpaired),
		competing: u32_of(others().filter(|&other| other + near >= score).count()),
	})
}

pub(crate) fn compatible(model: &InsertModel, anchor: &Candidate, mates: &[Candidate]) -> bool {
	mates
		.iter()
		.any(|mate| distance(anchor, mate).is_some_and(|span| model.contains(span)))
}
