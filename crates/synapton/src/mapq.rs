use crate::candidate::Candidate;
use crate::config::Scoring;
use crate::num::round_i32;

const MIN_SEED_LENGTH: i32 = 19;
const LONG_READ_LENGTH: f64 = 50.0;
const MAX_MAPQ: i32 = 60;
const MAPQ_SCALE: f64 = 6.02;
const REPEAT_SCALE: f64 = 4.343;
const ROUNDING: f64 = 0.499;
const OVERLAP_DIVISOR: u32 = 2;

#[derive(Debug, Clone, Copy)]
pub(crate) struct Competition {
	pub(crate) best_other: Option<i32>,
	pub(crate) near_best: u32,
}

pub(crate) fn competition(
	scoring: &Scoring,
	candidates: &[Candidate],
	chosen: usize,
	read_length: u32,
) -> Competition {
	let Some(primary) = candidates.get(chosen) else {
		return Competition {
			best_other: None,
			near_best: 0,
		};
	};
	let span = primary.read_span(read_length);
	let competes = |index: usize, other: &Candidate| {
		let other_span = other.read_span(read_length);
		let overlap = span
			.end
			.min(other_span.end)
			.saturating_sub(span.start.max(other_span.start));
		let shorter = (span.end - span.start).min(other_span.end - other_span.start);
		index != chosen && overlap * OVERLAP_DIVISOR >= shorter
	};
	let best_other = candidates
		.iter()
		.enumerate()
		.filter(|&(index, other)| competes(index, other))
		.map(|(_, other)| other.score)
		.max();
	let threshold = best_other.map_or(0, |best| {
		best - i32::from(scoring.matched + scoring.mismatch)
	});
	let near_best = candidates
		.iter()
		.enumerate()
		.filter(|&(index, other)| competes(index, other) && other.score >= threshold)
		.count();
	Competition {
		best_other,
		near_best: u32::try_from(near_best.saturating_sub(1)).unwrap_or(u32::MAX),
	}
}

pub(crate) fn single_end(
	scoring: &Scoring,
	candidate: &Candidate,
	competition: &Competition,
	repeat_fraction: f32,
) -> u8 {
	let matched = i32::from(scoring.matched);
	let sub = competition.best_other.unwrap_or(MIN_SEED_LENGTH * matched);
	if sub >= candidate.score {
		return 0;
	}
	let length = f64::from(candidate.query_length().max(candidate.reference_length()));
	let penalty = f64::from(matched + i32::from(scoring.mismatch));
	let identity =
		1.0 - (length * f64::from(matched) - f64::from(candidate.score)) / penalty / length;
	let coefficient = if length < LONG_READ_LENGTH {
		1.0
	} else {
		LONG_READ_LENGTH.ln() / length.ln()
	};
	let factor = coefficient * identity * identity;
	let gap = f64::from(candidate.score - sub) / f64::from(matched);
	let mut mapq = (MAPQ_SCALE * gap * factor * factor + ROUNDING).floor();
	if competition.near_best > 0 {
		mapq -= (REPEAT_SCALE * f64::from(competition.near_best + 1).ln() + ROUNDING).floor();
	}
	let clamped = mapq.clamp(0.0, f64::from(MAX_MAPQ));
	let scaled = (clamped * (1.0 - f64::from(repeat_fraction)) + ROUNDING).floor();
	to_mapq(scaled)
}

pub(crate) fn paired(
	scoring: &Scoring,
	best_pair: i32,
	runner_up: i32,
	competing: u32,
	repeat_fraction: f32,
) -> u8 {
	let matched = i32::from(scoring.matched);
	let gap = f64::from(best_pair - runner_up) / f64::from(matched);
	let mut mapq = (MAPQ_SCALE * gap + ROUNDING).floor();
	if competing > 0 {
		mapq -= (REPEAT_SCALE * f64::from(competing + 1).ln() + ROUNDING).floor();
	}
	let clamped = mapq.clamp(0.0, f64::from(MAX_MAPQ));
	to_mapq((clamped * (1.0 - f64::from(repeat_fraction)) + ROUNDING).floor())
}

fn to_mapq(value: f64) -> u8 {
	let clamped = value.clamp(0.0, f64::from(MAX_MAPQ));
	u8::try_from(round_i32(clamped)).unwrap_or(0)
}
