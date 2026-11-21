use multiversion::multiversion;

use crate::alignment::{CigarKind, CigarOp};
use crate::config::Scoring;
use crate::extend::{Extension, ExtensionBatch, ExtensionJob, ExtensionResults, Outcome};
use crate::nucleotide::CODE_PAD;
use crate::num::{u16_of, u32_of};

const LANES: usize = 32;
const NEGATIVE: i16 = -20_000;
const INACTIVE_QUERY: u8 = u8::MAX;
const PAD_PENALTY: i16 = 8;
const DIRECTION_MASK: u8 = 0b11;
const DIRECTION_START: u8 = 0;
const DIRECTION_DIAGONAL: u8 = 1;
const DIRECTION_INSERTION: u8 = 2;
const INSERTION_EXTENDS: u8 = 0b100;
const DELETION_EXTENDS: u8 = 0b1000;

type Lane<T> = [T; LANES];

#[derive(Debug, Clone, Copy)]
struct Penalties {
	matched: i16,
	mismatch: i16,
	ambiguous: i16,
	open: i16,
	extend: i16,
	clip: i16,
}

impl Penalties {
	fn new(scoring: &Scoring) -> Self {
		Self {
			matched: scoring.matched,
			mismatch: scoring.mismatch,
			ambiguous: scoring.ambiguous,
			open: scoring.gap_open + scoring.gap_extend,
			extend: scoring.gap_extend,
			clip: scoring.clip,
		}
	}

	fn substitution(&self, query: u8, reference: u8) -> i16 {
		if reference == CODE_PAD {
			-PAD_PENALTY
		} else if query > 3 || reference > 3 {
			-self.ambiguous
		} else if query == reference {
			self.matched
		} else {
			-self.mismatch
		}
	}
}

#[derive(Debug, Default)]
pub(crate) struct CpuScratch {
	queries: Vec<Lane<u8>>,
	references: Vec<Lane<u8>>,
	scores: Vec<Lane<i16>>,
	insertions: Vec<Lane<i16>>,
	trace: Vec<Lane<u8>>,
	order: Vec<u32>,
	path: Vec<CigarOp>,
	reversed: Vec<CigarOp>,
}

struct Tables<'a> {
	queries: &'a [Lane<u8>],
	references: &'a [Lane<u8>],
	lengths: Lane<u16>,
	bands: Lane<u16>,
	band: usize,
	penalties: Penalties,
	scores: &'a mut [Lane<i16>],
	insertions: &'a mut [Lane<i16>],
	trace: &'a mut [Lane<u8>],
}

struct Best {
	score: Lane<i16>,
	row: Lane<u16>,
	diagonal: Lane<u16>,
}

#[multiversion(targets(
	"x86_64+avx512f+avx512bw+avx512vl+avx512dq+avx512cd",
	"x86_64+avx2+bmi1+bmi2+lzcnt+fma",
	"aarch64+neon"
))]
fn fill(tables: &mut Tables<'_>) -> Best {
	let Tables {
		queries,
		references,
		lengths,
		bands,
		band,
		penalties,
		scores,
		insertions,
		trace,
	} = tables;
	let band = *band;
	let Penalties {
		matched,
		mismatch,
		ambiguous,
		open,
		extend,
		clip,
	} = *penalties;
	let mut best = Best {
		score: [0; LANES],
		row: [0; LANES],
		diagonal: [0; LANES],
	};
	for (row, query_row) in queries.iter().enumerate() {
		let mut row_bonus = [0i16; LANES];
		for (bonus, &length) in row_bonus.iter_mut().zip(lengths.iter()) {
			*bonus = match usize::from(length).cmp(&(row + 1)) {
				std::cmp::Ordering::Equal => clip,
				std::cmp::Ordering::Greater => 0,
				std::cmp::Ordering::Less => NEGATIVE,
			};
		}
		let mut left = [NEGATIVE; LANES];
		let mut deletion = [NEGATIVE; LANES];
		let mut row_best = [0i16; LANES];
		let mut row_diagonal = [0u16; LANES];
		let reference_row = &references[row..row + band];
		let trace_row = &mut trace[row * band..(row + 1) * band];
		for k in 0..band {
			let reference = reference_row[k];
			let previous = scores[k];
			let upper = scores[k + 1];
			let upper_insertion = insertions[k + 1];
			let mut value = [0i16; LANES];
			let mut insertion_out = [0i16; LANES];
			let mut deletion_out = [0i16; LANES];
			let mut flags = [0u8; LANES];
			for lane in 0..LANES {
				let q = query_row[lane];
				let r = reference[lane];
				let substitution = if r == CODE_PAD {
					-PAD_PENALTY
				} else if q > 3 || r > 3 {
					-ambiguous
				} else if q == r {
					matched
				} else {
					-mismatch
				};
				let diagonal = previous[lane].saturating_add(substitution);
				let insertion_open = upper[lane].saturating_sub(open);
				let insertion_extend = upper_insertion[lane].saturating_sub(extend);
				let insertion_extends = insertion_extend > insertion_open;
				let insertion = if insertion_extends {
					insertion_extend
				} else {
					insertion_open
				};
				let deletion_open = left[lane].saturating_sub(open);
				let deletion_extend = deletion[lane].saturating_sub(extend);
				let deletion_extends = deletion_extend > deletion_open;
				let deletion_value = if deletion_extends {
					deletion_extend
				} else {
					deletion_open
				};
				let mut cell = diagonal;
				let mut direction = DIRECTION_DIAGONAL;
				if insertion > cell {
					cell = insertion;
					direction = DIRECTION_INSERTION;
				}
				if deletion_value > cell {
					cell = deletion_value;
					direction = DIRECTION_INSERTION + 1;
				}
				if cell <= 0 {
					cell = 0;
					direction = DIRECTION_START;
				}
				let inside = k < usize::from(bands[lane]);
				let (cell, insertion, deletion_value) = if inside {
					(cell, insertion, deletion_value)
				} else {
					(NEGATIVE, NEGATIVE, NEGATIVE)
				};
				value[lane] = cell;
				insertion_out[lane] = insertion;
				deletion_out[lane] = deletion_value;
				flags[lane] = direction
					| (u8::from(insertion_extends) * INSERTION_EXTENDS)
					| (u8::from(deletion_extends) * DELETION_EXTENDS);
				let candidate = cell.saturating_add(row_bonus[lane]);
				if candidate > row_best[lane] {
					row_best[lane] = candidate;
					row_diagonal[lane] = u16_of(k);
				}
			}
			scores[k] = value;
			insertions[k] = insertion_out;
			left = value;
			deletion = deletion_out;
			trace_row[k] = flags;
		}
		for lane in 0..LANES {
			if row_best[lane] > best.score[lane] {
				best.score[lane] = row_best[lane];
				best.row[lane] = u16_of(row);
				best.diagonal[lane] = row_diagonal[lane];
			}
		}
	}
	best
}

pub(crate) fn extend(
	scoring: &Scoring,
	scratch: &mut CpuScratch,
	batch: &ExtensionBatch,
	results: &mut ExtensionResults,
) {
	results.reset(batch.len());
	let penalties = Penalties::new(scoring);
	let mut order = std::mem::take(&mut scratch.order);
	order.clear();
	order.extend(0..u32_of(batch.len()));
	order.sort_unstable_by_key(|&job| {
		let spec = batch.jobs().get(job as usize);
		(
			spec.map_or(0, ExtensionJob::band),
			spec.map_or(0, |spec| batch.query(spec).len()),
			job,
		)
	});
	for group in order.chunks(LANES) {
		run_group(scratch, &penalties, batch, group, results);
	}
	scratch.order = order;
}

fn run_group(
	scratch: &mut CpuScratch,
	penalties: &Penalties,
	batch: &ExtensionBatch,
	group: &[u32],
	results: &mut ExtensionResults,
) {
	let jobs: Vec<(u32, &ExtensionJob)> = group
		.iter()
		.filter_map(|&index| batch.jobs().get(index as usize).map(|job| (index, job)))
		.collect();
	let rows = jobs
		.iter()
		.map(|(_, job)| batch.query(job).len())
		.max()
		.unwrap_or(0);
	let band = jobs
		.iter()
		.map(|(_, job)| job.band())
		.max()
		.unwrap_or(1)
		.max(1);
	if rows == 0 {
		return;
	}
	scratch.queries.clear();
	scratch.queries.resize(rows, [INACTIVE_QUERY; LANES]);
	scratch.references.clear();
	scratch.references.resize(rows + band, [CODE_PAD; LANES]);
	let mut lengths = [0u16; LANES];
	let mut bands = [0u16; LANES];
	for (lane, (_, job)) in jobs.iter().enumerate() {
		let query = batch.query(job);
		lengths[lane] = u16_of(query.len());
		bands[lane] = job.band;
		for (row, &code) in scratch.queries.iter_mut().zip(query) {
			row[lane] = code;
		}
		for (row, &code) in scratch.references.iter_mut().zip(batch.window(job)) {
			row[lane] = code;
		}
	}
	scratch.scores.clear();
	scratch.scores.resize(band + 1, [penalties.clip; LANES]);
	for (lane, &limit) in bands.iter().enumerate() {
		for row in scratch.scores.iter_mut().skip(usize::from(limit)) {
			row[lane] = NEGATIVE;
		}
	}
	if let Some(sentinel) = scratch.scores.last_mut() {
		*sentinel = [NEGATIVE; LANES];
	}
	scratch.insertions.clear();
	scratch.insertions.resize(band + 1, [NEGATIVE; LANES]);
	scratch.trace.clear();
	scratch.trace.resize(rows * band, [0; LANES]);
	let best = fill(&mut Tables {
		queries: &scratch.queries,
		references: &scratch.references,
		lengths,
		bands,
		band,
		penalties: *penalties,
		scores: &mut scratch.scores,
		insertions: &mut scratch.insertions,
		trace: &mut scratch.trace,
	});
	for (lane, &(index, job)) in jobs.iter().enumerate() {
		let found = best.score.get(lane).copied().unwrap_or(0);
		if found <= 0 {
			continue;
		}
		let row = usize::from(best.row.get(lane).copied().unwrap_or(0));
		let diagonal = usize::from(best.diagonal.get(lane).copied().unwrap_or(0));
		let outcome = trace_lane(
			scratch,
			penalties,
			batch,
			job,
			lane,
			(row, diagonal),
			band,
			results,
		);
		results.set(index as usize, outcome);
	}
}

#[derive(Debug, Clone, Copy)]
enum State {
	Match,
	Insertion,
	Deletion,
}

fn push_run(path: &mut Vec<CigarOp>, kind: CigarKind) {
	match path.last_mut() {
		Some(last) if last.kind() == kind => *last = CigarOp::new(kind, last.length() + 1),
		_ => path.push(CigarOp::new(kind, 1)),
	}
}

#[allow(clippy::too_many_arguments)]
fn trace_lane(
	scratch: &mut CpuScratch,
	penalties: &Penalties,
	batch: &ExtensionBatch,
	job: &ExtensionJob,
	lane: usize,
	end: (usize, usize),
	band: usize,
	results: &mut ExtensionResults,
) -> Outcome {
	let (end_row, end_diagonal) = end;
	let (mut row, mut diagonal) = (end_row, end_diagonal);
	let mut state = State::Match;
	scratch.path.clear();
	let mut reference_end = None;
	let mut reference_span = 0usize;
	let mut query_start = end_row;
	while let Some(flags) = scratch
		.trace
		.get(row * band + diagonal)
		.and_then(|cell| cell.get(lane))
		.copied()
	{
		match state {
			State::Match => match flags & DIRECTION_MASK {
				DIRECTION_START => break,
				DIRECTION_DIAGONAL => {
					push_run(&mut scratch.path, CigarKind::Match);
					reference_end.get_or_insert(row + diagonal + 1);
					reference_span += 1;
					query_start = row;
					if row == 0 {
						break;
					}
					row -= 1;
				}
				DIRECTION_INSERTION => state = State::Insertion,
				_ => state = State::Deletion,
			},
			State::Insertion => {
				push_run(&mut scratch.path, CigarKind::Insertion);
				query_start = row;
				state = if flags & INSERTION_EXTENDS != 0 {
					State::Insertion
				} else {
					State::Match
				};
				if row == 0 || diagonal + 1 >= job.band() {
					break;
				}
				row -= 1;
				diagonal += 1;
			}
			State::Deletion => {
				push_run(&mut scratch.path, CigarKind::Deletion);
				reference_end.get_or_insert(row + diagonal + 1);
				reference_span += 1;
				state = if flags & DELETION_EXTENDS != 0 {
					State::Deletion
				} else {
					State::Match
				};
				if diagonal == 0 {
					break;
				}
				diagonal -= 1;
			}
		}
	}
	let Some(reference_end) = reference_end else {
		return Outcome::Unaligned;
	};
	let reference_start = reference_end.saturating_sub(reference_span);
	scratch.reversed.clear();
	scratch.reversed.extend(scratch.path.iter().rev());
	let query = batch.query(job);
	let window = batch.window(job);
	let score = path_score(
		penalties,
		&scratch.reversed,
		query.get(query_start..).unwrap_or_default(),
		window.get(reference_start..).unwrap_or_default(),
	);
	if score <= 0 {
		return Outcome::Unaligned;
	}
	let operations = results.push_operations(&scratch.reversed);
	Outcome::Aligned(Extension::new(
		score,
		u32_of(query_start)..u32_of(end_row + 1),
		u32_of(reference_start)..u32_of(reference_end),
		operations,
	))
}

fn path_score(
	penalties: &Penalties,
	operations: &[CigarOp],
	query: &[u8],
	reference: &[u8],
) -> i32 {
	let (mut q, mut r) = (0usize, 0usize);
	let mut score = 0i32;
	for operation in operations {
		let length = operation.length() as usize;
		match operation.kind() {
			CigarKind::Match => {
				let pairs = query
					.iter()
					.skip(q)
					.zip(reference.iter().skip(r))
					.take(length);
				score += pairs
					.map(|(&a, &b)| i32::from(penalties.substitution(a, b)))
					.sum::<i32>();
				q += length;
				r += length;
			}
			CigarKind::Insertion => {
				score -= i32::from(penalties.open)
					+ i32::from(penalties.extend) * (i32::from(u16_of(length)) - 1);
				q += length;
			}
			CigarKind::Deletion => {
				score -= i32::from(penalties.open)
					+ i32::from(penalties.extend) * (i32::from(u16_of(length)) - 1);
				r += length;
			}
			CigarKind::SoftClip => q += length,
		}
	}
	score
}
