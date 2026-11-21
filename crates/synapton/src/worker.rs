use std::ops::Range;

use gyrase_genome::ContigId;

use crate::alignment::{
	AlignedBatch, Alignment, CigarKind, CigarOp, End, Flag, FlagSet, Mate, MateLocation, Placement,
	ReadId, Record, Strand,
};
use crate::candidate::{Candidate, settle, ungapped};
use crate::chain::Chain;
use crate::config::AlignConfig;
use crate::describe::{self, describe};
use crate::extend::{ExtensionBatch, ExtensionResults, Outcome};
use crate::index::Index;
use crate::kernel::{self, CpuScratch};
use crate::mapper::Mapper;
use crate::mapq;
use crate::nucleotide::encode;
use crate::num::{i32_of, i64_len, i64_of, u32_of, u64_of, usize_of};
use crate::pairing::{self, InsertModel, PairChoice};
use crate::rescue::Voter;

const MAX_EXTENSIONS: usize = 6;
const EXTEND_FRACTION: f32 = 0.5;
const EXTENSION_PAD: usize = 20;
const MAX_BAND: usize = 160;
const RESCUE_PAD: usize = 16;
const RESCUE_ANCHORS: usize = 3;
const MAX_RESCUE_WINDOW: i64 = 4096;
const MIN_NEW_BASES: u32 = 20;
const UNIQUE_MARGIN: i32 = 10;
const MAX_MATE_BOOST: u8 = 40;
const MIN_PAIRED_SCORE: i32 = 19;

#[derive(Debug, Clone, Copy)]
pub(crate) enum Layout<'m> {
	Single,
	Paired(&'m InsertModel),
}

#[derive(Debug, Clone, Copy)]
struct Query {
	start: usize,
	length: usize,
}

#[derive(Debug, Clone, Copy)]
struct Pending {
	owner: u32,
	strand: Strand,
	contig: ContigId,
	span_start: u64,
	window_start: i64,
}

#[derive(Debug, Clone, Copy)]
struct Decision {
	chosen: usize,
	mapq: u8,
	suboptimal: Option<i32>,
}

#[derive(Debug, Clone, Copy)]
struct PairContext {
	mate: Mate,
	proper: bool,
	template_length: i32,
}

#[derive(Debug, Clone, Copy)]
struct Emission {
	id: ReadId,
	flags: FlagSet,
	context: PairContext,
	mapq: u8,
	suboptimal: Option<i32>,
	primary: bool,
}

#[derive(Debug)]
pub(crate) struct Worker<'a> {
	index: &'a Index,
	config: AlignConfig,
	scratch: CpuScratch,
	mapper: Mapper<'a>,
	voter: Voter,
	chains: Vec<Chain>,
	codes: Vec<u8>,
	queries: Vec<Query>,
	repeats: Vec<f32>,
	candidates: Vec<Candidate>,
	ranges: Vec<Range<usize>>,
	operations: Vec<CigarOp>,
	batch: ExtensionBatch,
	pending: Vec<Pending>,
	results: ExtensionResults,
	window: Vec<u8>,
	cigar: Vec<CigarOp>,
	md: Vec<u8>,
	accepted: Vec<usize>,
}

impl<'a> Worker<'a> {
	pub(crate) fn new(index: &'a Index, config: AlignConfig) -> Self {
		Self {
			index,
			config,
			scratch: CpuScratch::default(),
			mapper: Mapper::new(index, &config),
			voter: Voter::default(),
			chains: Vec::new(),
			codes: Vec::new(),
			queries: Vec::new(),
			repeats: Vec::new(),
			candidates: Vec::new(),
			ranges: Vec::new(),
			operations: Vec::new(),
			batch: ExtensionBatch::default(),
			pending: Vec::new(),
			results: ExtensionResults::default(),
			window: Vec::new(),
			cigar: Vec::new(),
			md: Vec::new(),
			accepted: Vec::new(),
		}
	}

	pub(crate) fn align(
		&mut self,
		sequences: &[&[u8]],
		layout: Layout<'_>,
		first_read: u32,
		out: &mut AlignedBatch,
	) {
		self.seed_and_extend(sequences);
		if let Layout::Paired(model) = layout {
			self.rescue(model);
		}
		self.emit(layout, first_read, out);
	}

	pub(crate) fn distances(&mut self, sequences: &[&[u8]], out: &mut Vec<i64>) {
		self.seed_and_extend(sequences);
		for pair in 0..sequences.len() / 2 {
			let (Some(first), Some(second)) = (self.list(2 * pair), self.list(2 * pair + 1)) else {
				continue;
			};
			let (Some(left), Some(right)) = (first.first(), second.first()) else {
				continue;
			};
			let unique = |list: &[Candidate], best: &Candidate| {
				list.get(1)
					.is_none_or(|other| other.score + UNIQUE_MARGIN < best.score)
			};
			if left.score < self.config.min_score
				|| right.score < self.config.min_score
				|| !unique(first, left)
				|| !unique(second, right)
			{
				continue;
			}
			if let Some(span) = pairing::distance(left, right) {
				out.push(span);
			}
		}
	}

	fn seed_and_extend(&mut self, sequences: &[&[u8]]) {
		self.load(sequences);
		self.candidates.clear();
		self.operations.clear();
		self.batch.clear();
		self.pending.clear();
		self.generate();
		self.extend_pending();
		settle(&mut self.candidates, &mut self.ranges, self.queries.len());
	}

	fn load(&mut self, sequences: &[&[u8]]) {
		self.codes.clear();
		self.queries.clear();
		for sequence in sequences {
			let start = self.codes.len();
			encode(sequence, &mut self.codes);
			let length = sequence.len();
			for offset in (0..length).rev() {
				let code = self.codes.get(start + offset).copied().unwrap_or(4);
				self.codes.push(if code < 4 { 3 - code } else { code });
			}
			self.queries.push(Query { start, length });
		}
	}

	fn list(&self, end: usize) -> Option<&[Candidate]> {
		self.candidates.get(self.ranges.get(end)?.clone())
	}

	fn generate(&mut self) {
		self.repeats.clear();
		for end in 0..self.queries.len() {
			let Some(&query) = self.queries.get(end) else {
				break;
			};
			let codes = std::mem::take(&mut self.codes);
			let start = query.start;
			let forward = codes.get(start..start + query.length).unwrap_or_default();
			let reverse = codes
				.get(start + query.length..start + 2 * query.length)
				.unwrap_or_default();
			let mapped = self.mapper.chains(forward, reverse, &mut self.chains);
			self.repeats.push(mapped.repeat_fraction);
			let best = self.chains.first().map_or(0.0, |chain| chain.score);
			let chains = std::mem::take(&mut self.chains);
			for chain in chains
				.iter()
				.take(MAX_EXTENSIONS)
				.filter(|chain| chain.score >= EXTEND_FRACTION * best)
			{
				let codes_for = match chain.strand {
					Strand::Forward => forward,
					Strand::Reverse => reverse,
				};
				self.place_chain(u32_of(end), chain, codes_for);
			}
			self.chains = chains;
			self.codes = codes;
		}
	}

	fn place_chain(&mut self, owner: u32, chain: &Chain, query: &[u8]) {
		let reference = self.index.reference();
		let Some((contig, span)) = reference.locate(chain.reference_start) else {
			return;
		};
		let span_start = span.start;
		let start = chain.diagonal_middle;
		reference.window(span, start, query.len(), &mut self.window);
		if let Some(found) = ungapped(&self.config.scoring, query, &self.window) {
			let local = u64_of(start - i64_of(span_start));
			let operations =
				self.push_operations(&[CigarOp::new(CigarKind::Match, u32_of(query.len()))]);
			self.candidates.push(Candidate {
				owner,
				strand: chain.strand,
				contig,
				start: u32_of(usize_of(local)),
				stop: u32_of(usize_of(local)) + u32_of(query.len()),
				query_start: 0,
				query_stop: u32_of(query.len()),
				score: found.score,
				operations,
			});
			return;
		}
		let band = chain.diagonal_spread() + 2 * EXTENSION_PAD + 1;
		if band > MAX_BAND {
			return;
		}
		let window_start = chain.diagonal_low - i64_len(EXTENSION_PAD);
		reference.window(span, window_start, query.len() + band - 1, &mut self.window);
		self.batch.push(query, &self.window, band);
		self.pending.push(Pending {
			owner,
			strand: chain.strand,
			contig,
			span_start,
			window_start,
		});
	}

	fn push_operations(&mut self, operations: &[CigarOp]) -> Range<u32> {
		let start = u32_of(self.operations.len());
		self.operations.extend_from_slice(operations);
		start..u32_of(self.operations.len())
	}

	fn extend_pending(&mut self) {
		if self.batch.len() == 0 {
			return;
		}
		kernel::extend(
			&self.config.scoring,
			&mut self.scratch,
			&self.batch,
			&mut self.results,
		);
		let results = std::mem::take(&mut self.results);
		let pending = std::mem::take(&mut self.pending);
		for (outcome, job) in results.outcomes().iter().zip(&pending) {
			let Outcome::Aligned(found) = outcome else {
				continue;
			};
			if found.score < self.config.min_score.min(self.config.scoring.unpaired) {
				continue;
			}
			let flat = job.window_start + i64::from(found.window_start);
			let Ok(local) = u32::try_from(flat - i64_of(job.span_start)) else {
				continue;
			};
			let operations = self.push_operations(results.operations(found));
			self.candidates.push(Candidate {
				owner: job.owner,
				strand: job.strand,
				contig: job.contig,
				start: local,
				stop: local + (found.window_end - found.window_start),
				query_start: found.query_start,
				query_stop: found.query_end,
				score: found.score,
				operations,
			});
		}
		self.results = results;
		self.pending = pending;
	}

	fn rescue(&mut self, model: &InsertModel) {
		self.batch.clear();
		self.pending.clear();
		let (low, high) = model.window();
		for pair in 0..self.queries.len() / 2 {
			for anchor_end in [2 * pair, 2 * pair + 1] {
				let mate_end = anchor_end ^ 1;
				self.plan_rescue(anchor_end, mate_end, model, (low, high));
			}
		}
		if self.batch.len() == 0 {
			return;
		}
		self.extend_pending();
		settle(&mut self.candidates, &mut self.ranges, self.queries.len());
	}

	fn plan_rescue(
		&mut self,
		anchor_end: usize,
		mate_end: usize,
		model: &InsertModel,
		bounds: (i64, i64),
	) {
		let (Some(anchors), Some(mates), Some(&mate_query)) = (
			self.list(anchor_end),
			self.list(mate_end),
			self.queries.get(mate_end),
		) else {
			return;
		};
		let Some(best) = anchors.first().map(|candidate| candidate.score) else {
			return;
		};
		let wanted: Vec<Candidate> = anchors
			.iter()
			.take(RESCUE_ANCHORS)
			.filter(|candidate| candidate.score >= best - self.config.scoring.unpaired)
			.filter(|candidate| !pairing::compatible(model, candidate, mates))
			.cloned()
			.collect();
		let mate_length = i64_len(mate_query.length);
		let (low, high) = bounds;
		let reference = self.index.reference();
		for anchor in wanted {
			let Some(span) = reference.span(anchor.contig) else {
				continue;
			};
			let span_start = span.start;
			let (strand, local_start) = match anchor.strand {
				Strand::Forward => (Strand::Reverse, i64::from(anchor.start) + low - mate_length),
				Strand::Reverse => (Strand::Forward, i64::from(anchor.stop) - high),
			};
			let length = (high - low + mate_length).clamp(1, MAX_RESCUE_WINDOW);
			let flat_start = i64_of(span_start) + local_start;
			let mate_codes = strand_slice(&self.codes, mate_query, strand);
			reference.window(span, flat_start, usize_of(u64_of(length)), &mut self.window);
			let Some(diagonal) = self.voter.best_diagonal(mate_codes, &self.window) else {
				continue;
			};
			let band = 2 * RESCUE_PAD + 1;
			let window_start = flat_start + i64::from(diagonal) - i64_len(RESCUE_PAD);
			reference.window(
				span,
				window_start,
				mate_query.length + band - 1,
				&mut self.window,
			);
			let mate_codes = strand_slice(&self.codes, mate_query, strand);
			self.batch.push(mate_codes, &self.window, band);
			self.pending.push(Pending {
				owner: u32_of(mate_end),
				strand,
				contig: anchor.contig,
				span_start,
				window_start,
			});
		}
	}

	fn emit(&mut self, layout: Layout<'_>, first_read: u32, out: &mut AlignedBatch) {
		match layout {
			Layout::Single => {
				for end in 0..self.queries.len() {
					let decision = self.decide_single(end);
					let context = PairContext {
						mate: Mate::Unpaired,
						proper: false,
						template_length: 0,
					};
					let id = ReadId {
						index: first_read + u32_of(end),
						end: End::Single,
					};
					self.emit_end(out, id, end, decision, FlagSet::default(), context);
				}
			}
			Layout::Paired(model) => {
				for pair in 0..self.queries.len() / 2 {
					self.emit_pair(out, first_read + u32_of(pair), pair, model);
				}
			}
		}
	}

	fn decide_single(&self, end: usize) -> Option<Decision> {
		let list = self.list(end)?;
		let best = list.first()?;
		if best.score < self.config.min_score {
			return None;
		}
		Some(self.decision(end, list, 0))
	}

	fn decision(&self, end: usize, list: &[Candidate], chosen: usize) -> Decision {
		let length = self
			.queries
			.get(end)
			.map_or(0, |query| u32_of(query.length));
		let competition = mapq::competition(&self.config.scoring, list, chosen, length);
		let repeat = self.repeats.get(end).copied().unwrap_or(0.0);
		let mapq = list.get(chosen).map_or(0, |candidate| {
			mapq::single_end(&self.config.scoring, candidate, &competition, repeat)
		});
		Decision {
			chosen,
			mapq,
			suboptimal: competition.best_other,
		}
	}

	fn emit_pair(&mut self, out: &mut AlignedBatch, read: u32, pair: usize, model: &InsertModel) {
		let (first_end, second_end) = (2 * pair, 2 * pair + 1);
		let choice = match (self.list(first_end), self.list(second_end)) {
			(Some(first), Some(second)) => {
				pairing::choose_pair(&self.config.scoring, model, first, second)
					.filter(|found| pair_scores_valid(first, second, found))
			}
			_ => None,
		};
		let (first, second, proper) = match choice {
			Some(found) => {
				let (first, second) = self.decide_pair(first_end, second_end, &found);
				(Some(first), Some(second), true)
			}
			None => (
				self.decide_single(first_end),
				self.decide_single(second_end),
				false,
			),
		};
		let location = |end: usize, decision: Option<Decision>| {
			let candidate = decision.and_then(|found| self.list(end)?.get(found.chosen))?;
			Some((
				MateLocation {
					contig: candidate.contig,
					position: candidate.start,
					strand: candidate.strand,
				},
				candidate.start,
				candidate.stop,
			))
		};
		let first_location = location(first_end, first);
		let second_location = location(second_end, second);
		let template = match (&first_location, &second_location) {
			(Some((a, a_start, a_stop)), Some((b, b_start, b_stop))) if a.contig == b.contig => {
				let left = (*a_start).min(*b_start);
				let right = (*a_stop).max(*b_stop);
				let span = i32_of(i64::from(right) - i64::from(left));
				let first_is_left = (a_start, a.strand.index()) <= (b_start, b.strand.index());
				Some(if first_is_left { span } else { -span })
			}
			_ => None,
		};
		let contexts = [
			PairContext {
				mate: second_location.map_or(Mate::Unmapped, |found| Mate::Mapped(found.0)),
				proper,
				template_length: template.unwrap_or(0),
			},
			PairContext {
				mate: first_location.map_or(Mate::Unmapped, |found| Mate::Mapped(found.0)),
				proper,
				template_length: template.map_or(0, |value| -value),
			},
		];
		let ids = [
			ReadId {
				index: read,
				end: End::First,
			},
			ReadId {
				index: read,
				end: End::Second,
			},
		];
		let base = [
			FlagSet::default().with(Flag::Paired).with(Flag::First),
			FlagSet::default().with(Flag::Paired).with(Flag::Second),
		];
		self.emit_end(out, ids[0], first_end, first, base[0], contexts[0]);
		self.emit_end(out, ids[1], second_end, second, base[1], contexts[1]);
	}

	fn decide_pair(
		&self,
		first_end: usize,
		second_end: usize,
		choice: &PairChoice,
	) -> (Decision, Decision) {
		let first_list = self.list(first_end).unwrap_or_default();
		let second_list = self.list(second_end).unwrap_or_default();
		let repeat = self
			.repeats
			.get(first_end)
			.copied()
			.unwrap_or(0.0)
			.midpoint(self.repeats.get(second_end).copied().unwrap_or(0.0));
		let pair_mapq = mapq::paired(
			&self.config.scoring,
			choice.score,
			choice.runner_up,
			choice.competing,
			repeat,
		);
		let adjust = |decision: Decision| Decision {
			mapq: decision
				.mapq
				.max(pair_mapq)
				.min(pair_mapq.saturating_add(MAX_MATE_BOOST)),
			..decision
		};
		(
			adjust(self.decision(first_end, first_list, choice.first)),
			adjust(self.decision(second_end, second_list, choice.second)),
		)
	}

	fn emit_end(
		&mut self,
		out: &mut AlignedBatch,
		id: ReadId,
		end: usize,
		decision: Option<Decision>,
		base: FlagSet,
		context: PairContext,
	) {
		let mut flags = base;
		if matches!(context.mate, Mate::Unmapped) {
			flags = flags.with(Flag::MateUnmapped);
		}
		if let Mate::Mapped(mate) = context.mate
			&& mate.strand == Strand::Reverse
		{
			flags = flags.with(Flag::MateReverse);
		}
		if context.proper {
			flags = flags.with(Flag::ProperPair);
		}
		let Some(decision) = decision else {
			out.push_record(Record {
				read: id,
				placement: Placement::Unmapped {
					flags: flags.with(Flag::Unmapped),
					mate: context.mate,
				},
			});
			return;
		};
		let Some(range) = self.ranges.get(end).cloned() else {
			return;
		};
		let Some(&query) = self.queries.get(end) else {
			return;
		};
		let primary = Emission {
			id,
			flags,
			context,
			mapq: decision.mapq,
			suboptimal: decision.suboptimal,
			primary: true,
		};
		self.push_alignment(out, &primary, range.clone(), query, decision.chosen);
		self.accepted.clear();
		self.accepted.push(decision.chosen);
		let mut supplementary = 0usize;
		let mut secondary = 0usize;
		let length = u32_of(query.length);
		for index in 0..range.len() {
			if index == decision.chosen {
				continue;
			}
			let Some(list) = self.candidates.get(range.clone()) else {
				break;
			};
			let Some(candidate) = list.get(index) else {
				break;
			};
			let span = candidate.read_span(length);
			let overlap_with = |other: &Candidate| {
				let other_span = other.read_span(length);
				other_span
					.end
					.min(span.end)
					.saturating_sub(other_span.start.max(span.start))
			};
			let candidate_length = span.end - span.start;
			let primary_overlap = list.get(decision.chosen).map_or(0, overlap_with);
			let disjoint = self
				.accepted
				.iter()
				.filter_map(|&accepted| list.get(accepted))
				.all(|other| overlap_with(other) * 2 < candidate_length);
			let kind = if candidate.score >= self.config.min_score
				&& disjoint && candidate_length - primary_overlap >= MIN_NEW_BASES
				&& supplementary < self.config.max_supplementary
			{
				supplementary += 1;
				Some(Flag::Supplementary)
			} else if !disjoint && secondary < self.config.max_secondary {
				secondary += 1;
				Some(Flag::Secondary)
			} else {
				None
			};
			let Some(kind) = kind else {
				continue;
			};
			self.accepted.push(index);
			let competition = mapq::competition(&self.config.scoring, list, index, length);
			let repeat = self.repeats.get(end).copied().unwrap_or(0.0);
			let own = match kind {
				Flag::Secondary => 0,
				_ => mapq::single_end(&self.config.scoring, candidate, &competition, repeat),
			};
			let extra = Emission {
				flags: extra_flags(flags, kind),
				mapq: own,
				suboptimal: competition.best_other,
				primary: false,
				..primary
			};
			self.push_alignment(out, &extra, range.clone(), query, index);
		}
	}

	fn push_alignment(
		&mut self,
		out: &mut AlignedBatch,
		emission: &Emission,
		range: Range<usize>,
		query: Query,
		index: usize,
	) {
		let Emission {
			id,
			flags,
			context,
			mapq,
			suboptimal,
			primary,
		} = *emission;
		let Some(candidate) = self.candidates.get(range).and_then(|list| list.get(index)) else {
			return;
		};
		let reference = self.index.reference();
		let Some(span) = reference.span(candidate.contig) else {
			return;
		};
		reference.window(
			span,
			i64_of(span.start) + i64::from(candidate.start),
			candidate.reference_length() as usize,
			&mut self.window,
		);
		let codes = strand_slice(&self.codes, query, candidate.strand);
		let aligned = codes
			.get(candidate.query_start as usize..candidate.query_stop as usize)
			.unwrap_or_default();
		let operations = self
			.operations
			.get(candidate.operations.start as usize..candidate.operations.end as usize)
			.unwrap_or_default();
		let edit_distance = describe(
			&describe::Alignment {
				operations,
				query: aligned,
				reference: &self.window,
				left_clip: candidate.query_start,
				right_clip: u32_of(query.length).saturating_sub(candidate.query_stop),
			},
			&mut self.cigar,
			&mut self.md,
		);
		let flags = match candidate.strand {
			Strand::Reverse => flags.with(Flag::Reverse),
			Strand::Forward => flags,
		};
		let cigar = out.push_cigar(&self.cigar);
		let md = out.push_md(&self.md);
		out.push_record(Record {
			read: id,
			placement: Placement::Mapped(Alignment {
				contig: candidate.contig,
				position: candidate.start,
				strand: candidate.strand,
				flags,
				mapq,
				score: candidate.score,
				suboptimal_score: suboptimal,
				edit_distance,
				mate: context.mate,
				template_length: if primary { context.template_length } else { 0 },
				cigar,
				md,
			}),
		});
	}
}

fn extra_flags(flags: FlagSet, kind: Flag) -> FlagSet {
	match kind {
		Flag::Secondary => flags.without(Flag::ProperPair).with(kind),
		_ => flags.with(kind),
	}
}

fn strand_slice(codes: &[u8], query: Query, strand: Strand) -> &[u8] {
	let from = query.start + strand.index() * query.length;
	codes.get(from..from + query.length).unwrap_or_default()
}

fn pair_scores_valid(first: &[Candidate], second: &[Candidate], choice: &PairChoice) -> bool {
	first
		.get(choice.first)
		.is_some_and(|candidate| candidate.score >= MIN_PAIRED_SCORE)
		&& second
			.get(choice.second)
			.is_some_and(|candidate| candidate.score >= MIN_PAIRED_SCORE)
}
