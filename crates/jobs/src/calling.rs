use std::collections::BTreeMap;
use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;
use std::thread;

use gyrase_call::{Calibration, Call, Caller, IndelRates, Ploidy, Site};
use gyrase_extract::{Sex, SexCall, SexEvidence, infer_sex};
use gyrase_genome::{Base, Build, Chromosome, ContigId, Interval};
use gyrase_io::Workers;
use gyrase_pileup::{Filter, Pileup, Scratch, Window};
use rayon::prelude::*;

use crate::alignments::EMPTY_COST;
use crate::{Alignments, Cancel, Error};

const SHARD_LEN: u64 = 1 << 18;
const MIN_SHARD_LEN: u64 = 1 << 14;
const SPLIT_FACTOR: u64 = 2;
const PAD: u64 = 2_000;
const SURVEY_STRIDE: u64 = 16;
const COVERAGE_MAPPING_QUALITY: u8 = 20;

#[derive(Debug, Clone)]
pub struct Shard {
	pub interval: Interval,
	pub chromosome: Chromosome,
	pub cost: u64,
	pseudoautosomal: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Survey {
	pub sex: SexCall,
	pub rates: IndelRates,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Target {
	pub contig: ContigId,
	pub position: u64,
	pub alternate: Option<Base>,
}

#[derive(Debug, Default, Clone, Copy)]
struct Coverage {
	autosome: Depth,
	x: Depth,
	y: Depth,
}

#[derive(Debug, Default, Clone, Copy)]
struct Depth {
	bases: u64,
	positions: u64,
}

#[derive(Debug, Clone, Copy)]
pub struct Calling<'a> {
	pub alignments: &'a Alignments,
	pub workers: &'a Workers,
	pub io: &'a Workers,
	pub cancel: &'a Cancel,
}

#[derive(Debug, Default)]
struct Calls {
	discovered: Vec<Call>,
	known: Vec<Call>,
}

pub fn shards(alignments: &Alignments, build: Build) -> Vec<Shard> {
	let pseudoautosomal = build.pseudoautosomal();
	let mut shards = Vec::new();
	for (index, reference) in alignments.header().references().iter().enumerate() {
		let Some(chromosome) = Chromosome::parse(reference.name()) else {
			continue;
		};
		let contig = ContigId::new(u32::try_from(index).unwrap_or(u32::MAX));
		if alignments.reference_contig(contig).is_none() {
			continue;
		}
		let cuts: Vec<u64> = if chromosome == Chromosome::X {
			pseudoautosomal
				.iter()
				.flat_map(|range| [range.start, range.end])
				.collect()
		} else {
			Vec::new()
		};
		let length = u64::from(reference.length());
		let mut start = 0;
		while start < length {
			let end = cuts
				.iter()
				.copied()
				.filter(|&cut| cut > start)
				.fold((start + SHARD_LEN).min(length), u64::min);
			let interval = Interval { contig, start, end };
			shards.push(Shard {
				interval,
				chromosome,
				cost: alignments.cost(interval),
				pseudoautosomal: chromosome == Chromosome::X
					&& pseudoautosomal.iter().any(|range| range.contains(&start)),
			});
			start = end;
		}
	}
	let mut costs: Vec<u64> = shards
		.iter()
		.map(|shard| shard.cost)
		.filter(|&cost| cost > EMPTY_COST)
		.collect();
	costs.sort_unstable();
	let limit = costs
		.get(costs.len() / 2)
		.map_or(u64::MAX, |median| median.saturating_mul(SPLIT_FACTOR));
	shards
		.into_iter()
		.flat_map(|shard| split(alignments, shard, limit))
		.collect()
}

fn split(alignments: &Alignments, shard: Shard, limit: u64) -> Vec<Shard> {
	let interval = shard.interval;
	if shard.cost <= limit || interval.len() <= MIN_SHARD_LEN {
		return vec![shard];
	}
	let middle = interval.start + interval.len() / 2;
	[(interval.start, middle), (middle, interval.end)]
		.into_iter()
		.flat_map(|(start, end)| {
			let interval = Interval {
				start,
				end,
				..interval
			};
			let half = Shard {
				interval,
				cost: alignments.cost(interval),
				..shard.clone()
			};
			split(alignments, half, limit)
		})
		.collect()
}

impl Shard {
	pub fn ploidy(&self, sex: SexCall) -> Option<Ploidy> {
		let male = sex == SexCall::Determined(Sex::Male);
		match self.chromosome {
			Chromosome::X if male && !self.pseudoautosomal => Some(Ploidy::Haploid),
			Chromosome::Y if sex == SexCall::Determined(Sex::Female) => None,
			Chromosome::Autosome(_) | Chromosome::X => Some(Ploidy::Diploid),
			Chromosome::Y | Chromosome::Mitochondrion => Some(Ploidy::Haploid),
		}
	}
}

impl Calling<'_> {
	pub fn survey(self, shards: &[Shard]) -> Result<Survey, Error> {
		let parts: Vec<Result<(Calibration, Coverage), Error>> = self.workers.install(|| {
			shards
				.par_iter()
				.filter(|shard| (shard.interval.start / SHARD_LEN).is_multiple_of(SURVEY_STRIDE))
				.map_init(Scratch::new, |scratch, shard| {
					self.cancel.check()?;
					self.with_pileup(shard, scratch, |pileup| {
						let mut calibration = Calibration::new();
						if matches!(shard.chromosome, Chromosome::Autosome(_)) {
							calibration.observe(pileup);
						}
						(calibration, coverage(shard, pileup))
					})
				})
				.collect()
		});
		let mut calibration = Calibration::new();
		let mut coverage = Coverage::default();
		for part in parts {
			let (part_calibration, part_coverage) = part?;
			calibration.merge(&part_calibration);
			coverage.merge(part_coverage);
		}
		Ok(Survey {
			sex: infer_sex(SexEvidence {
				x_depth_ratio: coverage.x.mean() / coverage.autosome.mean(),
				y_depth_ratio: coverage.y.mean() / coverage.autosome.mean(),
				x_panel_heterozygous_fraction: f64::NAN,
			}),
			rates: calibration.rates(),
		})
	}

	pub fn call(
		self,
		shards: &[Shard],
		survey: &Survey,
		targets: &[Target],
		mut emit: impl FnMut(&Shard, &[Call]) -> Result<(), Error>,
	) -> Result<Vec<Option<Call>>, Error> {
		let next = AtomicUsize::new(0);
		let (sender, receiver) = mpsc::channel();
		thread::scope(|scope| {
			scope.spawn(|| {
				self.workers.install(|| {
					rayon::broadcast(|_| {
						let mut scratch = Scratch::new();
						loop {
							let index = next.fetch_add(1, Ordering::Relaxed);
							let Some(shard) = shards.get(index) else {
								break;
							};
							let calls = self.shard(shard, survey, targets, &mut scratch);
							if sender.send((index, calls)).is_err() {
								break;
							}
						}
					});
				});
				drop(sender);
			});
			let mut genotyped: Vec<Option<Call>> = vec![None; targets.len()];
			let mut pending = BTreeMap::new();
			let mut ready = 0;
			for (index, calls) in receiver {
				pending.insert(index, calls);
				while let Some(calls) = pending.remove(&ready) {
					let Calls { discovered, known } = calls?;
					let shard = shards.get(ready).ok_or(Error::UnknownContig)?;
					emit(shard, &discovered)?;
					for (slot, call) in genotyped[target_range(targets, shard)]
						.iter_mut()
						.zip(known)
					{
						*slot = Some(call);
					}
					ready += 1;
				}
			}
			Ok(genotyped)
		})
	}

	fn shard(
		self,
		shard: &Shard,
		survey: &Survey,
		targets: &[Target],
		scratch: &mut Scratch,
	) -> Result<Calls, Error> {
		self.cancel.check()?;
		let Some(ploidy) = shard.ploidy(survey.sex) else {
			return Ok(Calls::default());
		};
		if shard.cost <= EMPTY_COST {
			return Ok(Calls::default());
		}
		let mut caller = Caller::new(ploidy);
		caller.indel_rates.clone_from(&survey.rates);
		let targets = &targets[target_range(targets, shard)];
		self.with_pileup(shard, scratch, |pileup| {
			let sites = sites(pileup.window, targets);
			let mut calls = Calls::default();
			caller.discover(pileup, &mut calls.discovered);
			caller.genotype(pileup, &sites, &mut calls.known);
			calls
		})
	}

	fn with_pileup<T>(
		self,
		shard: &Shard,
		scratch: &mut Scratch,
		work: impl FnOnce(&Pileup<'_>) -> T,
	) -> Result<T, Error> {
		let alignments = self.alignments;
		let interval = shard.interval;
		let contig = alignments
			.reference_contig(interval.contig)
			.ok_or(Error::UnknownContig)?;
		let length = alignments
			.reference()
			.contigs()
			.get(contig)
			.map_or(interval.end, |reference| reference.length);
		let window = Interval {
			contig,
			start: interval.start.saturating_sub(PAD),
			end: (interval.end + PAD).min(length),
		};
		let mut bases = Vec::new();
		alignments.reference().fetch(&window, &mut bases)?;
		let batches = alignments.fetch(interval, self.io)?;
		let pileup = scratch.pileup(
			interval,
			Window {
				start: window.start,
				bases: &bases,
			},
			Filter::default(),
			batches.records(),
		)?;
		Ok(work(&pileup))
	}
}

fn target_range(targets: &[Target], shard: &Shard) -> Range<usize> {
	let key = |contig: ContigId, position: u64| (contig.get(), position);
	let interval = shard.interval;
	let begin = targets.partition_point(|target| {
		key(target.contig, target.position) < key(interval.contig, interval.start)
	});
	let end = targets.partition_point(|target| {
		key(target.contig, target.position) < key(interval.contig, interval.end)
	});
	begin..end.max(begin)
}

fn sites(window: Window<'_>, targets: &[Target]) -> Vec<Site> {
	targets
		.iter()
		.map(|target| Site {
			position: target.position,
			reference: usize::try_from(target.position - window.start)
				.ok()
				.and_then(|offset| window.bases.get(offset..=offset))
				.unwrap_or(b"N")
				.to_vec(),
			alternates: target
				.alternate
				.map(|base| vec![vec![base.ascii()]])
				.unwrap_or_default(),
		})
		.collect()
}

fn coverage(shard: &Shard, pileup: &Pileup<'_>) -> Coverage {
	let mut depth = Depth::default();
	for column in pileup.columns() {
		depth.positions += 1;
		depth.bases += column
			.observations
			.iter()
			.filter(|observation| observation.mapping_quality >= COVERAGE_MAPPING_QUALITY)
			.count() as u64;
	}
	let mut coverage = Coverage::default();
	match shard.chromosome {
		Chromosome::Autosome(_) => coverage.autosome = depth,
		Chromosome::X if !shard.pseudoautosomal => coverage.x = depth,
		Chromosome::Y => coverage.y = depth,
		Chromosome::X | Chromosome::Mitochondrion => {}
	}
	coverage
}

impl Coverage {
	fn merge(&mut self, other: Self) {
		for (mine, theirs) in [
			(&mut self.autosome, other.autosome),
			(&mut self.x, other.x),
			(&mut self.y, other.y),
		] {
			mine.bases += theirs.bases;
			mine.positions += theirs.positions;
		}
	}
}

impl Depth {
	#[allow(clippy::cast_precision_loss)]
	fn mean(self) -> f64 {
		self.bases as f64 / self.positions as f64
	}
}
