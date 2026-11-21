use std::hint::black_box;

use crate::alignment::Strand;
use crate::chain::{Chain, Chainer, pack_anchor};
use crate::config::AlignConfig;
use crate::index::{Index, Probe, entry_offset, entry_position};
use crate::seeding::{Seed, Seeder};

const MAX_FIRST_STROBE_HITS: usize = 20;
const QUICK_SEEDS: usize = 6;
const QUICK_HITS: usize = 4;

#[derive(Debug, Clone, Copy)]
pub(crate) struct Mapped {
	pub(crate) repeat_fraction: f32,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Diagonal {
	pub(crate) strand: Strand,
	pub(crate) start: i64,
	pub(crate) votes: u32,
}

#[derive(Debug, Clone, Copy)]
struct SeedHits {
	first: usize,
	relaxed: usize,
}

#[derive(Debug)]
pub(crate) struct Mapper<'a> {
	index: &'a Index,
	hits: SeedHits,
	seeder: Seeder,
	seeds: [Vec<Seed>; 2],
	probes: Vec<Probe>,
	quick: Vec<(Strand, u32, Probe)>,
	misses: Vec<usize>,
	anchors: Vec<u64>,
	chainer: Chainer,
}

impl<'a> Mapper<'a> {
	pub(crate) fn new(index: &'a Index, config: &AlignConfig) -> Self {
		Self {
			index,
			hits: SeedHits {
				first: config.max_seed_hits.max(index.filter_cutoff()),
				relaxed: config.relaxed_seed_hits,
			},
			seeder: Seeder::new(index.config()),
			seeds: [Vec::new(), Vec::new()],
			probes: Vec::new(),
			quick: Vec::new(),
			misses: Vec::new(),
			anchors: Vec::new(),
			chainer: Chainer::default(),
		}
	}

	pub(crate) fn chains(&mut self, forward: &[u8], chains: &mut Vec<Chain>) -> Mapped {
		self.seeder.extract_both(forward, &mut self.seeds);
		let mut repeat_fraction = self.pass(self.hits.first, chains);
		if chains.is_empty() && repeat_fraction > 0.0 {
			repeat_fraction = self.pass(self.hits.relaxed, chains);
		}
		chains.sort_unstable_by(|left, right| {
			right
				.score
				.total_cmp(&left.score)
				.then(left.strand.index().cmp(&right.strand.index()))
				.then(left.reference_start.cmp(&right.reference_start))
		});
		Mapped { repeat_fraction }
	}

	pub(crate) fn locate(&mut self, forward: &[u8], found: &mut Vec<Diagonal>) {
		found.clear();
		self.seeder.extract_both(forward, &mut self.seeds);
		self.quick.clear();
		let mut touched = 0u64;
		for (strand, seeds) in [Strand::Forward, Strand::Reverse]
			.into_iter()
			.zip(&self.seeds)
		{
			let wanted = QUICK_SEEDS.min(seeds.len());
			for step in 0..wanted {
				let seed = &seeds[(2 * step + 1) * seeds.len() / (2 * wanted)];
				let probe = self.index.probe(seed);
				touched ^= self.index.touch(&probe);
				self.quick.push((strand, seed.position, probe));
			}
		}
		black_box(touched);
		for (strand, position, probe) in &self.quick {
			let hits = self.index.run(probe);
			if hits.len() > QUICK_HITS {
				continue;
			}
			for &hit in hits {
				let start = i64::from(entry_position(hit)) - i64::from(*position);
				match found
					.iter_mut()
					.find(|found| found.strand == *strand && found.start == start)
				{
					Some(existing) => existing.votes += 1,
					None => found.push(Diagonal {
						strand: *strand,
						start,
						votes: 1,
					}),
				}
			}
		}
		found.sort_unstable_by(|left, right| {
			right
				.votes
				.cmp(&left.votes)
				.then(left.strand.index().cmp(&right.strand.index()))
				.then(left.start.cmp(&right.start))
		});
	}

	fn pass(&mut self, cutoff: usize, chains: &mut Vec<Chain>) -> f32 {
		chains.clear();
		let mut repeat_fraction = 0.0f32;
		for strand in [Strand::Forward, Strand::Reverse] {
			repeat_fraction = repeat_fraction.max(self.chains_for(strand, cutoff, chains));
		}
		repeat_fraction
	}

	fn chains_for(&mut self, strand: Strand, cutoff: usize, chains: &mut Vec<Chain>) -> f32 {
		let seeds = std::mem::take(&mut self.seeds[strand.index()]);
		let seed_length = self.index.config().k;
		let mut probes = std::mem::take(&mut self.probes);
		probes.clear();
		let mut touched = 0u64;
		for seed in &seeds {
			let probe = self.index.probe(seed);
			touched ^= self.index.touch(&probe);
			probes.push(probe);
		}
		black_box(touched);
		self.anchors.clear();
		self.misses.clear();
		let mut skipped = 0usize;
		for (number, (seed, probe)) in seeds.iter().zip(&probes).enumerate() {
			let hits = self.index.run(probe);
			if hits.len() > cutoff {
				skipped += 1;
				continue;
			}
			if hits.is_empty() {
				self.misses.push(number);
				continue;
			}
			for &hit in hits {
				let position = u64::from(entry_position(hit));
				self.anchors.push(pack_anchor(position, seed.position));
				let offset = entry_offset(hit);
				if offset != 0 && seed.offset != 0 {
					self.anchors.push(pack_anchor(
						position + u64::from(offset),
						seed.position + seed.offset,
					));
				}
			}
		}
		for &number in &self.misses {
			let (Some(seed), Some(probe)) = (seeds.get(number), probes.get(number)) else {
				continue;
			};
			let first = self.anchors.len();
			for hit in self
				.index
				.first_strobe_run(probe)
				.take(MAX_FIRST_STROBE_HITS + 1)
			{
				self.anchors
					.push(pack_anchor(u64::from(entry_position(hit)), seed.position));
			}
			if self.anchors.len() - first > MAX_FIRST_STROBE_HITS {
				self.anchors.truncate(first);
			}
		}
		let chain_length = u32::try_from(seed_length).unwrap_or(u32::MAX);
		self.chainer.chain(
			&mut self.anchors,
			strand,
			chain_length,
			self.index.reference().spans(),
			chains,
		);
		let total = seeds.len();
		self.seeds[strand.index()] = seeds;
		self.probes = probes;
		if total == 0 {
			0.0
		} else {
			let total = u16::try_from(total).unwrap_or(u16::MAX);
			f32::from(u16::try_from(skipped).unwrap_or(u16::MAX)) / f32::from(total)
		}
	}
}
