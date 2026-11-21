use crate::alignment::Strand;
use crate::chain::{Chain, Chainer, pack_anchor};
use crate::config::AlignConfig;
use crate::index::{Index, Probe, entry_offset, entry_position};
use crate::seeding::{Seed, Seeder};

const MAX_FIRST_STROBE_HITS: usize = 20;

#[derive(Debug, Clone, Copy)]
pub(crate) struct Mapped {
	pub(crate) repeat_fraction: f32,
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
	seeds: Vec<Seed>,
	probes: Vec<Probe>,
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
			seeds: Vec::new(),
			probes: Vec::new(),
			anchors: Vec::new(),
			chainer: Chainer::default(),
		}
	}

	pub(crate) fn chains(
		&mut self,
		forward: &[u8],
		reverse: &[u8],
		chains: &mut Vec<Chain>,
	) -> Mapped {
		let mut repeat_fraction = self.pass(forward, reverse, self.hits.first, chains);
		if chains.is_empty() && repeat_fraction > 0.0 {
			repeat_fraction = self.pass(forward, reverse, self.hits.relaxed, chains);
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

	fn pass(
		&mut self,
		forward: &[u8],
		reverse: &[u8],
		cutoff: usize,
		chains: &mut Vec<Chain>,
	) -> f32 {
		chains.clear();
		let mut repeat_fraction = 0.0f32;
		for strand in [Strand::Forward, Strand::Reverse] {
			let codes = match strand {
				Strand::Forward => forward,
				Strand::Reverse => reverse,
			};
			repeat_fraction = repeat_fraction.max(self.chains_for(codes, strand, cutoff, chains));
		}
		repeat_fraction
	}

	fn chains_for(
		&mut self,
		codes: &[u8],
		strand: Strand,
		cutoff: usize,
		chains: &mut Vec<Chain>,
	) -> f32 {
		self.seeds.clear();
		self.seeder.extract(codes, u32::MAX, &mut self.seeds);
		self.probes.clear();
		self.probes
			.extend(self.seeds.iter().map(|seed| self.index.probe(seed)));
		self.anchors.clear();
		let mut skipped = 0usize;
		for (seed, probe) in self.seeds.iter().zip(&self.probes) {
			let hits = self.index.run(probe);
			if hits.len() > cutoff {
				skipped += 1;
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
			if hits.is_empty()
				&& self.index.first_strobe_run(probe).count() <= MAX_FIRST_STROBE_HITS
			{
				for hit in self.index.first_strobe_run(probe) {
					let position = u64::from(entry_position(hit));
					self.anchors.push(pack_anchor(position, seed.position));
				}
			}
		}
		let seed_length = u32::try_from(self.index.config().k).unwrap_or(u32::MAX);
		self.chainer.chain(
			&mut self.anchors,
			strand,
			seed_length,
			self.index.reference().spans(),
			chains,
		);
		if self.seeds.is_empty() {
			0.0
		} else {
			let total = u16::try_from(self.seeds.len()).unwrap_or(u16::MAX);
			f32::from(u16::try_from(skipped).unwrap_or(u16::MAX)) / f32::from(total)
		}
	}
}
