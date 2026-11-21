use crate::config::{MAX_STROBE_OFFSET, MAX_SYNCMER_WINDOW, SeedConfig};
use crate::hash::mix;
use crate::num::u32_of;

pub(crate) const AUX_BITS: u32 = 12;
const AUX_MASK: u64 = (1 << AUX_BITS) - 1;
const STROBE_SELECTION_MASK: u64 = 0xff;
const ORIENTATION_FIRST: u8 = 0b10;
const ORIENTATION_SECOND: u8 = 0b01;

#[derive(Debug, Clone, Copy)]
struct Syncmer {
	position: u32,
	hash: u64,
	forward: bool,
	reverse: bool,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Seed {
	pub(crate) position: u32,
	pub(crate) offset: u32,
	pub(crate) key: u64,
	pub(crate) orientation: u8,
}

#[derive(Debug)]
pub(crate) struct Seeder {
	config: SeedConfig,
	syncmers: Vec<Syncmer>,
	flipped: Vec<Syncmer>,
}

impl Seeder {
	pub(crate) fn new(config: SeedConfig) -> Self {
		Self {
			config,
			syncmers: Vec::new(),
			flipped: Vec::new(),
		}
	}

	pub(crate) fn extract(&mut self, codes: &[u8], limit: u32, seeds: &mut Vec<Seed>) {
		self.collect_syncmers(codes);
		collect_seeds(&self.config, &self.syncmers, limit, seeds);
	}

	pub(crate) fn syncmer_positions(
		&mut self,
		codes: &[u8],
		limit: u32,
	) -> impl Iterator<Item = u32> + '_ {
		self.collect_syncmers(codes);
		self.syncmers
			.iter()
			.take_while(move |syncmer| syncmer.position < limit)
			.map(|syncmer| syncmer.position)
	}

	pub(crate) fn extract_both(&mut self, codes: &[u8], seeds: &mut [Vec<Seed>; 2]) {
		self.collect_syncmers(codes);
		self.flipped.clear();
		let k = u32_of(self.config.k);
		let length = u32_of(codes.len());
		self.flipped
			.extend(self.syncmers.iter().rev().map(|syncmer| Syncmer {
				position: length - k - syncmer.position,
				hash: syncmer.hash,
				forward: syncmer.reverse,
				reverse: syncmer.forward,
			}));
		let [forward, reverse] = seeds;
		forward.clear();
		reverse.clear();
		collect_seeds(&self.config, &self.syncmers, u32::MAX, forward);
		collect_seeds(&self.config, &self.flipped, u32::MAX, reverse);
	}

	fn collect_syncmers(&mut self, codes: &[u8]) {
		self.syncmers.clear();
		let SeedConfig { k, s, .. } = self.config;
		let window = self.config.syncmer_window();
		let target = self.config.syncmer_target();
		let k_mask = mask(k);
		let s_mask = mask(s);
		let k_shift = 2 * (k - 1);
		let s_shift = 2 * (s - 1);
		let mut ring = [0u64; MAX_SYNCMER_WINDOW];
		let (mut forward_k, mut reverse_k, mut forward_s, mut reverse_s) = (0u64, 0u64, 0u64, 0u64);
		let mut valid = 0usize;
		let mut smers = 0usize;
		for (index, &code) in codes.iter().enumerate() {
			if code > 3 {
				valid = 0;
				smers = 0;
				(forward_k, reverse_k, forward_s, reverse_s) = (0, 0, 0, 0);
				continue;
			}
			let base = u64::from(code);
			let complement = 3 - base;
			forward_k = ((forward_k << 2) | base) & k_mask;
			reverse_k = (reverse_k >> 2) | (complement << k_shift);
			forward_s = ((forward_s << 2) | base) & s_mask;
			reverse_s = (reverse_s >> 2) | (complement << s_shift);
			valid += 1;
			if valid < s {
				continue;
			}
			if let Some(slot) = ring.get_mut(smers % window) {
				*slot = mix(forward_s.min(reverse_s));
			}
			smers += 1;
			if valid < k {
				continue;
			}
			let smallest = ring.iter().take(window).copied().min().unwrap_or(u64::MAX);
			let centre = ring
				.get((smers - window + target) % window)
				.copied()
				.unwrap_or(u64::MAX);
			if centre == smallest {
				self.syncmers.push(Syncmer {
					position: u32_of(index + 1 - k),
					hash: mix(forward_k.min(reverse_k)),
					forward: forward_k <= reverse_k,
					reverse: reverse_k <= forward_k,
				});
			}
		}
	}
}

fn collect_seeds(config: &SeedConfig, syncmers: &[Syncmer], limit: u32, seeds: &mut Vec<Seed>) {
	let SeedConfig { w_min, w_max, .. } = *config;
	for (index, first) in syncmers.iter().enumerate() {
		if first.position >= limit {
			break;
		}
		let reach = first.position + u32_of(MAX_STROBE_OFFSET);
		let candidates = syncmers
			.iter()
			.skip(index + w_min)
			.take(w_max - w_min + 1)
			.take_while(|candidate| candidate.position <= reach);
		let mut second = first;
		let mut fewest = u32::MAX;
		for candidate in candidates {
			let ones = ((first.hash ^ candidate.hash) & STROBE_SELECTION_MASK).count_ones();
			if ones < fewest {
				fewest = ones;
				second = candidate;
			}
		}
		seeds.push(Seed {
			position: first.position,
			offset: second.position - first.position,
			key: (first.hash & !AUX_MASK) | (second.hash & AUX_MASK),
			orientation: if first.forward { ORIENTATION_FIRST } else { 0 }
				| if second.forward {
					ORIENTATION_SECOND
				} else {
					0
				},
		});
	}
}

const fn mask(bases: usize) -> u64 {
	if bases >= 32 {
		u64::MAX
	} else {
		(1u64 << (2 * bases)) - 1
	}
}
