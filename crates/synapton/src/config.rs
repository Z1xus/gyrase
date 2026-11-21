use crate::Error;

pub(crate) const MAX_STROBE_OFFSET: usize = 127;
pub(crate) const MAX_SYNCMER_WINDOW: usize = 16;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeedConfig {
	pub(crate) k: usize,
	pub(crate) s: usize,
	pub(crate) w_min: usize,
	pub(crate) w_max: usize,
}

impl SeedConfig {
	pub const fn short_read() -> Self {
		Self {
			k: 20,
			s: 16,
			w_min: 5,
			w_max: 11,
		}
	}

	pub fn new(k: usize, s: usize, w_min: usize, w_max: usize) -> Result<Self, Error> {
		if !(8..=32).contains(&k) {
			return Err(Error::InvalidSeeding("k must be between 8 and 32"));
		}
		if s == 0 || s > k {
			return Err(Error::InvalidSeeding("s must be between 1 and k"));
		}
		if !(k - s).is_multiple_of(2) {
			return Err(Error::InvalidSeeding("k minus s must be even"));
		}
		if k - s + 1 > MAX_SYNCMER_WINDOW {
			return Err(Error::InvalidSeeding("k minus s is too large"));
		}
		if w_min == 0 || w_min > w_max || w_max > 64 {
			return Err(Error::InvalidSeeding(
				"window must satisfy 1 <= min <= max <= 64",
			));
		}
		Ok(Self { k, s, w_min, w_max })
	}

	pub(crate) const fn syncmer_window(&self) -> usize {
		self.k - self.s + 1
	}

	pub(crate) const fn syncmer_target(&self) -> usize {
		(self.k - self.s) / 2
	}

	pub(crate) const fn expected_spacing(&self) -> usize {
		self.syncmer_window()
	}
}

impl Default for SeedConfig {
	fn default() -> Self {
		Self::short_read()
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Scoring {
	pub matched: i16,
	pub mismatch: i16,
	pub gap_open: i16,
	pub gap_extend: i16,
	pub ambiguous: i16,
	pub clip: i16,
	pub unpaired: i32,
}

impl Default for Scoring {
	fn default() -> Self {
		Self {
			matched: 1,
			mismatch: 4,
			gap_open: 6,
			gap_extend: 1,
			ambiguous: 1,
			clip: 5,
			unpaired: 17,
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignConfig {
	pub scoring: Scoring,
	pub min_score: i32,
	pub max_secondary: usize,
	pub max_supplementary: usize,
	pub max_seed_hits: usize,
	pub relaxed_seed_hits: usize,
	pub calibration_pairs: usize,
	pub chunk_reads: usize,
	pub md: bool,
}

impl Default for AlignConfig {
	fn default() -> Self {
		Self {
			scoring: Scoring::default(),
			min_score: 30,
			max_secondary: 0,
			max_supplementary: 2,
			max_seed_hits: 300,
			relaxed_seed_hits: 1000,
			calibration_pairs: 8192,
			chunk_reads: 512,
			md: false,
		}
	}
}
