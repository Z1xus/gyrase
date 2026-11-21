use crate::Error;

const MAX_MIN_SHIFT: u32 = 31;
const MAX_DEPTH: u32 = 9;
const LEVEL_BITS: u32 = 3;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Geometry {
	min_shift: u32,
	depth: u32,
}

const fn first_bin(level: u32) -> u32 {
	((1 << (LEVEL_BITS * level)) - 1) / 7
}

impl Geometry {
	pub const BAI: Self = Self {
		min_shift: 14,
		depth: 5,
	};

	pub fn new(min_shift: u32, depth: u32) -> Result<Self, Error> {
		if (1..=MAX_MIN_SHIFT).contains(&min_shift) && (0..=MAX_DEPTH).contains(&depth) {
			Ok(Self { min_shift, depth })
		} else {
			Err(Error::InvalidGeometry { min_shift, depth })
		}
	}

	pub fn covering(min_shift: u32, longest: u64) -> Result<Self, Error> {
		let invalid = |depth| Error::InvalidGeometry { min_shift, depth };
		if !(1..=MAX_MIN_SHIFT).contains(&min_shift) {
			return Err(invalid(0));
		}
		let depth = (0..=MAX_DEPTH)
			.find(|depth| longest <= 1 << (min_shift + LEVEL_BITS * depth))
			.ok_or_else(|| invalid(MAX_DEPTH + 1))?;
		Self::new(min_shift, depth)
	}

	pub const fn min_shift(self) -> u32 {
		self.min_shift
	}

	pub const fn depth(self) -> u32 {
		self.depth
	}

	pub const fn max_position(self) -> u64 {
		1 << (self.min_shift + LEVEL_BITS * self.depth)
	}

	pub const fn bin_count(self) -> u32 {
		first_bin(self.depth + 1)
	}

	pub const fn metadata_bin(self) -> u32 {
		self.bin_count() + 1
	}

	pub(crate) const fn first_leaf(self) -> u32 {
		first_bin(self.depth)
	}

	const fn level_shift(self, level: u32) -> u32 {
		self.min_shift + LEVEL_BITS * (self.depth - level)
	}

	pub fn bin(self, start: u64, end: u64) -> u32 {
		let last = end.wrapping_sub(1);
		for level in (1..=self.depth).rev() {
			let shift = self.level_shift(level);
			if start >> shift == last >> shift {
				return u32::try_from(start >> shift).map_or(0, |offset| first_bin(level) + offset);
			}
		}
		0
	}

	pub fn bins(self, start: u64, end: u64, out: &mut Vec<u32>) {
		out.clear();
		let end = end.min(self.max_position());
		if start >= end {
			return;
		}
		let last = end - 1;
		for level in 0..=self.depth {
			let shift = self.level_shift(level);
			let (Ok(low), Ok(high)) = (u32::try_from(start >> shift), u32::try_from(last >> shift))
			else {
				return;
			};
			let first = first_bin(level);
			out.extend(first + low..=first + high);
		}
	}

	pub(crate) fn level(self, bin: u32) -> Option<u32> {
		(0..=self.depth).find(|&level| bin < first_bin(level + 1))
	}

	pub(crate) fn bottom_leaf(self, bin: u32) -> Option<u64> {
		let level = self.level(bin)?;
		Some(u64::from(bin - first_bin(level)) << (LEVEL_BITS * (self.depth - level)))
	}
}
