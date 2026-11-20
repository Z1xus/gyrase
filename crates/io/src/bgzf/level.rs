use libdeflater::CompressionLvl;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Level(u8);

impl Level {
	pub const STORED: Self = Self(0);
	pub const FASTEST: Self = Self(1);
	pub const DEFAULT: Self = Self(6);
	pub const SMALLEST: Self = Self(12);

	pub fn new(level: u8) -> Option<Self> {
		(level <= Self::SMALLEST.0).then_some(Self(level))
	}

	pub(crate) fn deflate(self) -> CompressionLvl {
		CompressionLvl::new(i32::from(self.0)).unwrap_or_default()
	}
}

impl Default for Level {
	fn default() -> Self {
		Self::DEFAULT
	}
}
