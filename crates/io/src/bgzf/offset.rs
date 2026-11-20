#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct VirtualOffset(u64);

impl VirtualOffset {
	pub const fn new(block: u64, within: u16) -> Self {
		Self(block << 16 | within as u64)
	}

	pub const fn from_raw(raw: u64) -> Self {
		Self(raw)
	}

	pub const fn raw(self) -> u64 {
		self.0
	}

	pub const fn block(self) -> u64 {
		self.0 >> 16
	}

	pub const fn within(self) -> u16 {
		(self.0 & 0xffff) as u16
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct BlockPosition {
	pub ordinal: u64,
	pub within: u16,
}
