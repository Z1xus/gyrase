use std::ops::{BitOr, BitOrAssign};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Flags(u16);

impl Flags {
	pub const PAIRED: Self = Self(0x1);
	pub const PROPER_PAIR: Self = Self(0x2);
	pub const UNMAPPED: Self = Self(0x4);
	pub const MATE_UNMAPPED: Self = Self(0x8);
	pub const REVERSE: Self = Self(0x10);
	pub const MATE_REVERSE: Self = Self(0x20);
	pub const FIRST_IN_PAIR: Self = Self(0x40);
	pub const SECOND_IN_PAIR: Self = Self(0x80);
	pub const SECONDARY: Self = Self(0x100);
	pub const QC_FAIL: Self = Self(0x200);
	pub const DUPLICATE: Self = Self(0x400);
	pub const SUPPLEMENTARY: Self = Self(0x800);

	pub const fn new(bits: u16) -> Self {
		Self(bits)
	}

	pub const fn bits(self) -> u16 {
		self.0
	}

	pub const fn contains(self, other: Self) -> bool {
		self.0 & other.0 == other.0
	}

	pub const fn without(self, other: Self) -> Self {
		Self(self.0 & !other.0)
	}

	pub const fn is_paired(self) -> bool {
		self.contains(Self::PAIRED)
	}

	pub const fn is_proper_pair(self) -> bool {
		self.contains(Self::PROPER_PAIR)
	}

	pub const fn is_unmapped(self) -> bool {
		self.contains(Self::UNMAPPED)
	}

	pub const fn is_mapped(self) -> bool {
		!self.is_unmapped()
	}

	pub const fn is_mate_unmapped(self) -> bool {
		self.contains(Self::MATE_UNMAPPED)
	}

	pub const fn is_reverse(self) -> bool {
		self.contains(Self::REVERSE)
	}

	pub const fn is_mate_reverse(self) -> bool {
		self.contains(Self::MATE_REVERSE)
	}

	pub const fn is_first_in_pair(self) -> bool {
		self.contains(Self::FIRST_IN_PAIR)
	}

	pub const fn is_second_in_pair(self) -> bool {
		self.contains(Self::SECOND_IN_PAIR)
	}

	pub const fn is_secondary(self) -> bool {
		self.contains(Self::SECONDARY)
	}

	pub const fn is_qc_fail(self) -> bool {
		self.contains(Self::QC_FAIL)
	}

	pub const fn is_duplicate(self) -> bool {
		self.contains(Self::DUPLICATE)
	}

	pub const fn is_supplementary(self) -> bool {
		self.contains(Self::SUPPLEMENTARY)
	}

	pub const fn is_primary(self) -> bool {
		!self.is_secondary() && !self.is_supplementary()
	}
}

impl BitOr for Flags {
	type Output = Self;

	fn bitor(self, other: Self) -> Self {
		Self(self.0 | other.0)
	}
}

impl BitOrAssign for Flags {
	fn bitor_assign(&mut self, other: Self) {
		self.0 |= other.0;
	}
}
