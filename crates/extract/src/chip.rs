#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Chip {
	TwentyThreeAndMeV3,
	TwentyThreeAndMeV4,
	TwentyThreeAndMeV5,
	AncestryV1,
	AncestryV2,
	FtdnaV3,
	FtdnaV4,
	MyHeritageV1,
	MyHeritageGsa,
	LivingDnaV1,
	LivingDnaV2,
}

impl Chip {
	pub const ALL: [Self; 11] = [
		Self::TwentyThreeAndMeV3,
		Self::TwentyThreeAndMeV4,
		Self::TwentyThreeAndMeV5,
		Self::AncestryV1,
		Self::AncestryV2,
		Self::FtdnaV3,
		Self::FtdnaV4,
		Self::MyHeritageV1,
		Self::MyHeritageGsa,
		Self::LivingDnaV1,
		Self::LivingDnaV2,
	];

	pub const fn name(self) -> &'static str {
		match self {
			Self::TwentyThreeAndMeV3 => "23andMe v3",
			Self::TwentyThreeAndMeV4 => "23andMe v4",
			Self::TwentyThreeAndMeV5 => "23andMe v5",
			Self::AncestryV1 => "AncestryDNA v1",
			Self::AncestryV2 => "AncestryDNA v2",
			Self::FtdnaV3 => "FTDNA v3",
			Self::FtdnaV4 => "FTDNA v4",
			Self::MyHeritageV1 => "MyHeritage v1",
			Self::MyHeritageGsa => "MyHeritage GSA",
			Self::LivingDnaV1 => "LivingDNA v1",
			Self::LivingDnaV2 => "LivingDNA v2",
		}
	}

	const fn bit(self) -> u16 {
		1 << self as u16
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct ChipSet(u16);

impl ChipSet {
	pub const EMPTY: Self = Self(0);

	pub const fn from_bits(bits: u16) -> Option<Self> {
		if bits >> Chip::ALL.len() == 0 {
			Some(Self(bits))
		} else {
			None
		}
	}

	pub const fn bits(self) -> u16 {
		self.0
	}

	pub const fn with(self, chip: Chip) -> Self {
		Self(self.0 | chip.bit())
	}

	pub const fn contains(self, chip: Chip) -> bool {
		self.0 & chip.bit() != 0
	}

	pub const fn union(self, other: Self) -> Self {
		Self(self.0 | other.0)
	}

	pub const fn is_empty(self) -> bool {
		self.0 == 0
	}
}
