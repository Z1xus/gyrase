use gyrase_genome::Base;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Allele {
	A,
	C,
	G,
	T,
	Insertion,
	Deletion,
}

impl Allele {
	pub const fn ascii(self) -> u8 {
		match self {
			Self::A => b'A',
			Self::C => b'C',
			Self::G => b'G',
			Self::T => b'T',
			Self::Insertion => b'I',
			Self::Deletion => b'D',
		}
	}

	pub const fn complement(self) -> Self {
		match self {
			Self::A => Self::T,
			Self::C => Self::G,
			Self::G => Self::C,
			Self::T => Self::A,
			Self::Insertion => Self::Insertion,
			Self::Deletion => Self::Deletion,
		}
	}

	pub const fn from_base(base: Base) -> Option<Self> {
		match base {
			Base::A => Some(Self::A),
			Base::C => Some(Self::C),
			Base::G => Some(Self::G),
			Base::T => Some(Self::T),
			Base::N => None,
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Genotype {
	Missing,
	Hemizygous(Allele),
	Pair(Allele, Allele),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SiteCall {
	pub genotype: Genotype,
	pub depth: u32,
	pub quality: u32,
}

impl SiteCall {
	pub const MISSING: Self = Self {
		genotype: Genotype::Missing,
		depth: 0,
		quality: 0,
	};
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CallThresholds {
	pub min_depth: u32,
	pub min_quality: u32,
}

impl CallThresholds {
	pub const DEFAULT: Self = Self {
		min_depth: 10,
		min_quality: 20,
	};

	pub(crate) const fn accepts(self, call: &SiteCall) -> bool {
		call.depth >= self.min_depth && call.quality >= self.min_quality
	}
}
