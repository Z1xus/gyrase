use std::ops::Range;

use crate::Build;
use crate::contig::canonical;

const LABELS: [&str; 25] = [
	"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
	"18", "19", "20", "21", "22", "X", "Y", "MT",
];
const AUTOSOMES: u8 = 22;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Chromosome {
	Autosome(u8),
	X,
	Y,
	Mitochondrion,
}

impl Chromosome {
	pub const COUNT: usize = LABELS.len();

	pub fn parse(name: &str) -> Option<Self> {
		let name = canonical(name);
		LABELS
			.iter()
			.position(|label| *label == name)
			.and_then(Self::from_index)
	}

	pub fn from_index(index: usize) -> Option<Self> {
		match u8::try_from(index).ok()? {
			number if number < AUTOSOMES => Some(Self::Autosome(number + 1)),
			AUTOSOMES => Some(Self::X),
			23 => Some(Self::Y),
			24 => Some(Self::Mitochondrion),
			_ => None,
		}
	}

	pub fn index(self) -> usize {
		usize::from(match self {
			Self::Autosome(number) => number - 1,
			Self::X => AUTOSOMES,
			Self::Y => AUTOSOMES + 1,
			Self::Mitochondrion => AUTOSOMES + 2,
		})
	}

	pub fn label(self) -> &'static str {
		LABELS.get(self.index()).copied().unwrap_or_default()
	}
}

impl Build {
	pub const fn pseudoautosomal(self) -> [Range<u64>; 2] {
		match self {
			Self::Grch37 => [60_000..2_699_520, 154_931_043..155_260_560],
			Self::Grch38 => [10_000..2_781_479, 155_701_382..156_030_895],
			Self::Chm13 => [0..2_394_410, 153_925_834..154_259_566],
		}
	}
}
