use crate::contig::canonical;
use crate::{Contig, Contigs};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Build {
	Grch37,
	Grch38,
	Chm13,
}

const MITOCHONDRION: &str = "MT";

const NAMES: [&str; 25] = [
	"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
	"18", "19", "20", "21", "22", "X", "Y", "MT",
];

const GRCH37: [u64; 25] = [
	249_250_621,
	243_199_373,
	198_022_430,
	191_154_276,
	180_915_260,
	171_115_067,
	159_138_663,
	146_364_022,
	141_213_431,
	135_534_747,
	135_006_516,
	133_851_895,
	115_169_878,
	107_349_540,
	102_531_392,
	90_354_753,
	81_195_210,
	78_077_248,
	59_128_983,
	63_025_520,
	48_129_895,
	51_304_566,
	155_270_560,
	59_373_566,
	16_569,
];

const GRCH38: [u64; 25] = [
	248_956_422,
	242_193_529,
	198_295_559,
	190_214_555,
	181_538_259,
	170_805_979,
	159_345_973,
	145_138_636,
	138_394_717,
	133_797_422,
	135_086_622,
	133_275_309,
	114_364_328,
	107_043_718,
	101_991_189,
	90_338_345,
	83_257_441,
	80_373_285,
	58_617_616,
	64_444_167,
	46_709_983,
	50_818_468,
	156_040_895,
	57_227_415,
	16_569,
];

const CHM13: [u64; 25] = [
	248_387_328,
	242_696_752,
	201_105_948,
	193_574_945,
	182_045_439,
	172_126_628,
	160_567_428,
	146_259_331,
	150_617_247,
	134_758_134,
	135_127_769,
	133_324_548,
	113_566_686,
	101_161_492,
	99_753_195,
	96_330_374,
	84_276_897,
	80_542_538,
	61_707_364,
	66_210_255,
	45_090_682,
	51_324_926,
	154_259_566,
	62_460_029,
	16_569,
];

impl Build {
	pub const ALL: [Self; 3] = [Self::Grch37, Self::Grch38, Self::Chm13];

	pub const fn name(self) -> &'static str {
		match self {
			Self::Grch37 => "GRCh37",
			Self::Grch38 => "GRCh38",
			Self::Chm13 => "CHM13",
		}
	}

	pub fn detect(contigs: &Contigs) -> Option<Self> {
		let votes = Self::ALL.map(|build| {
			(
				build,
				contigs
					.iter()
					.filter(|(_, contig)| build.identifies(contig))
					.count(),
			)
		});
		let best = votes
			.iter()
			.map(|(_, count)| *count)
			.max()
			.filter(|&count| count > 0)?;
		let mut leaders = votes.iter().filter(|(_, count)| *count == best);
		let (build, _) = leaders.next()?;
		leaders.next().is_none().then_some(*build)
	}

	fn identifies(self, contig: &Contig) -> bool {
		canonical(&contig.name) != MITOCHONDRION && self.length(&contig.name) == Some(contig.length)
	}

	pub fn length(self, chromosome: &str) -> Option<u64> {
		let name = canonical(chromosome);
		let index = NAMES.iter().position(|candidate| *candidate == name)?;
		Some(self.lengths()[index])
	}

	const fn lengths(self) -> &'static [u64; 25] {
		match self {
			Self::Grch37 => &GRCH37,
			Self::Grch38 => &GRCH38,
			Self::Chm13 => &CHM13,
		}
	}
}
