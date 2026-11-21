use crate::{Chrom, Sex};

const X_PAR_GRCH37: [(u32, u32); 2] = [(60_001, 2_699_520), (154_931_044, 155_260_560)];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ploidy {
	Diploid,
	Haploid,
	Absent,
}

pub fn in_x_par(grch37: u32) -> bool {
	X_PAR_GRCH37
		.iter()
		.any(|&(start, end)| (start..=end).contains(&grch37))
}

pub fn expected_ploidy(sex: Sex, chrom: Chrom, grch37: u32) -> Ploidy {
	match (chrom, sex) {
		(Chrom::MT, _) | (Chrom::Y, Sex::Male) => Ploidy::Haploid,
		(Chrom::X, Sex::Male) if !in_x_par(grch37) => Ploidy::Haploid,
		(Chrom::Y, Sex::Female) => Ploidy::Absent,
		_ => Ploidy::Diploid,
	}
}
