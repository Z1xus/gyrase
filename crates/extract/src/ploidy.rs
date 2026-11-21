use gyrase_genome::{Build, Chromosome};

use crate::Sex;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ploidy {
	Diploid,
	Haploid,
	Absent,
}

pub fn in_x_par(grch37: u32) -> bool {
	Build::Grch37
		.pseudoautosomal()
		.iter()
		.any(|range| range.contains(&(u64::from(grch37) - 1)))
}

pub fn expected_ploidy(sex: Sex, chromosome: Chromosome, grch37: u32) -> Ploidy {
	match (chromosome, sex) {
		(Chromosome::Mitochondrion, _) | (Chromosome::Y, Sex::Male) => Ploidy::Haploid,
		(Chromosome::X, Sex::Male) if !in_x_par(grch37) => Ploidy::Haploid,
		(Chromosome::Y, Sex::Female) => Ploidy::Absent,
		_ => Ploidy::Diploid,
	}
}
