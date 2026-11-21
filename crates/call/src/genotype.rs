#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Ploidy {
	Haploid,
	Diploid,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Genotype {
	Missing(Ploidy),
	Haploid(u8),
	Diploid(u8, u8),
}

pub(crate) const MAX_ALLELES: usize = 3;
pub(crate) const MAX_GENOTYPES: usize = MAX_ALLELES * (MAX_ALLELES + 1) / 2;

impl Ploidy {
	pub(crate) const fn genotype_count(self, alleles: usize) -> usize {
		match self {
			Self::Haploid => alleles,
			Self::Diploid => alleles * (alleles + 1) / 2,
		}
	}

	pub(crate) fn genotype(self, index: usize) -> Genotype {
		let allele = |index: usize| u8::try_from(index).unwrap_or(u8::MAX);
		match self.pair(index) {
			(first, None) => Genotype::Haploid(allele(first)),
			(first, Some(second)) => Genotype::Diploid(allele(first), allele(second)),
		}
	}

	pub(crate) fn pair(self, index: usize) -> (usize, Option<usize>) {
		match self {
			Self::Haploid => (index, None),
			Self::Diploid => {
				let mut second = 0;
				while (second + 1) * (second + 2) / 2 <= index {
					second += 1;
				}
				(index - second * (second + 1) / 2, Some(second))
			}
		}
	}
}

impl Genotype {
	pub fn is_reference(self) -> bool {
		matches!(self, Self::Haploid(0) | Self::Diploid(0, 0))
	}

	pub fn alleles(self) -> impl Iterator<Item = u8> {
		let pair = match self {
			Self::Missing(_) => [None, None],
			Self::Haploid(allele) => [Some(allele), None],
			Self::Diploid(first, second) => [Some(first), Some(second)],
		};
		pair.into_iter().flatten()
	}
}
