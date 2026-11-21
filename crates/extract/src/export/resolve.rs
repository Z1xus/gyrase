use super::layout::Resolved;
use super::{CallBuild, Settings};
use crate::{Allele, Genotype, Ploidy, Site, SiteCall, expected_ploidy};

pub(super) fn resolve(call: &SiteCall, site: &Site, settings: &Settings) -> Resolved {
	if !settings.thresholds.accepts(call) {
		return Resolved::NoCall;
	}
	let flip = match settings.build {
		CallBuild::Grch37 => false,
		CallBuild::Grch38 => {
			if site.grch38().is_none() {
				return Resolved::NoCall;
			}
			site.orientation() == crate::Orientation::Opposite
		}
	};
	let ploidy = expected_ploidy(settings.sex, site.chrom(), site.grch37());
	apply_ploidy(call.genotype, ploidy, flip)
}

fn apply_ploidy(genotype: Genotype, ploidy: Ploidy, flip: bool) -> Resolved {
	let orient = |allele: Allele| if flip { allele.complement() } else { allele };
	match (ploidy, genotype) {
		(Ploidy::Diploid, Genotype::Pair(first, second)) => {
			Resolved::Two(orient(first), orient(second))
		}
		(Ploidy::Diploid, Genotype::Hemizygous(allele)) => {
			Resolved::Two(orient(allele), orient(allele))
		}
		(Ploidy::Haploid, Genotype::Hemizygous(allele)) => Resolved::One(orient(allele)),
		(Ploidy::Haploid, Genotype::Pair(first, second)) if first == second => {
			Resolved::One(orient(first))
		}
		(Ploidy::Absent, _) | (_, Genotype::Missing) | (Ploidy::Haploid, Genotype::Pair(..)) => {
			Resolved::NoCall
		}
	}
}
