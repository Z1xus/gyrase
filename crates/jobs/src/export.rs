use gyrase_call::{Call, Genotype as Called};
use gyrase_extract::{Allele, CallBuild, Genotype, Panel, SiteCall};
use gyrase_genome::{Base, Chromosome, ContigId};

use crate::{Alignments, Target};

#[derive(Debug, Clone, Default)]
pub struct Targets {
	pub targets: Vec<Target>,
	sites: Vec<usize>,
}

impl Targets {
	pub fn new(panel: &Panel, alignments: &Alignments, build: CallBuild) -> Self {
		let mut contigs = [None; Chromosome::COUNT];
		for (index, reference) in alignments.header().references().iter().enumerate() {
			if let (Some(chromosome), Ok(index)) =
				(Chromosome::parse(reference.name()), u32::try_from(index))
			{
				contigs[chromosome.index()].get_or_insert(ContigId::new(index));
			}
		}
		let mut located: Vec<(Target, usize)> = panel
			.sites()
			.iter()
			.enumerate()
			.filter_map(|(index, site)| {
				site.reference()?;
				let position = match build {
					CallBuild::Grch37 => Some(site.grch37()),
					CallBuild::Grch38 => site.grch38(),
				}?;
				let target = Target {
					contig: contigs[site.chromosome().index()]?,
					position: u64::from(position.checked_sub(1)?),
					alternate: site.alternate(),
				};
				Some((target, index))
			})
			.collect();
		located.sort_by_key(|(target, _)| (target.contig.get(), target.position));
		let (targets, sites) = located.into_iter().unzip();
		Self { targets, sites }
	}

	pub fn site_calls(&self, panel: &Panel, genotyped: &[Option<Call>]) -> Vec<SiteCall> {
		let mut calls = vec![SiteCall::MISSING; panel.len()];
		for (&site, call) in self.sites.iter().zip(genotyped) {
			if let (Some(slot), Some(call)) = (calls.get_mut(site), call) {
				*slot = site_call(call);
			}
		}
		calls
	}
}

fn site_call(call: &Call) -> SiteCall {
	let allele = |index: u8| {
		let text = match index {
			0 => Some(&call.reference),
			_ => call.alternates.get(usize::from(index) - 1),
		}?;
		match text.as_slice() {
			&[base] => Allele::from_base(Base::from_ascii(base)),
			_ => None,
		}
	};
	let genotype = match call.genotype {
		Called::Missing(_) => None,
		Called::Haploid(only) => allele(only).map(Genotype::Hemizygous),
		Called::Diploid(first, second) => allele(first)
			.zip(allele(second))
			.map(|(first, second)| Genotype::Pair(first, second)),
	};
	genotype.map_or(SiteCall::MISSING, |genotype| SiteCall {
		genotype,
		depth: call.depth,
		quality: u32::from(call.genotype_quality),
	})
}
