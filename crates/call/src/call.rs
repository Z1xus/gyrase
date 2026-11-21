use gyrase_genome::Base;
use gyrase_pileup::{Column, Pileup};

use crate::calibration::IndelRates;
use crate::genotype::MAX_ALLELES;
use crate::indel::{self, Allele, Candidate, Locus};
use crate::likelihood::{Emissions, Evidence, Likelihoods, Priors};
use crate::realign;
use crate::{Genotype, Model, Ploidy, snp};

const SHADOW_DISTANCE: u64 = 5;
const SHADOW_RATIO: u32 = 4;

#[derive(Debug, Clone, PartialEq)]
pub struct Call {
	pub position: u64,
	pub reference: Vec<u8>,
	pub alternates: Vec<Vec<u8>>,
	pub genotype: Genotype,
	pub quality: f64,
	pub genotype_quality: u8,
	pub depth: u32,
	pub allele_depths: Vec<u32>,
	pub likelihoods: Vec<u32>,
}

impl Call {
	pub fn is_indel(&self) -> bool {
		self.reference.len() != 1 || self.alternates.iter().any(|alternate| alternate.len() != 1)
	}

	pub fn alternate_depth(&self) -> u32 {
		self.allele_depths.iter().skip(1).sum()
	}
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
	pub position: u64,
	pub reference: Vec<u8>,
	pub alternates: Vec<Vec<u8>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Caller {
	pub model: Model,
	pub ploidy: Ploidy,
	pub indel_rates: IndelRates,
	emissions: Emissions,
}

struct Scored {
	likelihoods: Likelihoods,
	genotype: Genotype,
	quality: f64,
	genotype_quality: u8,
}

impl Caller {
	pub fn new(ploidy: Ploidy) -> Self {
		Self {
			model: Model::default(),
			ploidy,
			indel_rates: IndelRates::default(),
			emissions: Emissions::default(),
		}
	}

	pub fn discover(&self, pileup: &Pileup<'_>, calls: &mut Vec<Call>) {
		let first = calls.len();
		let snp_priors = Priors::Novel {
			heterozygosity: self.model.heterozygosity,
		};
		let candidates = indel::candidates(pileup);
		for column in pileup.columns() {
			if let Some(call) = self.discover_snp(pileup, &candidates, &column, &snp_priors) {
				calls.push(call);
			}
		}
		for candidate in &candidates {
			let alternates = candidate
				.alleles
				.iter()
				.take(MAX_ALLELES - 1)
				.map(|&(allele, _)| allele)
				.collect();
			let locus = indel::locus(
				pileup,
				&self.indel_rates,
				&self.emissions,
				&candidates,
				candidate.anchor,
				alternates,
			);
			if let Some(call) = self.indel_call(pileup, locus, Purpose::Discover) {
				calls.push(call);
			}
		}
		calls[first..].sort_by_key(|call| call.position);
		drop_shadowed_indels(calls, first);
	}

	pub fn genotype(&self, pileup: &Pileup<'_>, sites: &[Site], calls: &mut Vec<Call>) {
		let candidates = indel::candidates(pileup);
		for site in sites.iter().filter(|site| {
			site.position >= pileup.interval.start && site.position < pileup.interval.end
		}) {
			calls.push(self.genotype_site(pileup, &candidates, site));
		}
	}

	fn discover_snp<'a>(
		&self,
		pileup: &Pileup<'a>,
		candidates: &[Candidate<'a>],
		column: &Column<'a>,
		priors: &Priors,
	) -> Option<Call> {
		let mut alleles = snp::alleles(column, &self.model)?;
		let (mut evidence, support) = self.snp_evidence(pileup, candidates, column, &alleles);
		let mut scored = self.score_evidence(alleles.len(), &evidence, priors);
		let used = used_alleles(scored.genotype, alleles.len());
		if used.len() < alleles.len() {
			for item in &mut evidence {
				item.select(&used);
			}
			alleles = used.iter().map(|&index| alleles[index]).collect();
			scored = self.score_evidence(alleles.len(), &evidence, priors);
		}
		if scored.genotype.is_reference() || scored.quality < self.model.min_quality {
			return None;
		}
		Some(snp_call(column, &alleles, &evidence, support, &scored))
	}

	fn snp_evidence<'a>(
		&self,
		pileup: &Pileup<'a>,
		candidates: &[Candidate<'a>],
		column: &Column<'a>,
		alleles: &[Base],
	) -> (Vec<Evidence>, Support) {
		let competitors = indel::competitors(candidates, column.position, &[]);
		if competitors.is_empty() {
			(
				snp::evidence(column, alleles, &self.model, &self.emissions),
				Support::Column,
			)
		} else {
			(
				snp::realigned(
					pileup,
					&self.model,
					&self.emissions,
					column.position,
					alleles,
					&competitors,
				),
				Support::Realigned,
			)
		}
	}

	fn score_evidence(&self, alleles: usize, evidence: &[Evidence], priors: &Priors) -> Scored {
		self.score(
			Likelihoods::compute(self.ploidy, alleles, evidence.iter().copied()),
			priors,
		)
	}

	fn indel_call<'a>(
		&self,
		pileup: &Pileup<'a>,
		mut locus: Locus<'a>,
		purpose: Purpose,
	) -> Option<Call> {
		let priors = &match purpose {
			Purpose::Discover => Priors::Novel {
				heterozygosity: self.indel_rates.heterozygosity(locus.tract),
			},
			Purpose::Genotype => Priors::Flat,
		};
		let mut scored = self.score_indel(&locus, priors);
		let used = used_alleles(scored.genotype, locus.alternates.len() + 1);
		if purpose == Purpose::Discover && used.len() < locus.alternates.len() + 1 {
			locus.retain(&used);
			scored = self.score_indel(&locus, priors);
		}
		if purpose == Purpose::Discover
			&& (scored.genotype.is_reference() || scored.quality < self.model.min_quality)
		{
			return None;
		}
		let (reference, alternates) =
			indel::alleles_text(pileup.window, locus.anchor, &locus.alternates)?;
		let allele_depths = locus.depths();
		Some(Call {
			position: locus.anchor,
			reference,
			alternates,
			genotype: scored.genotype,
			quality: scored.quality,
			genotype_quality: scored.genotype_quality,
			depth: allele_depths.iter().sum(),
			allele_depths,
			likelihoods: scored.likelihoods.phred_scaled().collect(),
		})
	}

	fn score_indel(&self, locus: &Locus<'_>, priors: &Priors) -> Scored {
		self.score_evidence(locus.alternates.len() + 1, &locus.evidence, priors)
	}

	fn score(&self, likelihoods: Likelihoods, priors: &Priors) -> Scored {
		let posterior = likelihoods.posterior(priors);
		Scored {
			likelihoods,
			genotype: self.ploidy.genotype(posterior.best),
			quality: posterior.quality,
			genotype_quality: posterior.genotype_quality,
		}
	}

	fn genotype_site<'a>(
		&self,
		pileup: &Pileup<'a>,
		candidates: &[Candidate<'a>],
		site: &'a Site,
	) -> Call {
		let missing = Call {
			position: site.position,
			reference: site.reference.clone(),
			alternates: site.alternates.clone(),
			genotype: Genotype::Missing(self.ploidy),
			quality: 0.0,
			genotype_quality: 0,
			depth: 0,
			allele_depths: Vec::new(),
			likelihoods: Vec::new(),
		};
		let is_snp = site.reference.len() == 1
			&& site.alternates.iter().all(|alternate| alternate.len() == 1);
		let known = if is_snp {
			self.genotype_snp_site(pileup, candidates, site)
		} else {
			self.genotype_indel_site(pileup, candidates, site)
		};
		known.filter(|call| call.depth > 0).unwrap_or(missing)
	}

	fn genotype_snp_site<'a>(
		&self,
		pileup: &Pileup<'a>,
		candidates: &[Candidate<'a>],
		site: &Site,
	) -> Option<Call> {
		let column = pileup.column(site.position)?;
		let reference = Base::from_ascii(*site.reference.first()?);
		reference.index()?;
		if reference != column.reference {
			return None;
		}
		let mut alleles = vec![reference];
		for base in site
			.alternates
			.iter()
			.filter_map(|alternate| alternate.first().map(|&base| Base::from_ascii(base)))
			.chain(snp::alleles(&column, &self.model).unwrap_or_default())
			.chain(Base::ACGT)
		{
			if alleles.len() < MAX_ALLELES && base.index().is_some() && !alleles.contains(&base) {
				alleles.push(base);
			}
		}
		let (evidence, support) = self.snp_evidence(pileup, candidates, &column, &alleles);
		let scored = self.score_evidence(alleles.len(), &evidence, &Priors::Flat);
		let mut call = snp_call(&column, &alleles, &evidence, support, &scored);
		call.reference.clone_from(&site.reference);
		Some(call)
	}

	fn genotype_indel_site<'a>(
		&self,
		pileup: &Pileup<'a>,
		candidates: &[Candidate<'a>],
		site: &'a Site,
	) -> Option<Call> {
		let alternates: Vec<Allele<'a>> = site
			.alternates
			.iter()
			.filter_map(|alternate| indel::from_text(&site.reference, alternate))
			.take(MAX_ALLELES - 1)
			.collect();
		if alternates.is_empty() {
			return None;
		}
		let locus = indel::locus(
			pileup,
			&self.indel_rates,
			&self.emissions,
			candidates,
			site.position,
			alternates,
		);
		self.indel_call(pileup, locus, Purpose::Genotype)
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Purpose {
	Discover,
	Genotype,
}

fn drop_shadowed_indels(calls: &mut Vec<Call>, first: usize) {
	let mut shard = calls.split_off(first);
	let indels: Vec<(u64, u32)> = shard
		.iter()
		.filter(|call| call.is_indel())
		.map(|call| (call.position, call.alternate_depth()))
		.collect();
	shard.retain(|call| {
		!call.is_indel()
			|| !indels.iter().any(|&(position, support)| {
				position != call.position
					&& position.abs_diff(call.position) <= SHADOW_DISTANCE
					&& support >= call.alternate_depth() * SHADOW_RATIO
			})
	});
	calls.append(&mut shard);
}

fn used_alleles(genotype: Genotype, count: usize) -> Vec<usize> {
	(0..count)
		.filter(|&index| {
			index == 0
				|| genotype
					.alleles()
					.any(|allele| usize::from(allele) == index)
		})
		.collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Support {
	Column,
	Realigned,
}

fn snp_call(
	column: &Column<'_>,
	alleles: &[Base],
	evidence: &[Evidence],
	support: Support,
	scored: &Scored,
) -> Call {
	let (depth, allele_depths) = match support {
		Support::Column => (column.observations.len(), snp::depths(column, alleles)),
		Support::Realigned => (evidence.len(), realign::depths(evidence, alleles.len())),
	};
	Call {
		position: column.position,
		reference: vec![column.reference.ascii()],
		alternates: alleles
			.iter()
			.skip(1)
			.map(|base| vec![base.ascii()])
			.collect(),
		genotype: scored.genotype,
		quality: scored.quality,
		genotype_quality: scored.genotype_quality,
		depth: u32::try_from(depth).unwrap_or(u32::MAX),
		allele_depths,
		likelihoods: scored.likelihoods.phred_scaled().collect(),
	}
}
