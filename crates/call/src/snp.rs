use std::cmp::Reverse;

use gyrase_genome::Base;
use gyrase_pileup::{Column, Observation, Pileup};

use crate::Model;
use crate::genotype::MAX_ALLELES;
use crate::likelihood::{Emissions, Evidence};
use crate::realign::{Edit, realign};

const MIN_ALTERNATE_QUALITY_SUM: u32 = 20;

pub(crate) fn alleles(column: &Column<'_>, model: &Model) -> Option<Vec<Base>> {
	column.reference.index()?;
	let mut sums = [0u32; 4];
	for observation in column.observations {
		if let Some(index) = observation.base.index() {
			sums[index] += u32::from(model.base_quality(observation));
		}
	}
	let mut alternates: Vec<(u32, Base)> = Base::ACGT
		.into_iter()
		.zip(sums)
		.filter(|&(base, sum)| base != column.reference && sum >= MIN_ALTERNATE_QUALITY_SUM)
		.map(|(base, sum)| (sum, base))
		.collect();
	if alternates.is_empty() {
		return None;
	}
	alternates.sort_by(|left, right| {
		right
			.0
			.cmp(&left.0)
			.then(left.1.ascii().cmp(&right.1.ascii()))
	});
	let mut alleles = vec![column.reference];
	alleles.extend(
		alternates
			.into_iter()
			.take(MAX_ALLELES - 1)
			.map(|(_, base)| base),
	);
	Some(alleles)
}

pub(crate) fn evidence(
	column: &Column<'_>,
	alleles: &[Base],
	model: &Model,
	emissions: &Emissions,
) -> Vec<Evidence> {
	let mut ranked: Vec<(&Observation, u8)> = column
		.observations
		.iter()
		.map(|observation| (observation, model.base_quality(observation)))
		.collect();
	ranked.sort_by_key(|&(observation, quality)| {
		(
			observation.base.ascii(),
			observation.strand,
			Reverse(quality),
		)
	});
	let mut evidence = Vec::with_capacity(ranked.len());
	let mut repeat = 0;
	let mut previous = None;
	for (observation, quality) in ranked {
		let group = Some((observation.base, observation.strand));
		repeat = if group == previous { repeat + 1 } else { 0 };
		previous = group;
		evidence.push(emissions.base(
			|allele| alleles.get(allele) == Some(&observation.base),
			model.repeated_quality(quality, repeat),
		));
	}
	evidence
}

pub(crate) fn depths(column: &Column<'_>, alleles: &[Base]) -> Vec<u32> {
	alleles
		.iter()
		.map(|&allele| {
			column
				.observations
				.iter()
				.filter(|observation| observation.base == allele)
				.count()
		})
		.map(|count| u32::try_from(count).unwrap_or(u32::MAX))
		.collect()
}

const SUBSTITUTIONS: [u8; 4] = *b"ACGT";

pub(crate) fn realigned(
	pileup: &Pileup<'_>,
	model: &Model,
	emissions: &Emissions,
	position: u64,
	alleles: &[Base],
	competitors: &[Edit<'_>],
) -> Vec<Evidence> {
	let alternates: Vec<Edit<'_>> = alleles
		.iter()
		.skip(1)
		.filter_map(|base| base.index())
		.map(|index| Edit {
			start: position,
			deleted: 1,
			inserted: &SUBSTITUTIONS[index..=index],
		})
		.collect();
	realign(
		pileup,
		emissions,
		position..position + 1,
		&alternates,
		competitors,
		|read| {
			let (_, bases, qualities) = read.slice(position, position + 1)?;
			let (&base, &quality) = bases.first().zip(qualities.first())?;
			let flank = (position - read.start).min(read.end - 1 - position);
			(base != b'N' && quality >= pileup.filter.min_base_quality)
				.then(|| model.quality(quality, read.mapping_quality, flank))
		},
	)
}
