use std::num::NonZeroU32;

use gyrase_pileup::{Change, Indel, Pileup, Window};

use crate::calibration::{IndelRates, Tract};
use crate::likelihood::Emissions;
use crate::likelihood::Evidence;
use crate::realign::{self, Edit, realign};

const MIN_SUPPORT: u32 = 2;
const COMPETITOR_DISTANCE: u64 = 48;
const MAX_COMPETITORS: usize = 4;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Allele<'a> {
	Insertion(&'a [u8]),
	Deletion(u32),
}

impl<'a> Allele<'a> {
	pub fn deletion(self) -> u64 {
		match self {
			Self::Deletion(len) => u64::from(len),
			Self::Insertion(_) => 0,
		}
	}

	pub fn inserted(self) -> &'a [u8] {
		match self {
			Self::Insertion(bases) => bases,
			Self::Deletion(_) => &[],
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Variant<'a> {
	pub anchor: u64,
	pub allele: Allele<'a>,
}

impl<'a> Variant<'a> {
	pub fn edit(self) -> Edit<'a> {
		Edit {
			start: self.anchor + 1,
			deleted: self.allele.deletion(),
			inserted: self.allele.inserted(),
		}
	}
}

#[derive(Debug, Clone)]
pub(crate) struct Candidate<'a> {
	pub anchor: u64,
	pub alleles: Vec<(Allele<'a>, u32)>,
}

#[derive(Debug, Clone)]
pub(crate) struct Locus<'a> {
	pub anchor: u64,
	pub tract: Tract,
	pub alternates: Vec<Allele<'a>>,
	pub evidence: Vec<Evidence>,
}

impl Locus<'_> {
	pub fn retain(&mut self, used: &[usize]) {
		for item in &mut self.evidence {
			let values = item.likelihoods;
			let reference = (0..=self.alternates.len())
				.filter(|index| !used.contains(index))
				.map(|index| values[index])
				.fold(values[0], f64::max);
			item.likelihoods = std::array::from_fn(|slot| match slot {
				0 => reference,
				_ => used
					.get(slot)
					.map_or(f64::NEG_INFINITY, |&index| values[index]),
			});
		}
		self.alternates = used
			.iter()
			.skip(1)
			.map(|&index| self.alternates[index - 1])
			.collect();
	}

	pub fn depths(&self) -> Vec<u32> {
		realign::depths(&self.evidence, self.alternates.len() + 1)
	}
}

pub(crate) fn candidates<'a>(pileup: &Pileup<'a>) -> Vec<Candidate<'a>> {
	pileup
		.indels()
		.chunk_by(|left, right| left.anchor == right.anchor)
		.filter_map(|group| candidate(pileup, group))
		.collect()
}

fn candidate<'a>(pileup: &Pileup<'a>, group: &'a [Indel]) -> Option<Candidate<'a>> {
	let anchor = group.first()?.anchor;
	let mut alleles: Vec<(Allele<'a>, u32)> = Vec::new();
	for allele in distinct(pileup, group) {
		match alleles.iter_mut().find(|(known, _)| *known == allele) {
			Some((_, count)) => *count += 1,
			None => alleles.push((allele, 1)),
		}
	}
	alleles.retain(|&(_, count)| count >= MIN_SUPPORT);
	alleles.sort_by(|left, right| {
		right
			.1
			.cmp(&left.1)
			.then_with(|| order(left.0).cmp(&order(right.0)))
	});
	(!alleles.is_empty()).then_some(Candidate { anchor, alleles })
}

pub(crate) fn competitors<'a>(
	candidates: &[Candidate<'a>],
	anchor: u64,
	alternates: &[Allele<'a>],
) -> Vec<Edit<'a>> {
	let begin =
		candidates.partition_point(|candidate| candidate.anchor + COMPETITOR_DISTANCE < anchor);
	let end =
		candidates.partition_point(|candidate| candidate.anchor <= anchor + COMPETITOR_DISTANCE);
	let mut nearby: Vec<(Variant<'a>, u32)> = candidates[begin..end.max(begin)]
		.iter()
		.flat_map(|candidate| {
			candidate.alleles.iter().map(|&(allele, support)| {
				(
					Variant {
						anchor: candidate.anchor,
						allele,
					},
					support,
				)
			})
		})
		.filter(|(variant, _)| variant.anchor != anchor || !alternates.contains(&variant.allele))
		.collect();
	nearby.sort_by(|left, right| {
		right
			.1
			.cmp(&left.1)
			.then_with(|| left.0.anchor.cmp(&right.0.anchor))
	});
	nearby
		.into_iter()
		.take(MAX_COMPETITORS)
		.map(|(variant, _)| variant.edit())
		.collect()
}

pub(crate) fn locus<'a>(
	pileup: &Pileup<'a>,
	rates: &IndelRates,
	emissions: &Emissions,
	candidates: &[Candidate<'a>],
	anchor: u64,
	alternates: Vec<Allele<'a>>,
) -> Locus<'a> {
	let tract = Tract::at(pileup.window, anchor);
	let quality = rates.error_quality(tract);
	let edits: Vec<Edit<'a>> = alternates
		.iter()
		.map(|&allele| Variant { anchor, allele }.edit())
		.collect();
	let evidence = realign(
		pileup,
		emissions,
		anchor..anchor + 1 + tract.span as u64,
		&edits,
		&competitors(candidates, anchor, &alternates),
		|read| Some(quality.min(read.mapping_quality)),
	);
	Locus {
		anchor,
		tract,
		alternates,
		evidence,
	}
}

pub(crate) fn alleles_text(
	window: Window<'_>,
	anchor: u64,
	alternates: &[Allele<'_>],
) -> Option<(Vec<u8>, Vec<Vec<u8>>)> {
	let longest = alternates
		.iter()
		.map(|allele| match allele {
			Allele::Deletion(len) => *len as usize,
			Allele::Insertion(_) => 0,
		})
		.max()
		.unwrap_or(0);
	let start = usize::try_from(anchor.checked_sub(window.start)?).ok()?;
	let reference = window.bases.get(start..=start + longest)?.to_vec();
	let alternates = alternates
		.iter()
		.map(|allele| {
			let mut text = vec![reference[0]];
			match allele {
				Allele::Deletion(len) => text.extend_from_slice(&reference[1 + *len as usize..]),
				Allele::Insertion(bases) => {
					text.extend_from_slice(bases);
					text.extend_from_slice(&reference[1..]);
				}
			}
			text
		})
		.collect();
	Some((reference, alternates))
}

pub(crate) fn from_text<'a>(reference: &[u8], alternate: &'a [u8]) -> Option<Allele<'a>> {
	match (reference, alternate) {
		([first, ..], [only]) if reference.len() > 1 && first == only => {
			u32::try_from(reference.len() - 1)
				.ok()
				.map(Allele::Deletion)
		}
		([only], [first, inserted @ ..]) if only == first && !inserted.is_empty() => {
			Some(Allele::Insertion(inserted))
		}
		_ => None,
	}
}

fn distinct<'a>(pileup: &Pileup<'a>, group: &'a [Indel]) -> impl Iterator<Item = Allele<'a>> {
	let mut seen: Vec<NonZeroU32> = Vec::new();
	let pileup = *pileup;
	group.iter().filter_map(move |indel| {
		if let Some(fragment) = indel.fragment {
			if seen.contains(&fragment) {
				return None;
			}
			seen.push(fragment);
		}
		Some(allele(&pileup, indel))
	})
}

fn allele<'a>(pileup: &Pileup<'a>, indel: &Indel) -> Allele<'a> {
	match indel.change {
		Change::Insertion { start, len } => Allele::Insertion(pileup.inserted(start, len)),
		Change::Deletion { len } => Allele::Deletion(len),
	}
}

fn order(allele: Allele<'_>) -> (u8, u32, &[u8]) {
	match allele {
		Allele::Deletion(len) => (0, len, &[]),
		Allele::Insertion(bases) => (1, u32::try_from(bases.len()).unwrap_or(u32::MAX), bases),
	}
}
