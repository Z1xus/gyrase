use crate::genotype::{MAX_ALLELES, MAX_GENOTYPES, Ploidy};

const LN_TEN: f64 = std::f64::consts::LN_10;
const MAX_GENOTYPE_QUALITY: f64 = 99.0;
const RANDOM_BASE_ERROR: f64 = 0.75;
const EMISSION_QUALITIES: usize = 64;

#[derive(Debug, Clone, Copy)]
pub(crate) struct Evidence {
	pub likelihoods: [f64; MAX_ALLELES],
}

impl Evidence {
	pub fn select(&mut self, alleles: &[usize]) {
		let values = self.likelihoods;
		self.likelihoods = std::array::from_fn(|slot| {
			alleles
				.get(slot)
				.map_or(f64::NEG_INFINITY, |&index| values[index])
		});
	}
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Likelihoods {
	pub ploidy: Ploidy,
	pub alleles: usize,
	pub values: [f64; MAX_GENOTYPES],
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Posterior {
	pub best: usize,
	pub quality: f64,
	pub genotype_quality: u8,
}

pub(crate) enum Priors {
	Novel { heterozygosity: f64 },
	Flat,
}

impl Likelihoods {
	pub fn compute(
		ploidy: Ploidy,
		alleles: usize,
		evidence: impl Iterator<Item = Evidence>,
	) -> Self {
		let alleles = alleles.min(MAX_ALLELES);
		let mut values = [0.0; MAX_GENOTYPES];
		for item in evidence {
			let read = item.likelihoods;
			for (index, value) in values
				.iter_mut()
				.enumerate()
				.take(ploidy.genotype_count(alleles))
			{
				*value += match ploidy.pair(index) {
					(first, None) => read[first],
					(first, Some(second)) => log_mean(read[first], read[second]),
				};
			}
		}
		Self {
			ploidy,
			alleles,
			values,
		}
	}

	pub fn count(&self) -> usize {
		self.ploidy.genotype_count(self.alleles)
	}

	pub fn phred_scaled(&self) -> impl Iterator<Item = u32> + '_ {
		let values = &self.values[..self.count()];
		let best = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
		values.iter().map(move |value| {
			whole((-10.0 * (value - best) / LN_TEN).clamp(0.0, f64::from(u16::MAX)))
		})
	}

	pub fn posterior(&self, priors: &Priors) -> Posterior {
		let count = self.count();
		let mut posterior = [f64::NEG_INFINITY; MAX_GENOTYPES];
		for (index, slot) in posterior.iter_mut().enumerate().take(count) {
			*slot = self.values[index] + priors.ln(self.ploidy, self.alleles, index);
		}
		let peak = posterior[..count]
			.iter()
			.copied()
			.fold(f64::NEG_INFINITY, f64::max);
		let total: f64 = posterior[..count]
			.iter()
			.map(|value| (value - peak).exp())
			.sum();
		let probability = |index: usize| (posterior[index] - peak).exp() / total;
		let best = (0..count)
			.max_by(|&left, &right| posterior[left].total_cmp(&posterior[right]))
			.unwrap_or(0);
		let reference = probability(0);
		let wrong = (1.0 - probability(best)).max(0.0);
		Posterior {
			best,
			quality: phred(reference),
			genotype_quality: u8::try_from(whole(phred(wrong).min(MAX_GENOTYPE_QUALITY)))
				.unwrap_or(u8::MAX),
		}
	}
}

impl Priors {
	fn ln(&self, ploidy: Ploidy, alleles: usize, index: usize) -> f64 {
		match self {
			Self::Flat => 0.0,
			Self::Novel { heterozygosity } => {
				let alternates = alternate_count(alleles);
				let probability = match ploidy.pair(index) {
					(0, None | Some(0)) => 1.0 - 1.5 * heterozygosity,
					(_, None) | (0, Some(_)) => heterozygosity / alternates,
					(first, Some(second)) if first == second => heterozygosity / 2.0 / alternates,
					_ => heterozygosity * heterozygosity,
				};
				probability.ln()
			}
		}
	}
}

pub(crate) fn alternate_count(alleles: usize) -> f64 {
	f64::from(u8::try_from(alleles.max(2) - 1).unwrap_or(u8::MAX))
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub(crate) fn whole(value: f64) -> u32 {
	value.round().max(0.0) as u32
}

fn log_mean(first: f64, second: f64) -> f64 {
	if first.total_cmp(&second).is_eq() {
		return first;
	}
	let peak = first.max(second);
	peak + f64::midpoint((first - peak).exp(), (second - peak).exp()).ln()
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Emissions {
	matched: [f64; EMISSION_QUALITIES],
	mismatched: [f64; EMISSION_QUALITIES],
}

impl Default for Emissions {
	fn default() -> Self {
		let error = |quality: usize| phred_error(u8::try_from(quality).unwrap_or(u8::MAX));
		Self {
			matched: std::array::from_fn(|quality| (1.0 - error(quality)).ln()),
			mismatched: std::array::from_fn(|quality| (error(quality) / 3.0).ln()),
		}
	}
}

impl Emissions {
	pub fn matched(&self, quality: u8) -> f64 {
		self.matched[usize::from(quality).min(EMISSION_QUALITIES - 1)]
	}

	pub fn mismatched(&self, quality: u8) -> f64 {
		self.mismatched[usize::from(quality).min(EMISSION_QUALITIES - 1)]
	}

	pub fn base(&self, base_matches: impl Fn(usize) -> bool, quality: u8) -> Evidence {
		let (matched, mismatched) = (self.matched(quality), self.mismatched(quality));
		Evidence {
			likelihoods: std::array::from_fn(|allele| {
				if base_matches(allele) {
					matched
				} else {
					mismatched
				}
			}),
		}
	}
}

pub(crate) fn phred_to_ln(quality: u8) -> f64 {
	-f64::from(quality) * LN_TEN / 10.0
}

fn phred_error(quality: u8) -> f64 {
	10f64
		.powf(-f64::from(quality) / 10.0)
		.min(RANDOM_BASE_ERROR)
}

fn phred(probability: f64) -> f64 {
	-10.0 * probability.max(1e-300).log10()
}
