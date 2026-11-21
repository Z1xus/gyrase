use gyrase_pileup::{Pileup, Window};

const MAX_PERIOD: usize = 6;
const MAX_SPAN: usize = 40;
const PSEUDO_DEPTH: f64 = 20_000.0;
const VARIANT_FRACTION: f64 = 0.1;
const CLEAN_QUALITY: f64 = 55.0;
const MIN_QUALITY: f64 = 10.0;
const NOVEL_HETEROZYGOSITY: f64 = 1.25e-4;
const MAX_HETEROZYGOSITY: f64 = 0.05;
const PSEUDO_LOCI: f64 = 2_000.0;
const MIN_LOCUS_DEPTH: u64 = 10;
const MIN_VARIANT_EVENTS: u64 = 3;
const VARIANT_LOCUS_FRACTION: f64 = 0.25;
const HETEROZYGOSITY_FOLD_PER_QUALITY: f64 = 20.0;
const VARIANT_SHARE_OF_HETEROZYGOSITY: f64 = 1.5;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Tract {
	pub period: usize,
	pub span: usize,
}

type Table<T> = [[T; MAX_SPAN + 1]; MAX_PERIOD];

#[derive(Debug, Clone, PartialEq)]
pub struct IndelRates {
	qualities: Table<u8>,
	heterozygosities: Table<f64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Calibration {
	depth: Table<u64>,
	events: Table<u64>,
	loci: Table<u64>,
	variants: Table<u64>,
}

impl Tract {
	pub fn at(window: Window<'_>, anchor: u64) -> Self {
		let after = anchor
			.checked_sub(window.start)
			.and_then(|offset| usize::try_from(offset + 1).ok())
			.and_then(|offset| window.bases.get(offset..))
			.unwrap_or_default();
		(1..=MAX_PERIOD)
			.map(|period| Self {
				period,
				span: period + periodic_run(after, period),
			})
			.filter(|tract| tract.period == 1 || tract.span >= 2 * tract.period)
			.max_by_key(|tract| (tract.span, std::cmp::Reverse(tract.period)))
			.unwrap_or(Self { period: 1, span: 1 })
	}

	fn cell(self) -> (usize, usize) {
		(
			self.period.clamp(1, MAX_PERIOD) - 1,
			self.span.min(MAX_SPAN),
		)
	}
}

impl Default for IndelRates {
	fn default() -> Self {
		Calibration::default().rates()
	}
}

impl IndelRates {
	pub fn error_quality(&self, tract: Tract) -> u8 {
		let (period, span) = tract.cell();
		self.qualities[period][span]
	}

	pub fn heterozygosity(&self, tract: Tract) -> f64 {
		let (period, span) = tract.cell();
		self.heterozygosities[period][span]
	}
}

impl Default for Calibration {
	fn default() -> Self {
		Self {
			depth: [[0; MAX_SPAN + 1]; MAX_PERIOD],
			events: [[0; MAX_SPAN + 1]; MAX_PERIOD],
			loci: [[0; MAX_SPAN + 1]; MAX_PERIOD],
			variants: [[0; MAX_SPAN + 1]; MAX_PERIOD],
		}
	}
}

impl Calibration {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn observe(&mut self, pileup: &Pileup<'_>) {
		for column in pileup.columns() {
			let tract = Tract::at(pileup.window, column.position);
			let required = u8::try_from(tract.span + 1).unwrap_or(u8::MAX);
			let depth = column
				.observations
				.iter()
				.filter(|observation| observation.flank >= required)
				.count() as u64;
			let events = pileup
				.indels_at(column.position)
				.iter()
				.filter(|indel| indel.flank >= required)
				.count() as u64;
			if depth == 0 {
				continue;
			}
			let (period, span) = tract.cell();
			if depth >= MIN_LOCUS_DEPTH {
				self.loci[period][span] += 1;
				if events >= MIN_VARIANT_EVENTS && ratio(events, depth) >= VARIANT_LOCUS_FRACTION {
					self.variants[period][span] += 1;
				}
			}
			if events >= 2 && ratio(events, depth) >= VARIANT_FRACTION {
				continue;
			}
			self.depth[period][span] += depth;
			self.events[period][span] += events.min(depth);
		}
	}

	pub fn merge(&mut self, other: &Self) {
		for period in 0..MAX_PERIOD {
			for span in 0..=MAX_SPAN {
				self.depth[period][span] += other.depth[period][span];
				self.events[period][span] += other.events[period][span];
				self.loci[period][span] += other.loci[period][span];
				self.variants[period][span] += other.variants[period][span];
			}
		}
	}

	pub fn rates(&self) -> IndelRates {
		let mut qualities = [[0; MAX_SPAN + 1]; MAX_PERIOD];
		let mut heterozygosities = [[0.0; MAX_SPAN + 1]; MAX_PERIOD];
		for period in 0..MAX_PERIOD {
			let mut ceiling = CLEAN_QUALITY;
			let mut floor = 0.0_f64;
			for span in 0..=MAX_SPAN {
				let prior = prior_quality(period + 1, span);
				let error = 10f64.powf(-prior / 10.0);
				let depth = real(self.depth[period][span]);
				let events = real(self.events[period][span]);
				let rate = (events + error * PSEUDO_DEPTH) / (depth + PSEUDO_DEPTH);
				ceiling = ceiling.min(-10.0 * rate.log10());
				qualities[period][span] = to_quality(ceiling);
				let expected = NOVEL_HETEROZYGOSITY
					* 10f64.powf((CLEAN_QUALITY - prior) / HETEROZYGOSITY_FOLD_PER_QUALITY);
				let loci = real(self.loci[period][span]);
				let variants = real(self.variants[period][span]) / VARIANT_SHARE_OF_HETEROZYGOSITY;
				floor = floor.max((variants + expected * PSEUDO_LOCI) / (loci + PSEUDO_LOCI));
				heterozygosities[period][span] =
					floor.clamp(NOVEL_HETEROZYGOSITY, MAX_HETEROZYGOSITY);
			}
		}
		IndelRates {
			qualities,
			heterozygosities,
		}
	}
}

fn prior_quality(period: usize, span: usize) -> f64 {
	let span = real(span.max(period) as u64);
	let quality = if period == 1 {
		if span <= 10.0 {
			54.0 - 3.2 * (span - 1.0)
		} else {
			25.2 - 0.6 * (span - 10.0)
		}
	} else {
		52.0 - 2.5 * (span / real(period as u64) - 2.0).max(0.0)
	};
	quality.clamp(MIN_QUALITY, CLEAN_QUALITY)
}

fn periodic_run(bases: &[u8], period: usize) -> usize {
	bases
		.iter()
		.zip(bases.iter().skip(period))
		.take_while(|(left, right)| left == right)
		.count()
}

fn ratio(numerator: u64, denominator: u64) -> f64 {
	real(numerator) / real(denominator)
}

#[allow(clippy::cast_precision_loss)]
fn real(value: u64) -> f64 {
	value as f64
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn to_quality(value: f64) -> u8 {
	value.round().clamp(MIN_QUALITY, CLEAN_QUALITY) as u8
}
