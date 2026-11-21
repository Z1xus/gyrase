use gyrase_pileup::Observation;

use crate::likelihood::whole;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Model {
	pub heterozygosity: f64,
	pub max_base_quality: u8,
	pub read_end_bases: u8,
	pub read_end_quality: u8,
	pub min_quality: f64,
	pub repeated_error_decay: f64,
	pub repeated_error_floor: f64,
}

impl Default for Model {
	fn default() -> Self {
		Self {
			heterozygosity: 1e-3,
			max_base_quality: 40,
			read_end_bases: 12,
			read_end_quality: 25,
			min_quality: 10.0,
			repeated_error_decay: 0.83,
			repeated_error_floor: 0.03,
		}
	}
}

impl Model {
	pub(crate) fn base_quality(&self, observation: &Observation) -> u8 {
		self.quality(
			observation.quality,
			observation.mapping_quality,
			u64::from(observation.flank),
		)
	}

	pub(crate) fn quality(&self, quality: u8, mapping_quality: u8, flank: u64) -> u8 {
		let capped = quality.min(mapping_quality).min(self.max_base_quality);
		if flank < u64::from(self.read_end_bases) {
			capped.min(self.read_end_quality)
		} else {
			capped
		}
	}

	pub(crate) fn repeated_quality(&self, quality: u8, repeat: i32) -> u8 {
		let weight = self.repeated_error_decay.powi(repeat) * (1.0 - self.repeated_error_floor)
			+ self.repeated_error_floor;
		u8::try_from(whole(f64::from(quality) * weight)).unwrap_or(u8::MAX)
	}
}
