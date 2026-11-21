use gyrase_reads::Flags;

const SKIPPED: [Flags; 4] = [
	Flags::UNMAPPED,
	Flags::SECONDARY,
	Flags::QC_FAIL,
	Flags::DUPLICATE,
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Filter {
	pub min_mapping_quality: u8,
	pub min_base_quality: u8,
	pub max_depth: u32,
	pub anomalous_pairs: bool,
	pub mismatch_penalty: u8,
}

impl Default for Filter {
	fn default() -> Self {
		Self {
			min_mapping_quality: 0,
			min_base_quality: 13,
			max_depth: 1000,
			anomalous_pairs: false,
			mismatch_penalty: 10,
		}
	}
}

impl Filter {
	pub(crate) fn accepts(self, flags: Flags, mapping_quality: u8) -> bool {
		mapping_quality >= self.min_mapping_quality
			&& !SKIPPED.iter().any(|&skipped| flags.contains(skipped))
			&& (self.anomalous_pairs || !flags.is_paired() || flags.is_proper_pair())
	}
}
