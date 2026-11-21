use crate::filter::variant_type;
use crate::record::Record;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct Kinds(u8);

const SNP: u8 = 1;
const MNP: u8 = 2;
const INDEL: u8 = 4;
const OTHER: u8 = 8;
const BREAKEND: u8 = 16;
const OVERLAP: u8 = 32;

impl Kinds {
	pub(crate) fn of(record: &Record) -> Self {
		let reference = record.allele(0).unwrap_or_default();
		let bits = record
			.alleles()
			.skip(1)
			.map(|alt| match variant_type(reference, alt) {
				"snp" => SNP,
				"mnp" => MNP,
				"indel" => INDEL,
				"bnd" => BREAKEND,
				"overlap" => OVERLAP,
				"other" => OTHER,
				_ => 0,
			})
			.fold(0, |all, bit| all | bit);
		Self(bits)
	}

	pub(crate) const fn has_breakend(self) -> bool {
		self.0 & BREAKEND != 0
	}

	pub(crate) const fn has_overlap(self) -> bool {
		self.0 & OVERLAP != 0
	}

	pub(crate) const fn has_substitution(self) -> bool {
		self.0 & (SNP | MNP) != 0
	}

	pub(crate) const fn has_snp(self) -> bool {
		self.0 & SNP != 0
	}

	pub(crate) const fn has_indel(self) -> bool {
		self.0 & INDEL != 0
	}

	pub(crate) const fn union(self, other: Self) -> Self {
		Self(self.0 | other.0)
	}
}
