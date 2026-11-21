use std::fmt;

use gyrase_genome::Base;

use crate::{ChipSet, Chrom};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
	Same,
	Opposite,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NameCode {
	Rs(u32),
	Named(u32),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Site {
	pub(crate) chrom: Chrom,
	pub(crate) grch37: u32,
	pub(crate) grch38: u32,
	pub(crate) name: NameCode,
	pub(crate) reference: Option<Base>,
	pub(crate) alternate: Option<Base>,
	pub(crate) orientation: Orientation,
	pub(crate) chips: ChipSet,
}

impl Site {
	pub const fn chrom(&self) -> Chrom {
		self.chrom
	}

	pub const fn grch37(&self) -> u32 {
		self.grch37
	}

	pub const fn grch38(&self) -> Option<u32> {
		if self.grch38 == 0 {
			None
		} else {
			Some(self.grch38)
		}
	}

	pub const fn reference(&self) -> Option<Base> {
		self.reference
	}

	pub const fn alternate(&self) -> Option<Base> {
		self.alternate
	}

	pub const fn orientation(&self) -> Orientation {
		self.orientation
	}

	pub const fn chips(&self) -> ChipSet {
		self.chips
	}
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SiteEntry {
	pub name: String,
	pub chrom: Chrom,
	pub grch37: u32,
	pub grch38: Option<u32>,
	pub orientation: Orientation,
	pub reference: Option<Base>,
	pub alternate: Option<Base>,
	pub chips: ChipSet,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiteName<'a> {
	Rs(u32),
	Named(&'a str),
}

impl fmt::Display for SiteName<'_> {
	fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
		match self {
			Self::Rs(number) => write!(formatter, "rs{number}"),
			Self::Named(name) => formatter.write_str(name),
		}
	}
}

pub(crate) fn rs_number(name: &str) -> Option<u32> {
	let digits = name.strip_prefix("rs")?;
	if digits.starts_with('0') || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
		return None;
	}
	digits.parse().ok()
}
