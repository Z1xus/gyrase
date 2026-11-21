use std::num::NonZeroU32;

use gyrase_genome::Base;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strand {
	Forward,
	Reverse,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Observation {
	pub base: Base,
	pub quality: u8,
	pub mapping_quality: u8,
	pub strand: Strand,
	pub flank: u8,
	pub fragment: Option<NonZeroU32>,
}
