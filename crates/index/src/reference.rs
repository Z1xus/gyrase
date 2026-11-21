use std::collections::HashMap;

use gyrase_io::{Chunk, VirtualOffset};

use crate::{Error, Geometry};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bin {
	pub id: u32,
	pub loff: VirtualOffset,
	pub chunks: Vec<Chunk>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Metadata {
	pub begin: VirtualOffset,
	pub end: VirtualOffset,
	pub mapped: u64,
	pub unmapped: u64,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RefIndex {
	bins: Vec<Bin>,
	slots: HashMap<u32, usize>,
	linear: Vec<VirtualOffset>,
}

impl RefIndex {
	pub fn new(bins: Vec<Bin>, linear: Vec<VirtualOffset>) -> Result<Self, Error> {
		let mut slots = HashMap::with_capacity(bins.len());
		for (position, bin) in bins.iter().enumerate() {
			if slots.insert(bin.id, position).is_some() {
				return Err(Error::DuplicateBin { bin: bin.id });
			}
		}
		Ok(Self {
			bins,
			slots,
			linear,
		})
	}

	pub fn bins(&self) -> &[Bin] {
		&self.bins
	}

	pub fn bin(&self, id: u32) -> Option<&Bin> {
		self.bins.get(*self.slots.get(&id)?)
	}

	pub fn linear(&self) -> &[VirtualOffset] {
		&self.linear
	}

	pub fn metadata(&self, geometry: Geometry) -> Option<Metadata> {
		let [positions, counts, ..] = self.bin(geometry.metadata_bin())?.chunks.as_slice() else {
			return None;
		};
		Some(Metadata {
			begin: positions.begin,
			end: positions.end,
			mapped: counts.begin.raw(),
			unmapped: counts.end.raw(),
		})
	}
}

pub(crate) fn fill_loff(geometry: Geometry, bins: &mut [Bin], linear: &[VirtualOffset]) {
	for bin in bins {
		bin.loff = geometry
			.bottom_leaf(bin.id)
			.and_then(|bottom| linear.get(usize::try_from(bottom).ok()?))
			.copied()
			.unwrap_or_default();
	}
}
