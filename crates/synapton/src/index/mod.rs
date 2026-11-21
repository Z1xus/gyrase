mod build;
mod entry;
mod file;

use std::ops::Range;
use std::path::Path;
use std::sync::atomic::AtomicBool;

use gyrase_genome::{Contigs, Reference};
use gyrase_io::Workers;

use crate::Error;
use crate::config::SeedConfig;
use crate::reference::PackedReference;
use crate::seeding::Seed;

pub(crate) use entry::{offset_of as entry_offset, position_of as entry_position};

#[derive(Debug)]
pub struct Index {
	config: SeedConfig,
	bucket_bits: u32,
	filter_cutoff: u32,
	buckets: Vec<u32>,
	entries: Vec<u64>,
	reference: PackedReference,
}

#[derive(Debug, Clone)]
pub(crate) struct Probe {
	range: Range<usize>,
	prefix: u64,
}

impl Index {
	pub fn build(
		reference: &Reference,
		config: SeedConfig,
		workers: &Workers,
		stop: &AtomicBool,
	) -> Result<Self, Error> {
		let packed = PackedReference::from_fasta(reference, workers)?;
		build::build(packed, config, workers, stop)
	}

	pub fn load_or_build(
		path: &Path,
		reference: &Reference,
		config: SeedConfig,
		workers: &Workers,
		stop: &AtomicBool,
	) -> Result<Self, Error> {
		let packed = PackedReference::from_fasta(reference, workers)?;
		let packed = match file::reuse(path, packed, config, workers) {
			file::Reuse::Loaded(index) => return Ok(*index),
			file::Reuse::Stale(packed) => packed,
		};
		let index = build::build(packed, config, workers, stop)?;
		index.save(path, workers)?;
		Ok(index)
	}

	pub fn load(path: &Path, workers: &Workers) -> Result<Self, Error> {
		file::load(path, workers)
	}

	pub fn save(&self, path: &Path, workers: &Workers) -> Result<(), Error> {
		file::save(self, path, workers)
	}

	pub fn contigs(&self) -> &Contigs {
		self.reference.contigs()
	}

	pub(crate) fn reference(&self) -> &PackedReference {
		&self.reference
	}

	pub fn config(&self) -> SeedConfig {
		self.config
	}

	pub fn seed_count(&self) -> usize {
		self.entries.len()
	}

	pub fn filter_cutoff(&self) -> usize {
		self.filter_cutoff as usize
	}

	pub fn heap_bytes(&self) -> usize {
		size_of_val(self.entries.as_slice())
			+ size_of_val(self.buckets.as_slice())
			+ size_of_val(self.reference.words())
	}

	pub(crate) fn probe(&self, seed: &Seed) -> Probe {
		let (bucket, tag) = entry::split_key(seed.key, self.bucket_bits);
		let start = self.buckets.get(bucket).copied().unwrap_or(0) as usize;
		let end = self.buckets.get(bucket + 1).copied().unwrap_or(0) as usize;
		Probe {
			range: start..end,
			prefix: entry::prefix(tag, seed.orientation),
		}
	}

	pub(crate) fn first_strobe_run(&self, probe: &Probe) -> impl Iterator<Item = u64> + '_ {
		let main = entry::main_tag_of(probe.prefix);
		let orientation = entry::first_orientation_of(probe.prefix);
		self.entries
			.get(probe.range.clone())
			.unwrap_or_default()
			.iter()
			.copied()
			.filter(move |&found| {
				entry::main_tag_of(entry::prefix_of(found)) == main
					&& entry::first_orientation_of(entry::prefix_of(found)) == orientation
			})
	}

	pub(crate) fn run(&self, probe: &Probe) -> &[u64] {
		let bucket = self.entries.get(probe.range.clone()).unwrap_or_default();
		let first = bucket
			.iter()
			.position(|&found| entry::prefix_of(found) >= probe.prefix)
			.unwrap_or(bucket.len());
		let tail = bucket.get(first..).unwrap_or_default();
		let length = tail
			.iter()
			.take_while(|&&found| entry::prefix_of(found) == probe.prefix)
			.count();
		tail.get(..length).unwrap_or_default()
	}
}
