use gyrase_io::Workers;
use rayon::prelude::*;

use crate::Error;
use crate::alignment::AlignedBatch;
use crate::config::AlignConfig;
use crate::extend::Extender;
use crate::index::Index;
use crate::kernel::CpuExtender;
use crate::num::u32_of;
use crate::pairing::InsertModel;
use crate::worker::{Layout, Worker};

pub const MAX_READ_LENGTH: usize = 1024;

#[derive(Debug, Clone, Copy)]
pub struct Read<'a> {
	pub name: &'a [u8],
	pub sequence: &'a [u8],
	pub qualities: &'a [u8],
}

impl Read<'_> {
	fn validate(&self, index: usize) -> Result<(), Error> {
		if self.sequence.len() > MAX_READ_LENGTH {
			return Err(Error::ReadTooLong {
				index,
				length: self.sequence.len(),
				limit: MAX_READ_LENGTH,
			});
		}
		if !self.qualities.is_empty() && self.qualities.len() != self.sequence.len() {
			return Err(Error::QualityLength {
				index,
				bases: self.sequence.len(),
				qualities: self.qualities.len(),
			});
		}
		Ok(())
	}
}

#[derive(Debug, Clone, Copy)]
pub struct ReadPair<'a> {
	pub first: Read<'a>,
	pub second: Read<'a>,
}

#[derive(Debug)]
pub struct Aligner<'a, E: Extender = CpuExtender> {
	index: &'a Index,
	config: AlignConfig,
	workers: &'a Workers,
	extender: E,
}

impl<'a> Aligner<'a> {
	pub fn new(index: &'a Index, config: AlignConfig, workers: &'a Workers) -> Self {
		Self::with_extender(index, config, workers, CpuExtender)
	}
}

impl<'a, E: Extender> Aligner<'a, E> {
	pub fn with_extender(
		index: &'a Index,
		config: AlignConfig,
		workers: &'a Workers,
		extender: E,
	) -> Self {
		Self {
			index,
			config,
			workers,
			extender,
		}
	}

	pub fn align_reads(&self, reads: &[Read<'_>]) -> Result<AlignedBatch, Error> {
		for (index, read) in reads.iter().enumerate() {
			read.validate(index)?;
		}
		let sequences: Vec<&[u8]> = reads.iter().map(|read| read.sequence).collect();
		Ok(self.run(&sequences, 1, Layout::Single))
	}

	pub fn align_pairs(&self, pairs: &[ReadPair<'_>]) -> Result<AlignedBatch, Error> {
		for (index, pair) in pairs.iter().enumerate() {
			pair.first.validate(index)?;
			pair.second.validate(index)?;
		}
		let sequences: Vec<&[u8]> = pairs
			.iter()
			.flat_map(|pair| [pair.first.sequence, pair.second.sequence])
			.collect();
		let model = self.estimate_insert_model(&sequences);
		Ok(self.run(&sequences, 2, Layout::Paired(&model)))
	}

	fn estimate_insert_model(&self, sequences: &[&[u8]]) -> InsertModel {
		let sample = sequences
			.get(..sequences.len().min(2 * self.config.calibration_pairs))
			.unwrap_or_default();
		let chunk = 2 * self.config.chunk_reads.max(1);
		let mut distances: Vec<i64> = self.workers.install(|| {
			sample
				.par_chunks(chunk)
				.map_init(
					|| Worker::new(self.index, self.config, &self.extender),
					|worker, part| {
						let mut found = Vec::new();
						worker.distances(part, &mut found);
						found
					},
				)
				.flatten()
				.collect()
		});
		InsertModel::estimate(&mut distances)
	}

	fn run(&self, sequences: &[&[u8]], ends: usize, layout: Layout<'_>) -> AlignedBatch {
		let chunk = ends * self.config.chunk_reads.max(1);
		let parts: Vec<AlignedBatch> = self.workers.install(|| {
			sequences
				.par_chunks(chunk)
				.enumerate()
				.map_init(
					|| Worker::new(self.index, self.config, &self.extender),
					|worker, (number, part)| {
						let mut out = AlignedBatch::default();
						let first_read = u32_of(number * self.config.chunk_reads.max(1));
						worker.align(part, layout, first_read, &mut out);
						out
					},
				)
				.collect()
		});
		let mut merged = AlignedBatch::default();
		for part in parts {
			merged.absorb(part);
		}
		merged
	}
}
