use gyrase_io::Workers;
use rayon::prelude::*;

use crate::Error;
use crate::alignment::AlignedBatch;
use crate::config::AlignConfig;
use crate::index::Index;
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

#[derive(Debug, Clone, Copy)]
pub struct InsertSizes(InsertModel);

fn pair_sequences<'a>(pairs: &[ReadPair<'a>]) -> Result<Vec<&'a [u8]>, Error> {
	let mut sequences = Vec::with_capacity(2 * pairs.len());
	for (index, pair) in pairs.iter().enumerate() {
		pair.first.validate(index)?;
		pair.second.validate(index)?;
		sequences.push(pair.first.sequence);
		sequences.push(pair.second.sequence);
	}
	Ok(sequences)
}

#[derive(Debug)]
pub struct Aligner<'a> {
	index: &'a Index,
	config: AlignConfig,
	workers: &'a Workers,
}

impl<'a> Aligner<'a> {
	pub fn new(index: &'a Index, config: AlignConfig, workers: &'a Workers) -> Self {
		Self {
			index,
			config,
			workers,
		}
	}

	pub fn align_reads(&self, reads: &[Read<'_>]) -> Result<AlignedBatch, Error> {
		Ok(merge(self.align_reads_each(
			reads,
			|| (),
			|(), part| part,
		)?))
	}

	pub fn align_reads_each<S, T: Send>(
		&self,
		reads: &[Read<'_>],
		init: impl Fn() -> S + Sync,
		each: impl Fn(&mut S, AlignedBatch) -> T + Sync,
	) -> Result<Vec<T>, Error> {
		for (index, read) in reads.iter().enumerate() {
			read.validate(index)?;
		}
		let sequences: Vec<&[u8]> = reads.iter().map(|read| read.sequence).collect();
		Ok(self.run(&sequences, 1, Layout::Single, init, each))
	}

	pub fn align_pairs(&self, pairs: &[ReadPair<'_>]) -> Result<AlignedBatch, Error> {
		let sizes = self.estimate_insert_sizes(pairs)?;
		self.align_pairs_with(&sizes, pairs)
	}

	pub fn estimate_insert_sizes(&self, pairs: &[ReadPair<'_>]) -> Result<InsertSizes, Error> {
		let sequences = pair_sequences(pairs)?;
		Ok(InsertSizes(self.estimate_insert_model(&sequences)))
	}

	pub fn align_pairs_with(
		&self,
		sizes: &InsertSizes,
		pairs: &[ReadPair<'_>],
	) -> Result<AlignedBatch, Error> {
		Ok(merge(self.align_pairs_each(
			sizes,
			pairs,
			|| (),
			|(), part| part,
		)?))
	}

	pub fn align_pairs_each<S, T: Send>(
		&self,
		sizes: &InsertSizes,
		pairs: &[ReadPair<'_>],
		init: impl Fn() -> S + Sync,
		each: impl Fn(&mut S, AlignedBatch) -> T + Sync,
	) -> Result<Vec<T>, Error> {
		let sequences = pair_sequences(pairs)?;
		Ok(self.run(&sequences, 2, Layout::Paired(&sizes.0), init, each))
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
					|| Worker::new(self.index, self.config),
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

	fn run<S, T: Send>(
		&self,
		sequences: &[&[u8]],
		ends: usize,
		layout: Layout<'_>,
		init: impl Fn() -> S + Sync,
		each: impl Fn(&mut S, AlignedBatch) -> T + Sync,
	) -> Vec<T> {
		let chunk = ends * self.config.chunk_reads.max(1);
		self.workers.install(|| {
			sequences
				.par_chunks(chunk)
				.enumerate()
				.map_init(
					|| (Worker::new(self.index, self.config), init()),
					|(worker, state), (number, part)| {
						let mut out = AlignedBatch::default();
						let first_read = u32_of(number * self.config.chunk_reads.max(1));
						worker.align(part, layout, first_read, &mut out);
						each(state, out)
					},
				)
				.collect()
		})
	}
}

fn merge(parts: Vec<AlignedBatch>) -> AlignedBatch {
	let mut merged = AlignedBatch::default();
	for part in parts {
		merged.absorb(part);
	}
	merged
}
