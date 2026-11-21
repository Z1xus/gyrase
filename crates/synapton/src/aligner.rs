#[cfg(feature = "gpu")]
use std::sync::atomic::{AtomicUsize, Ordering};
#[cfg(feature = "gpu")]
use std::sync::{Arc, Condvar, Mutex, PoisonError};

use gyrase_io::Workers;
use rayon::prelude::*;

use crate::Error;
use crate::alignment::AlignedBatch;
use crate::config::AlignConfig;
#[cfg(feature = "gpu")]
use crate::front::FrontBatch;
#[cfg(feature = "gpu")]
use crate::gpu::Gpu;
use crate::index::Index;
use crate::num::u32_of;
use crate::pairing::InsertModel;
#[cfg(feature = "gpu")]
use crate::worker::Front;
use crate::worker::{Layout, Worker};

#[cfg(feature = "gpu")]
const FIRST_PIECE_CHUNKS: usize = 8;
#[cfg(feature = "gpu")]
const MAX_PIECE_CHUNKS: usize = 32;

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
	#[cfg(feature = "gpu")]
	gpu: Option<Gpu>,
}

impl<'a> Aligner<'a> {
	pub fn new(index: &'a Index, config: AlignConfig, workers: &'a Workers) -> Self {
		Self {
			index,
			config,
			workers,
			#[cfg(feature = "gpu")]
			gpu: config
				.gpu
				.then(|| {
					Gpu::new(
						index,
						&config,
						2 * MAX_PIECE_CHUNKS * config.chunk_reads.max(1),
					)
				})
				.flatten(),
		}
	}

	pub fn uses_gpu(&self) -> bool {
		#[cfg(feature = "gpu")]
		return self.gpu.is_some();
		#[cfg(not(feature = "gpu"))]
		false
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
		#[cfg(feature = "gpu")]
		if let Some(gpu) = &self.gpu {
			return self.run_gpu(gpu, sequences, ends, layout, init, each);
		}
		self.workers.install(|| {
			sequences
				.par_chunks(chunk)
				.enumerate()
				.map_init(
					|| (Worker::new(self.index, self.config), init()),
					|(worker, state), (number, part)| {
						let mut out = AlignedBatch::default();
						let first_read = u32_of(number * self.config.chunk_reads.max(1));
						worker.align(part, layout, first_read, None, &mut out);
						each(state, out)
					},
				)
				.collect()
		})
	}

	#[cfg(feature = "gpu")]
	fn run_gpu<S, T: Send>(
		&self,
		gpu: &Gpu,
		sequences: &[&[u8]],
		ends: usize,
		layout: Layout<'_>,
		init: impl Fn() -> S + Sync,
		each: impl Fn(&mut S, AlignedBatch) -> T + Sync,
	) -> Vec<T> {
		let chunk = ends * self.config.chunk_reads.max(1);
		let parts: Vec<&[&[u8]]> = sequences.chunks(chunk).collect();
		let mut bounds = vec![0];
		let mut size = FIRST_PIECE_CHUNKS;
		while let Some(&last) = bounds.last().filter(|&&last| last < parts.len()) {
			bounds.push((last + size).min(parts.len()));
			size = (2 * size).min(MAX_PIECE_CHUNKS);
		}
		let pieces = Pieces::new(bounds.len() - 1);
		let results: Vec<Mutex<Option<T>>> = parts.iter().map(|_| Mutex::new(None)).collect();
		let next = AtomicUsize::new(0);
		std::thread::scope(|scope| {
			scope.spawn(|| {
				for (number, edge) in bounds.windows(2).enumerate() {
					let piece = sequences
						.get(edge[0] * chunk..(edge[1] * chunk).min(sequences.len()))
						.unwrap_or_default();
					pieces.publish(number, gpu.front(piece).map(Arc::new));
				}
			});
			self.workers.install(|| {
				rayon::scope(|workers| {
					for _ in 0..rayon::current_num_threads() {
						workers.spawn(|_| {
							let mut worker = Worker::new(self.index, self.config);
							let mut state = init();
							loop {
								let number = next.fetch_add(1, Ordering::Relaxed);
								let (Some(part), Some(slot)) =
									(parts.get(number), results.get(number))
								else {
									break;
								};
								let piece = bounds.partition_point(|&edge| edge <= number) - 1;
								let front = pieces.wait(piece);
								let mut out = AlignedBatch::default();
								worker.align(
									part,
									layout,
									u32_of(number * self.config.chunk_reads.max(1)),
									front.as_deref().map(|batch| Front {
										batch,
										first: (number - bounds[piece]) * chunk,
									}),
									&mut out,
								);
								*slot.lock().unwrap_or_else(PoisonError::into_inner) =
									Some(each(&mut state, out));
							}
						});
					}
				});
			});
		});
		results
			.into_iter()
			.filter_map(|slot| slot.into_inner().unwrap_or_else(PoisonError::into_inner))
			.collect()
	}
}

#[cfg(feature = "gpu")]
enum Piece {
	Pending,
	Ready(Option<Arc<FrontBatch>>),
}

#[cfg(feature = "gpu")]
struct Pieces {
	ready: Mutex<Vec<Piece>>,
	signal: Condvar,
}

#[cfg(feature = "gpu")]
impl Pieces {
	fn new(count: usize) -> Self {
		Self {
			ready: Mutex::new((0..count).map(|_| Piece::Pending).collect()),
			signal: Condvar::new(),
		}
	}

	fn publish(&self, piece: usize, batch: Option<Arc<FrontBatch>>) {
		let mut ready = self.ready.lock().unwrap_or_else(PoisonError::into_inner);
		if let Some(slot) = ready.get_mut(piece) {
			*slot = Piece::Ready(batch);
		}
		self.signal.notify_all();
	}

	fn wait(&self, piece: usize) -> Option<Arc<FrontBatch>> {
		let mut ready = self.ready.lock().unwrap_or_else(PoisonError::into_inner);
		loop {
			if let Some(Piece::Ready(batch)) = ready.get(piece) {
				return batch.clone();
			}
			ready = self
				.signal
				.wait(ready)
				.unwrap_or_else(PoisonError::into_inner);
		}
	}
}

fn merge(parts: Vec<AlignedBatch>) -> AlignedBatch {
	let mut merged = AlignedBatch::default();
	for part in parts {
		merged.absorb(part);
	}
	merged
}
