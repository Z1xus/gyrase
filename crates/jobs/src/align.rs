use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{self, SyncSender};
use std::thread;
use std::time::Instant;

use gyrase_genome::Reference;
use gyrase_io::{Decoded, Workers, decode};
use gyrase_reads::{
	FastqBatch, FastqReader, Fields, Header, InterleavedFastq, MarkSettings, MarkStats,
	PairedBatch, PairedFastq, ReadRef, Sorter, interleaved, interleaved_pairs,
};
use synapton::{
	AlignConfig, AlignedBatch, Aligner, BamEmitter, EmitError, Index, InsertSizes, ReadPair,
	SeedConfig, Source,
};

use crate::bam::write_sorted;
use crate::counted::Counted;
use crate::{Cancel, Error, Estimate, Input, Lane, Progress, Stage};

const BATCH_READS: usize = 1 << 16;
const QUEUE: usize = 2;

enum Chunk {
	Pairs(PairedBatch),
	Interleaved(FastqBatch),
	Single(FastqBatch),
}

enum Views<'a> {
	Pairs(Vec<ReadPair<'a>>),
	Reads(Vec<synapton::Read<'a>>),
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Alignment<'a> {
	pub index: &'a Index,
	pub compute: &'a Workers,
	pub io: &'a Workers,
	pub cancel: &'a Cancel,
}

pub(crate) fn index(reference: &Reference, path: &Path, workers: &Workers) -> Result<Index, Error> {
	if let Some(parent) = path.parent() {
		std::fs::create_dir_all(parent).map_err(|source| Error::Write {
			path: parent.to_owned(),
			source,
		})?;
	}
	Ok(Index::load_or_build(
		path,
		reference,
		SeedConfig::short_read(),
		workers,
	)?)
}

impl Alignment<'_> {
	pub(crate) fn run(
		self,
		lanes: &[Lane],
		sample: &str,
		output: &Path,
		budget: usize,
		report: &mut impl FnMut(Progress),
	) -> Result<Option<MarkStats>, Error> {
		let mut header =
			Header::from_contigs(self.index.contigs()).map_err(gyrase_reads::Error::from)?;
		header.push_read_group(Fields::new().with(*b"ID", sample).with(*b"SM", sample));
		let total = lanes.iter().flat_map(inputs).map(|input| input.size).sum();
		let read = Arc::new(AtomicU64::new(0));
		let mut estimate = Estimate::new(Stage::Align, total, Instant::now());
		report(estimate.progress());
		let stats = write_sorted(
			&header,
			output,
			budget,
			Some(MarkSettings::default()),
			self.compute,
			self.cancel,
			|sorter| {
				self.align(lanes, sample, sorter, &read, |done| {
					report(estimate.advance(done, Instant::now()));
				})?;
				report(Estimate::new(Stage::Sort, 1, Instant::now()).progress());
				Ok(())
			},
		)?;
		report(Estimate::new(Stage::Sort, 1, Instant::now()).advance(1, Instant::now()));
		Ok(stats)
	}

	fn align(
		self,
		lanes: &[Lane],
		sample: &str,
		sorter: &mut Sorter<'_>,
		read: &Arc<AtomicU64>,
		mut advance: impl FnMut(u64),
	) -> Result<(), Error> {
		let aligner = Aligner::new(self.index, AlignConfig::default(), self.compute);
		thread::scope(|scope| {
			let (chunk_sender, chunks) = mpsc::sync_channel(QUEUE);
			let (batch_sender, batches) = mpsc::sync_channel::<(Chunk, AlignedBatch)>(QUEUE);
			let reader = scope.spawn(move || read_lanes(lanes, self.io, read, &chunk_sender));
			let writer = scope.spawn(move || {
				let mut emitter = BamEmitter::new(self.index.contigs(), Some(sample.as_bytes()));
				for (chunk, batch) in batches {
					emitter
						.emit(&batch, chunk.views().source(), |record| sorter.push(record))
						.map_err(|error| match error {
							EmitError::Convert(error) => Error::Align(error),
							EmitError::Sink(error) => Error::Reads(error),
						})?;
				}
				Ok::<_, Error>(())
			});
			let mut sizes = None;
			for chunk in chunks {
				self.cancel.check()?;
				let batch = align_chunk(&aligner, &mut sizes, &chunk.views())?;
				if batch_sender.send((chunk, batch)).is_err() {
					break;
				}
				advance(read.load(Ordering::Relaxed));
			}
			drop(batch_sender);
			writer.join().map_err(|_| Error::Worker)??;
			reader.join().map_err(|_| Error::Worker)?
		})
	}
}

impl Chunk {
	fn views(&self) -> Views<'_> {
		match self {
			Self::Pairs(batch) => Views::Pairs(
				batch
					.pairs()
					.map(|pair| ReadPair {
						first: unaligned(pair.first),
						second: unaligned(pair.second),
					})
					.collect(),
			),
			Self::Interleaved(batch) => Views::Pairs(
				interleaved_pairs(batch)
					.map(|pair| ReadPair {
						first: unaligned(pair.first),
						second: unaligned(pair.second),
					})
					.collect(),
			),
			Self::Single(batch) => Views::Reads(batch.iter().map(unaligned).collect()),
		}
	}
}

impl Views<'_> {
	fn source(&self) -> Source<'_> {
		match self {
			Self::Pairs(pairs) => Source::Pairs(pairs),
			Self::Reads(reads) => Source::Reads(reads),
		}
	}
}

fn align_chunk(
	aligner: &Aligner<'_>,
	sizes: &mut Option<InsertSizes>,
	views: &Views<'_>,
) -> Result<AlignedBatch, Error> {
	Ok(match views {
		Views::Reads(reads) => aligner.align_reads(reads)?,
		Views::Pairs(pairs) => {
			let sizes = match sizes {
				Some(sizes) => *sizes,
				None => *sizes.insert(aligner.estimate_insert_sizes(pairs)?),
			};
			aligner.align_pairs_with(&sizes, pairs)?
		}
	})
}

fn read_lanes(
	lanes: &[Lane],
	io: &Workers,
	read: &Arc<AtomicU64>,
	sender: &SyncSender<Chunk>,
) -> Result<(), Error> {
	for lane in lanes {
		match lane {
			Lane::Paired(first, second) => {
				let mut reader = PairedFastq::new(open(first, read, io)?, open(second, read, io)?);
				pump(sender, || {
					let mut batch = PairedBatch::new();
					let pairs = reader.read_batch(&mut batch, BATCH_READS / 2)?;
					Ok((pairs > 0).then_some(Chunk::Pairs(batch)))
				})?;
			}
			Lane::Single(input) if paired_within(input, io)? => {
				let mut reader = InterleavedFastq::new(open(input, read, io)?);
				pump(sender, || {
					let mut batch = FastqBatch::new();
					let pairs = reader.read_batch(&mut batch, BATCH_READS / 2)?;
					Ok((pairs > 0).then_some(Chunk::Interleaved(batch)))
				})?;
			}
			Lane::Single(input) => {
				let mut reader = FastqReader::new(open(input, read, io)?);
				pump(sender, || {
					let mut batch = FastqBatch::new();
					let reads = reader.read_batch(&mut batch, BATCH_READS)?;
					Ok((reads > 0).then_some(Chunk::Single(batch)))
				})?;
			}
		}
	}
	Ok(())
}

fn pump(
	sender: &SyncSender<Chunk>,
	mut next: impl FnMut() -> Result<Option<Chunk>, Error>,
) -> Result<(), Error> {
	while let Some(chunk) = next()? {
		if sender.send(chunk).is_err() {
			break;
		}
	}
	Ok(())
}

fn paired_within(input: &Input, io: &Workers) -> Result<bool, Error> {
	let mut batch = FastqBatch::new();
	FastqReader::new(open(input, &Arc::default(), io)?).read_batch(&mut batch, 2)?;
	Ok(interleaved(&batch))
}

fn open(input: &Input, read: &Arc<AtomicU64>, io: &Workers) -> Result<Decoded, Error> {
	let inner = File::open(&input.path).map_err(|source| Error::Open {
		path: input.path.clone(),
		source,
	})?;
	decode(Box::new(Counted::new(inner, read)), io).map_err(|source| Error::Read {
		path: input.path.clone(),
		source,
	})
}

fn inputs(lane: &Lane) -> Vec<&Input> {
	match lane {
		Lane::Single(input) => vec![input],
		Lane::Paired(first, second) => vec![first, second],
	}
}

fn unaligned(read: ReadRef<'_>) -> synapton::Read<'_> {
	synapton::Read {
		name: read.name,
		sequence: read.bases,
		qualities: read.qualities,
	}
}
