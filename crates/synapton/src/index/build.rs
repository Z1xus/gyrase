use std::sync::atomic::{AtomicBool, AtomicU32, AtomicU64, Ordering};

use gyrase_io::Workers;
use rayon::prelude::*;

use super::Index;
use super::entry::{self, MAX_BUCKET_BITS};
use crate::Error;
use crate::config::{MAX_STROBE_OFFSET, SeedConfig};
use crate::num::{i64_of, u32_from, u32_of, usize_of};
use crate::reference::{ContigSpan, PackedReference};
use crate::seeding::{Seed, Seeder};

const BLOCK_BASES: u64 = 1 << 21;
const MIN_BUCKET_BITS: u32 = 10;
const TARGET_BUCKET_LOAD: u64 = 6;
const SORT_SEGMENT_ENTRIES: u64 = 1 << 20;
const HISTOGRAM_LEN: usize = 256;
const FILTER_DENOMINATOR: u64 = 5000;
const FILTER_MIN: u32 = 30;
const FILTER_MAX: u32 = 100;

struct Block {
	span: ContigSpan,
	start: u64,
	length: u64,
}

struct Extractor {
	seeder: Seeder,
	codes: Vec<u8>,
	seeds: Vec<Seed>,
}

impl Extractor {
	fn new(config: SeedConfig) -> Self {
		Self {
			seeder: Seeder::new(config),
			codes: Vec::new(),
			seeds: Vec::new(),
		}
	}

	fn seeds(&mut self, reference: &PackedReference, block: &Block, overlap: u64) -> &[Seed] {
		let wanted = (block.length + overlap).min(block.span.length - block.start);
		reference.window(
			&block.span,
			i64_of(block.span.start + block.start),
			usize_of(wanted),
			&mut self.codes,
		);
		self.seeds.clear();
		self.seeder
			.extract(&self.codes, u32_from(block.length), &mut self.seeds);
		&self.seeds
	}
}

pub(super) fn build(
	reference: PackedReference,
	config: SeedConfig,
	workers: &Workers,
	stop: &AtomicBool,
) -> Result<Index, Error> {
	let stopped = || {
		if stop.load(Ordering::Relaxed) {
			Err(Error::Cancelled)
		} else {
			Ok(())
		}
	};
	let blocks = blocks(&reference);
	let overlap = (MAX_STROBE_OFFSET + config.k + 1) as u64;
	let estimate = estimate_seeds(&reference, &config);
	if estimate == 0 {
		return Err(Error::EmptyReference);
	}
	let bucket_bits = (estimate / TARGET_BUCKET_LOAD)
		.max(1)
		.ilog2()
		.clamp(MIN_BUCKET_BITS, MAX_BUCKET_BITS);
	let bucket_count = 1usize << bucket_bits;

	let counts: Vec<AtomicU32> = (0..bucket_count).map(|_| AtomicU32::new(0)).collect();
	workers.install(|| {
		blocks.par_iter().for_each_init(
			|| Extractor::new(config),
			|extractor, block| {
				if stop.load(Ordering::Relaxed) {
					return;
				}
				for seed in extractor.seeds(&reference, block, overlap) {
					let (bucket, _) = entry::split_key(seed.key, bucket_bits);
					if let Some(count) = counts.get(bucket) {
						count.fetch_add(1, Ordering::Relaxed);
					}
				}
			},
		);
	});
	stopped()?;

	let mut buckets = Vec::with_capacity(bucket_count + 1);
	let mut total = 0u64;
	for count in counts {
		buckets.push(u32::try_from(total).map_err(|_| Error::CorruptIndex("too many seeds"))?);
		total += u64::from(count.into_inner());
	}
	buckets.push(u32::try_from(total).map_err(|_| Error::CorruptIndex("too many seeds"))?);

	let cursors: Vec<AtomicU32> = buckets
		.iter()
		.take(bucket_count)
		.map(|&start| AtomicU32::new(start))
		.collect();
	let slots: Vec<AtomicU64> = (0..total).map(|_| AtomicU64::new(0)).collect();
	workers.install(|| {
		blocks.par_iter().for_each_init(
			|| Extractor::new(config),
			|extractor, block| {
				if stop.load(Ordering::Relaxed) {
					return;
				}
				let flat_start = block.span.start + block.start;
				for seed in extractor.seeds(&reference, block, overlap) {
					let (bucket, tag) = entry::split_key(seed.key, bucket_bits);
					let position = u32_from(flat_start + u64::from(seed.position));
					let packed =
						entry::pack(entry::prefix(tag, seed.orientation), seed.offset, position);
					if let Some(cursor) = cursors.get(bucket) {
						let slot = cursor.fetch_add(1, Ordering::Relaxed) as usize;
						if let Some(target) = slots.get(slot) {
							target.store(packed, Ordering::Relaxed);
						}
					}
				}
			},
		);
	});
	drop(cursors);
	stopped()?;
	let mut entries: Vec<u64> = slots.into_iter().map(AtomicU64::into_inner).collect();

	let histogram = workers.install(|| sort_buckets(&mut entries, &buckets));
	stopped()?;
	let filter_cutoff = filter_cutoff(&histogram);
	Ok(Index {
		config,
		bucket_bits,
		filter_cutoff,
		buckets,
		entries,
		reference,
	})
}

fn blocks(reference: &PackedReference) -> Vec<Block> {
	let mut blocks = Vec::new();
	for span in reference.spans() {
		let mut start = 0;
		while start < span.length {
			let length = BLOCK_BASES.min(span.length - start);
			blocks.push(Block {
				span: *span,
				start,
				length,
			});
			start += length;
		}
	}
	blocks
}

fn estimate_seeds(reference: &PackedReference, config: &SeedConfig) -> u64 {
	let ambiguous: u64 = reference
		.ambiguous()
		.iter()
		.map(|run| run.end - run.start)
		.sum();
	reference.length().saturating_sub(ambiguous) / config.expected_spacing() as u64
}

fn sort_buckets(entries: &mut [u64], buckets: &[u32]) -> Vec<u64> {
	let bucket_count = buckets.len() - 1;
	let mut segments: Vec<(usize, &mut [u64])> = Vec::new();
	let mut rest = entries;
	let mut bucket = 0usize;
	while bucket < bucket_count {
		let first = bucket;
		let begin = u64::from(buckets.get(bucket).copied().unwrap_or(0));
		let mut end = bucket + 1;
		while end < bucket_count
			&& u64::from(buckets.get(end + 1).copied().unwrap_or(0)) - begin < SORT_SEGMENT_ENTRIES
		{
			end += 1;
		}
		let last = u64::from(buckets.get(end).copied().unwrap_or(0));
		let (head, tail) = rest.split_at_mut(usize_of(last - begin));
		rest = tail;
		segments.push((first, head));
		bucket = end;
	}
	segments
		.par_iter_mut()
		.map(|(first, segment)| sort_segment(*first, segment, buckets))
		.reduce(
			|| vec![0u64; HISTOGRAM_LEN],
			|mut left, right| {
				for (total, count) in left.iter_mut().zip(right) {
					*total += count;
				}
				left
			},
		)
}

fn sort_segment(first_bucket: usize, segment: &mut [u64], buckets: &[u32]) -> Vec<u64> {
	let mut histogram = vec![0u64; HISTOGRAM_LEN];
	let base = buckets.get(first_bucket).copied().unwrap_or(0) as usize;
	let mut bucket = first_bucket;
	let mut offset = 0usize;
	while offset < segment.len() {
		let end = (buckets.get(bucket + 1).copied().unwrap_or(0) as usize).saturating_sub(base);
		if let Some(slice) = segment.get_mut(offset..end) {
			slice.sort_unstable();
			count_runs(slice, &mut histogram);
		}
		offset = end;
		bucket += 1;
	}
	histogram
}

fn count_runs(sorted: &[u64], histogram: &mut [u64]) {
	for run in sorted.chunk_by(|a, b| entry::prefix_of(*a) == entry::prefix_of(*b)) {
		if let Some(slot) = histogram.get_mut(run.len().min(HISTOGRAM_LEN - 1)) {
			*slot += 1;
		}
	}
}

fn filter_cutoff(histogram: &[u64]) -> u32 {
	let distinct: u64 = histogram.iter().sum();
	let wanted = distinct / FILTER_DENOMINATOR;
	let mut total = 0;
	let mut cutoff = 1u32;
	for (count, &seeds) in histogram.iter().enumerate().skip(1).rev() {
		total += seeds;
		if total >= wanted {
			cutoff = u32_of(count);
			break;
		}
	}
	cutoff.clamp(FILTER_MIN, FILTER_MAX)
}
