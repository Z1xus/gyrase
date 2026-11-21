use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use gyrase_io::Workers;
use rayon::prelude::*;

use crate::Error;
use crate::config::SeedConfig;
use crate::num::{i64_of, u32_from, usize_of};
use crate::reference::{PackedReference, WORD_BASES};
use crate::seeding::Seeder;

const BLOCK_BASES: u64 = 1 << 21;
const KEY_BASES: usize = 28;
const PASS_BASES: u64 = 1 << 27;

struct Block {
	flat: u64,
	length: u64,
	start: u64,
	span_end: u64,
}

fn blocks(reference: &PackedReference) -> Vec<Block> {
	let mut blocks = Vec::new();
	for span in reference.spans() {
		let mut start = 0;
		while start < span.length {
			let length = BLOCK_BASES.min(span.length - start);
			blocks.push(Block {
				flat: span.start + start,
				length,
				start,
				span_end: span.length,
			});
			start += length;
		}
	}
	blocks
}

pub(super) fn build(
	reference: &PackedReference,
	config: SeedConfig,
	workers: &Workers,
	stop: &AtomicBool,
) -> Result<Vec<u64>, Error> {
	let tiles = reference.words().len();
	let repeated: Vec<AtomicU64> = (0..tiles.div_ceil(64)).map(|_| AtomicU64::new(0)).collect();
	let blocks = blocks(reference);
	let passes = reference.length().div_ceil(PASS_BASES).max(1);
	let overlap = KEY_BASES as u64;
	for pass in 0..passes {
		let mut found: Vec<(u64, u32)> = workers.install(|| {
			let parts: Vec<Vec<(u64, u32)>> = blocks
				.par_iter()
				.map_init(
					|| (Seeder::new(config), Vec::new()),
					|(seeder, codes), block| {
						if stop.load(Ordering::Relaxed) {
							return Vec::new();
						}
						let wanted = (block.length + overlap).min(block.span_end - block.start);
						let span = reference.locate(block.flat).map(|(_, span)| *span);
						let Some(span) = span else {
							return Vec::new();
						};
						reference.window(&span, i64_of(block.flat), usize_of(wanted), codes);
						let positions: Vec<u32> = seeder
							.syncmer_positions(codes, u32_from(block.length))
							.collect();
						positions
							.into_iter()
							.filter_map(|position| {
								let window = codes.get(
									usize_of(u64::from(position))
										..usize_of(u64::from(position)) + KEY_BASES,
								)?;
								let key = canonical(window);
								(key % passes == pass)
									.then_some((key, u32_from(block.flat) + position))
							})
							.collect()
					},
				)
				.collect();
			parts.into_iter().flatten().collect()
		});
		if stop.load(Ordering::Relaxed) {
			return Err(Error::Cancelled);
		}
		workers.install(|| found.par_sort_unstable());
		for run in found.chunk_by(|left, right| left.0 == right.0) {
			if run.len() > 1 {
				for &(_, position) in run {
					let start = u64::from(position);
					for tile in start / WORD_BASES..=(start + overlap - 1) / WORD_BASES {
						mark(&repeated, tile);
					}
				}
			}
		}
	}
	for run in reference.ambiguous() {
		for tile in run.start / WORD_BASES..run.end.div_ceil(WORD_BASES) {
			mark(&repeated, tile);
		}
	}
	Ok(repeated
		.into_iter()
		.map(|word| !word.into_inner())
		.collect())
}

fn mark(bits: &[AtomicU64], tile: u64) {
	if let Some(word) = bits.get(usize_of(tile / 64)) {
		word.fetch_or(1 << (tile % 64), Ordering::Relaxed);
	}
}

fn canonical(window: &[u8]) -> u64 {
	let mut forward = 0u64;
	let mut reverse = 0u64;
	for &code in window {
		forward = forward << 2 | u64::from(code & 3);
		reverse = reverse >> 2 | u64::from(!code & 3) << (2 * (KEY_BASES - 1));
	}
	forward.min(reverse).wrapping_mul(0x9E37_79B9_7F4A_7C15)
}
