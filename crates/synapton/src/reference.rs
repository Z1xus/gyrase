use std::ops::Range;

use gyrase_genome::{ContigId, Contigs, Interval, Reference};
use gyrase_io::Workers;
use rayon::prelude::*;

use crate::Error;
use crate::hash::{checksum_words, mix};
use crate::nucleotide::{CODE_N, CODE_PAD, encode};
use crate::num::{i64_len, i64_of, len_u64, u32_of, u64_of, usize_of};

pub(crate) const WORD_BASES: u64 = 32;
const BITS_PER_BASE: u64 = 2;
const PACK_TILE: u64 = 1 << 23;
pub(crate) const MAX_FLAT_LENGTH: u64 = u32::MAX as u64 - 4 * WORD_BASES;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ContigSpan {
	pub(crate) start: u64,
	pub(crate) length: u64,
}

impl ContigSpan {
	pub(crate) const fn end(&self) -> u64 {
		self.start + self.length
	}
}

#[derive(Debug)]
pub struct PackedReference {
	words: Vec<u64>,
	length: u64,
	ambiguous: Vec<Range<u64>>,
	spans: Vec<ContigSpan>,
	contigs: Contigs,
	checksum: u64,
}

struct PackTile {
	contig: ContigId,
	offset: u64,
	length: u64,
	flat_start: u64,
}

struct PackedTile {
	first_word: usize,
	words: Vec<u64>,
	ambiguous: Vec<Range<u64>>,
}

const fn padded(length: u64) -> u64 {
	length.div_ceil(WORD_BASES) * WORD_BASES
}

impl PackedReference {
	pub fn from_fasta(reference: &Reference, workers: &Workers) -> Result<Self, Error> {
		let mut spans = Vec::new();
		let mut tiles = Vec::new();
		let mut flat = 0u64;
		for (contig_id, contig) in reference.contigs().iter() {
			spans.push(ContigSpan {
				start: flat,
				length: contig.length,
			});
			let mut offset = 0;
			while offset < contig.length {
				let length = PACK_TILE.min(contig.length - offset);
				tiles.push(PackTile {
					contig: contig_id,
					offset,
					length,
					flat_start: flat + offset,
				});
				offset += length;
			}
			flat += padded(contig.length);
		}
		if flat > MAX_FLAT_LENGTH {
			return Err(Error::ReferenceTooLong {
				length: flat,
				limit: MAX_FLAT_LENGTH,
			});
		}
		let packed: Vec<Result<PackedTile, Error>> = workers.install(|| {
			tiles
				.par_iter()
				.map_init(Vec::new, |buffer, tile| pack_tile(reference, tile, buffer))
				.collect()
		});
		let mut words = vec![0u64; usize::try_from(flat / WORD_BASES)? + 1];
		let mut ambiguous: Vec<Range<u64>> = Vec::new();
		for tile in packed {
			let tile = tile?;
			let target = words
				.get_mut(tile.first_word..tile.first_word + tile.words.len())
				.ok_or(Error::CorruptIndex("tile outside reference"))?;
			target.copy_from_slice(&tile.words);
			for run in tile.ambiguous {
				push_run(&mut ambiguous, run);
			}
		}
		for span in &spans {
			push_run(&mut ambiguous, span.end()..span.start + padded(span.length));
		}
		ambiguous.retain(|run| run.start < run.end);
		Self::from_parts(words, flat, ambiguous, spans, reference.contigs().clone())
	}

	pub(crate) fn from_parts(
		words: Vec<u64>,
		length: u64,
		ambiguous: Vec<Range<u64>>,
		spans: Vec<ContigSpan>,
		contigs: Contigs,
	) -> Result<Self, Error> {
		if spans.len() != contigs.len() {
			return Err(Error::CorruptIndex("contig count"));
		}
		let needed = usize::try_from(length / WORD_BASES + 1)?;
		if words.len() != needed {
			return Err(Error::CorruptIndex("reference size"));
		}
		let mut checksum = checksum_words(&words);
		for (span, (_, contig)) in spans.iter().zip(contigs.iter()) {
			checksum = mix(checksum ^ span.start) ^ mix(span.length);
			for &byte in contig.name.as_bytes() {
				checksum = mix(checksum ^ u64::from(byte));
			}
		}
		for run in &ambiguous {
			checksum = mix(checksum ^ run.start) ^ mix(run.end);
		}
		Ok(Self {
			words,
			length,
			ambiguous,
			spans,
			contigs,
			checksum,
		})
	}

	pub fn contigs(&self) -> &Contigs {
		&self.contigs
	}

	pub fn checksum(&self) -> u64 {
		self.checksum
	}

	pub(crate) fn length(&self) -> u64 {
		self.length
	}

	pub(crate) fn words(&self) -> &[u64] {
		&self.words
	}

	pub(crate) fn ambiguous(&self) -> &[Range<u64>] {
		&self.ambiguous
	}

	pub(crate) fn spans(&self) -> &[ContigSpan] {
		&self.spans
	}

	pub(crate) fn span(&self, contig: ContigId) -> Option<&ContigSpan> {
		self.spans.get(contig.index())
	}

	pub(crate) fn locate(&self, position: u64) -> Option<(ContigId, &ContigSpan)> {
		let after = self.spans.partition_point(|span| span.start <= position);
		let index = after.checked_sub(1)?;
		let span = self.spans.get(index)?;
		(position < span.end()).then(|| (ContigId::new(u32_of(index)), span))
	}

	pub(crate) fn window(&self, span: &ContigSpan, start: i64, length: usize, out: &mut Vec<u8>) {
		out.clear();
		let end = start + i64_len(length);
		let inner_start = start.max(i64_of(span.start));
		let inner_end = end.min(i64_of(span.end()));
		if inner_start >= inner_end {
			out.resize(length, CODE_PAD);
			return;
		}
		out.resize(usize_of(u64_of(inner_start - start)), CODE_PAD);
		let (first, last) = (u64_of(inner_start), u64_of(inner_end));
		out.extend((first..last).map(|position| self.base_code(position)));
		out.resize(length, CODE_PAD);
		let mut run = self.ambiguous.partition_point(|run| run.end <= first);
		while let Some(ambiguous) = self.ambiguous.get(run) {
			if ambiguous.start >= last {
				break;
			}
			let from = usize_of(u64_of(i64_of(ambiguous.start.max(first)) - start));
			let to = usize_of(u64_of(i64_of(ambiguous.end.min(last)) - start));
			if let Some(target) = out.get_mut(from..to) {
				target.fill(CODE_N);
			}
			run += 1;
		}
	}

	fn base_code(&self, position: u64) -> u8 {
		let word = self
			.words
			.get(usize::try_from(position / WORD_BASES).unwrap_or(usize::MAX))
			.copied()
			.unwrap_or(0);
		((word >> ((position % WORD_BASES) * BITS_PER_BASE)) & 3) as u8
	}
}

fn push_run(runs: &mut Vec<Range<u64>>, run: Range<u64>) {
	if run.start >= run.end {
		return;
	}
	match runs.last_mut() {
		Some(last) if last.end >= run.start => last.end = last.end.max(run.end),
		_ => runs.push(run),
	}
}

fn pack_tile(
	reference: &Reference,
	tile: &PackTile,
	buffer: &mut Vec<u8>,
) -> Result<PackedTile, Error> {
	reference.fetch(
		&Interval {
			contig: tile.contig,
			start: tile.offset,
			end: tile.offset + tile.length,
		},
		buffer,
	)?;
	let mut codes = Vec::with_capacity(buffer.len());
	encode(buffer, &mut codes);
	let mut words = vec![0u64; usize::try_from(tile.length.div_ceil(WORD_BASES))?];
	let mut ambiguous = Vec::new();
	for (index, &code) in codes.iter().enumerate() {
		let position = tile.flat_start + len_u64(index);
		if code > 3 {
			push_run(&mut ambiguous, position..position + 1);
		} else if let Some(word) = words.get_mut(usize_of(len_u64(index) / WORD_BASES)) {
			*word |= u64::from(code) << ((len_u64(index) % WORD_BASES) * BITS_PER_BASE);
		}
	}
	Ok(PackedTile {
		first_word: usize::try_from(tile.flat_start / WORD_BASES)?,
		words,
		ambiguous,
	})
}
