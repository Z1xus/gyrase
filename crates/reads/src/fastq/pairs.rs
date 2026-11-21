use std::io::BufRead;

use super::batch::{FastqBatch, ReadRef};
use super::reader::FastqReader;
use crate::error::{Error, FastqError};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadPair<'a> {
	pub first: ReadRef<'a>,
	pub second: ReadRef<'a>,
}

#[derive(Debug, Default)]
pub struct PairedBatch {
	first: FastqBatch,
	second: FastqBatch,
}

impl PairedBatch {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn len(&self) -> usize {
		self.first.len()
	}

	pub fn is_empty(&self) -> bool {
		self.first.is_empty()
	}

	pub fn pairs(&self) -> impl Iterator<Item = ReadPair<'_>> {
		self.first
			.iter()
			.zip(self.second.iter())
			.map(|(first, second)| ReadPair { first, second })
	}
}

#[derive(Debug)]
pub struct PairedFastq<R1: BufRead, R2: BufRead> {
	first: FastqReader<R1>,
	second: FastqReader<R2>,
	pairs_read: u64,
}

impl<R1: BufRead, R2: BufRead> PairedFastq<R1, R2> {
	pub fn new(first: R1, second: R2) -> Self {
		Self {
			first: FastqReader::new(first),
			second: FastqReader::new(second),
			pairs_read: 0,
		}
	}

	pub fn read_batch(
		&mut self,
		batch: &mut PairedBatch,
		max_pairs: usize,
	) -> Result<usize, Error> {
		let first = self.first.read_records(&mut batch.first, max_pairs)?;
		let second = self.second.read_records(&mut batch.second, max_pairs)?;
		if first != second {
			return Err(FastqError::PairCount.into());
		}
		check_names(batch.pairs(), self.pairs_read)?;
		self.pairs_read += u64::try_from(first).unwrap_or(0);
		Ok(first)
	}
}

#[derive(Debug)]
pub struct InterleavedFastq<R: BufRead> {
	reader: FastqReader<R>,
	pairs_read: u64,
}

impl<R: BufRead> InterleavedFastq<R> {
	pub fn new(source: R) -> Self {
		Self {
			reader: FastqReader::new(source),
			pairs_read: 0,
		}
	}

	pub fn read_batch(&mut self, batch: &mut FastqBatch, max_pairs: usize) -> Result<usize, Error> {
		let mut reads = self.reader.read_batch(batch, max_pairs.saturating_mul(2))?;
		if reads % 2 != 0 {
			reads = self.reader.read_more(batch, reads + 1)?;
		}
		if reads % 2 != 0 {
			return Err(FastqError::PairCount.into());
		}
		check_names(interleaved_pairs(batch), self.pairs_read)?;
		let pairs = reads / 2;
		self.pairs_read += u64::try_from(pairs).unwrap_or(0);
		Ok(pairs)
	}
}

pub fn interleaved_pairs(batch: &FastqBatch) -> impl Iterator<Item = ReadPair<'_>> {
	(0..batch.len() / 2).filter_map(|index| {
		Some(ReadPair {
			first: batch.get(index * 2)?,
			second: batch.get(index * 2 + 1)?,
		})
	})
}

pub fn interleaved(batch: &FastqBatch) -> bool {
	match (batch.get(0), batch.get(1)) {
		(Some(first), Some(second)) => mate_stem(first.name) == mate_stem(second.name),
		_ => false,
	}
}

fn check_names<'a>(
	pairs: impl Iterator<Item = ReadPair<'a>>,
	already_read: u64,
) -> Result<(), FastqError> {
	for (index, pair) in (1..).zip(pairs) {
		if mate_stem(pair.first.name) != mate_stem(pair.second.name) {
			return Err(FastqError::PairNames {
				record: already_read + index,
			});
		}
	}
	Ok(())
}

fn mate_stem(name: &[u8]) -> &[u8] {
	name.strip_suffix(b"/1")
		.or_else(|| name.strip_suffix(b"/2"))
		.unwrap_or(name)
}
