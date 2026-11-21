use std::cmp::Reverse;
use std::collections::BinaryHeap;

use gyrase_io::Workers;

use super::buffer::Buffer;
use super::run::{RunFile, RunReader};
use crate::error::Error;

enum Source<'a> {
	Memory { buffer: &'a Buffer, at: usize },
	File(RunReader),
}

impl Source<'_> {
	fn key(&self) -> u64 {
		match self {
			Self::Memory { buffer, at } => buffer.key_at(*at).unwrap_or(u64::MAX),
			Self::File(reader) => reader.key(),
		}
	}

	fn record(&self) -> &[u8] {
		match self {
			Self::Memory { buffer, at } => buffer.record_at(*at).unwrap_or_default(),
			Self::File(reader) => reader.record(),
		}
	}

	fn advance(&mut self) -> Result<bool, Error> {
		match self {
			Self::Memory { buffer, at } => {
				*at += 1;
				Ok(*at < buffer.len())
			}
			Self::File(reader) => reader.advance(),
		}
	}
}

pub(super) fn merge(
	buffer: &Buffer,
	runs: &[RunFile],
	workers: &Workers,
	mut emit: impl FnMut(&[u8]) -> Result<(), Error>,
) -> Result<(), Error> {
	let mut sources = Vec::with_capacity(runs.len() + 1);
	let mut heap = BinaryHeap::new();
	for run in runs {
		let mut reader = run.open(workers)?;
		if reader.advance()? {
			heap.push(Reverse((reader.key(), sources.len())));
			sources.push(Source::File(reader));
		}
	}
	if buffer.len() > 0 {
		let source = Source::Memory { buffer, at: 0 };
		heap.push(Reverse((source.key(), sources.len())));
		sources.push(source);
	}
	while let Some(Reverse((_, index))) = heap.pop() {
		let Some(source) = sources.get_mut(index) else {
			continue;
		};
		emit(source.record())?;
		if source.advance()? {
			heap.push(Reverse((source.key(), index)));
		}
	}
	Ok(())
}
