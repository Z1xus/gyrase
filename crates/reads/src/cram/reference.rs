use std::sync::Arc;

use gyrase_genome::{ContigId, Interval, Reference};

use crate::Header;
use crate::cram::Error;

#[derive(Debug)]
struct Entry {
	name: String,
	length: u64,
	id: Option<ContigId>,
}

#[derive(Debug)]
pub(crate) struct References {
	reference: Option<Arc<Reference>>,
	entries: Vec<Entry>,
}

impl References {
	pub(crate) fn new(header: &Header, reference: Option<Arc<Reference>>) -> Self {
		let entries = header
			.references()
			.iter()
			.map(|entry| Entry {
				name: entry.name().to_owned(),
				length: u64::from(entry.length()),
				id: reference
					.as_ref()
					.and_then(|source| source.contigs().resolve(entry.name())),
			})
			.collect();
		Self { reference, entries }
	}

	pub(crate) fn count(&self) -> usize {
		self.entries.len()
	}

	pub(crate) fn length(&self, index: usize) -> Option<u64> {
		self.entries.get(index).map(|entry| entry.length)
	}

	pub(crate) fn has_source(&self) -> bool {
		self.reference.is_some()
	}

	fn interval(
		&self,
		index: usize,
		start: u64,
		end: u64,
	) -> Result<(Interval, &Reference), Error> {
		let source = self.reference.as_deref().ok_or(Error::ReferenceRequired)?;
		let entry = self
			.entries
			.get(index)
			.ok_or(Error::Malformed("reference id"))?;
		let contig = entry
			.id
			.ok_or_else(|| Error::UnknownReference(entry.name.clone()))?;
		let end = end.min(entry.length);
		let start = start.min(end);
		Ok((Interval { contig, start, end }, source))
	}

	pub(crate) fn window(&self, index: usize, start: u64, end: u64) -> Result<Window, Error> {
		let (interval, source) = self.interval(index, start, end)?;
		let mut bases = Vec::new();
		source.fetch(&interval, &mut bases)?;
		Ok(Window {
			contig_length: self.length(index).unwrap_or(0),
			start: interval.start,
			bases,
		})
	}

	pub(crate) fn checksum(&self, index: usize, start: u64, end: u64) -> Result<[u8; 16], Error> {
		let (interval, source) = self.interval(index, start, end)?;
		Ok(source.checksum(&interval)?)
	}
}

#[derive(Debug, Default)]
pub(crate) struct Window {
	contig_length: u64,
	start: u64,
	bases: Vec<u8>,
}

impl Window {
	pub(crate) fn embedded(contig_length: u64, start: u64, bases: Vec<u8>) -> Self {
		Self {
			contig_length,
			start,
			bases,
		}
	}

	pub(crate) fn at(&self, position: u64) -> Result<u8, Error> {
		if position >= self.contig_length {
			return Ok(b'N');
		}
		let offset = position
			.checked_sub(self.start)
			.and_then(|offset| usize::try_from(offset).ok())
			.ok_or(Error::Malformed("reference position"))?;
		self.bases
			.get(offset)
			.copied()
			.ok_or(Error::Malformed("reference position"))
	}

	pub(crate) fn copy(&self, position: u64, out: &mut [u8]) -> Result<(), Error> {
		let wanted = out.len() as u64;
		let inside = self.contig_length.saturating_sub(position).min(wanted);
		let (head, tail) =
			out.split_at_mut(usize::try_from(inside).map_err(|_| Error::IntegerOverflow)?);
		if !head.is_empty() {
			let offset = position
				.checked_sub(self.start)
				.and_then(|offset| usize::try_from(offset).ok())
				.ok_or(Error::Malformed("reference position"))?;
			let source = self
				.bases
				.get(offset..offset + head.len())
				.ok_or(Error::Malformed("reference position"))?;
			head.copy_from_slice(source);
		}
		tail.fill(b'N');
		Ok(())
	}
}
