use std::io::BufRead;

use crate::Error;
use crate::header::Header;
use crate::read::{BcfReader, Reader, VcfReader};
use crate::record::Record;
use crate::remap::Remap;

pub(super) type Place = (u32, u32);

#[derive(Debug)]
enum Source<R> {
	Vcf(VcfReader<R>),
	Bcf(BcfReader<R>),
}

#[derive(Debug)]
pub(super) struct Input<R> {
	source: Source<R>,
	header: Header,
	remap: Remap,
	head: Option<Record>,
}

#[derive(Debug, Default)]
pub(super) struct Taken {
	pool: Vec<Record>,
	used: usize,
}

impl Taken {
	pub(super) fn record(&self, slot: usize) -> Option<&Record> {
		self.pool.get(..self.used)?.get(slot)
	}

	pub(super) fn count(&self) -> usize {
		self.used
	}
}

fn place(record: &Record) -> Place {
	(record.contig(), record.start())
}

impl<R: BufRead> Input<R> {
	pub(super) fn new(reader: Reader<R>, output: &Header) -> Result<Self, Error> {
		let header = reader.header().clone();
		let remap = Remap::new(&header, output)?;
		let source = match reader {
			Reader::Vcf(vcf) => Source::Vcf(vcf),
			Reader::Bcf(bcf) => Source::Bcf(bcf),
		};
		let mut input = Self {
			source,
			header,
			remap,
			head: None,
		};
		input.head = input.read(Record::new())?;
		Ok(input)
	}

	pub(super) fn sample_count(&self) -> usize {
		self.header.sample_count()
	}

	pub(super) fn position(&self) -> Option<Place> {
		self.head.as_ref().map(place)
	}

	fn read(&mut self, mut record: Record) -> Result<Option<Record>, Error> {
		match &mut self.source {
			Source::Vcf(reader) => match reader.next_line()? {
				Some(line) => record.read_vcf(&line, &self.header)?,
				None => return Ok(None),
			},
			Source::Bcf(reader) => match reader.next_record()? {
				Some(source) => record.read_bcf(&source)?,
				None => return Ok(None),
			},
		}
		self.remap.apply(&mut record)?;
		Ok(Some(record))
	}

	pub(super) fn take(&mut self, at: Place, taken: &mut Taken) -> Result<(), Error> {
		taken.used = 0;
		while let Some(mut head) = self.head.take() {
			if place(&head) != at {
				self.head = Some(head);
				break;
			}
			if taken.used == taken.pool.len() {
				taken.pool.push(Record::new());
			}
			let slot = taken
				.pool
				.get_mut(taken.used)
				.ok_or(Error::Limits("records at one position"))?;
			std::mem::swap(&mut head, slot);
			taken.used += 1;
			self.head = self.read(head)?;
			if self.position().is_some_and(|next| next < at) {
				return Err(Error::Unsorted);
			}
		}
		Ok(())
	}
}
