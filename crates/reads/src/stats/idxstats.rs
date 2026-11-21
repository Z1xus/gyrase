use std::fmt::Write;

use gyrase_genome::ContigId;

use crate::bam::BamReader;
use crate::error::{Error, RecordError, StatsError};
use crate::header::Header;
use crate::record::RecordRef;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Last {
	Start,
	At(Option<ContigId>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Idxstats {
	placed: Vec<[u64; 2]>,
	unplaced: [u64; 2],
	last: Last,
}

impl Idxstats {
	pub fn new(reference_count: usize) -> Self {
		Self {
			placed: vec![[0; 2]; reference_count],
			unplaced: [0; 2],
			last: Last::Start,
		}
	}

	pub fn add(&mut self, record: RecordRef<'_>) -> Result<(), Error> {
		let contig = record.contig();
		let class = usize::from(record.flags().is_unmapped());
		let counts = match contig {
			Some(id) => self
				.placed
				.get_mut(id.index())
				.ok_or(RecordError::ContigOutOfRange { contig: id.get() })?,
			None => &mut self.unplaced,
		};
		if self.last != Last::At(contig) && self.last != Last::Start && counts[0] + counts[1] > 0 {
			return Err(StatsError::NotSorted.into());
		}
		if let Some(slot) = counts.get_mut(class) {
			*slot += 1;
		}
		self.last = Last::At(contig);
		Ok(())
	}

	pub fn render(&self, header: &Header) -> String {
		let mut out = String::new();
		for (reference, counts) in header.references().iter().zip(&self.placed) {
			let _ = writeln!(
				out,
				"{}\t{}\t{}\t{}",
				reference.name(),
				reference.length(),
				counts[0],
				counts[1]
			);
		}
		let _ = writeln!(out, "*\t0\t{}\t{}", self.unplaced[0], self.unplaced[1]);
		out
	}
}

pub fn idxstats(reader: &mut BamReader) -> Result<Idxstats, Error> {
	let mut stats = Idxstats::new(reader.header().references().len());
	while let Some(located) = reader.next_record()? {
		stats.add(located.record)?;
	}
	Ok(stats)
}
