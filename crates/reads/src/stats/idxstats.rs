use std::fmt::Write;

use crate::bam::BamReader;
use crate::error::{Error, RecordError};
use crate::header::Header;
use crate::record::RecordRef;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Idxstats {
	placed: Vec<[u64; 2]>,
	unplaced: [u64; 2],
}

impl Idxstats {
	pub fn new(reference_count: usize) -> Self {
		Self {
			placed: vec![[0; 2]; reference_count],
			unplaced: [0; 2],
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
		if let Some(slot) = counts.get_mut(class) {
			*slot += 1;
		}
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
