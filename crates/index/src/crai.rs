use std::io::{Read, Write};

use gyrase_genome::{ContigId, Interval};
use gyrase_io::{BgzfWriter, Level, Workers};

use crate::Error;
use crate::wire::read_all;

const UNMAPPED: i64 = -1;
const MULTI_REFERENCE: i64 = -2;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CraiReference {
	Unmapped,
	Multi,
	Contig(ContigId),
}

impl CraiReference {
	fn from_code(code: i64) -> Option<Self> {
		match code {
			UNMAPPED => Some(Self::Unmapped),
			MULTI_REFERENCE => Some(Self::Multi),
			_ => u32::try_from(code)
				.ok()
				.map(|id| Self::Contig(ContigId::new(id))),
		}
	}

	fn code(self) -> i64 {
		match self {
			Self::Unmapped => UNMAPPED,
			Self::Multi => MULTI_REFERENCE,
			Self::Contig(contig) => i64::from(contig.get()),
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CraiEntry {
	pub reference: CraiReference,
	pub start: u64,
	pub span: u64,
	pub container_offset: u64,
	pub slice_offset: u64,
	pub slice_size: u64,
}

impl CraiEntry {
	fn parse(line: &[u8]) -> Option<Self> {
		let mut fields = std::str::from_utf8(line).ok()?.split('\t');
		let entry = Self {
			reference: CraiReference::from_code(fields.next()?.parse().ok()?)?,
			start: fields.next()?.parse().ok()?,
			span: fields.next()?.parse().ok()?,
			container_offset: fields.next()?.parse().ok()?,
			slice_offset: fields.next()?.parse().ok()?,
			slice_size: fields.next()?.parse().ok()?,
		};
		fields.next().is_none().then_some(entry)
	}

	fn overlaps(&self, interval: Interval) -> bool {
		match self.reference {
			CraiReference::Unmapped => false,
			CraiReference::Multi => true,
			CraiReference::Contig(contig) => {
				let begin = self.start.saturating_sub(1);
				contig == interval.contig
					&& begin < interval.end
					&& begin.saturating_add(self.span.max(1)) > interval.start
			}
		}
	}
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Crai {
	entries: Vec<CraiEntry>,
}

impl Crai {
	pub fn new(entries: Vec<CraiEntry>) -> Self {
		Self { entries }
	}

	pub fn entries(&self) -> &[CraiEntry] {
		&self.entries
	}

	pub fn overlapping(&self, interval: Interval) -> impl Iterator<Item = &CraiEntry> {
		self.entries
			.iter()
			.filter(move |entry| entry.overlaps(interval))
	}

	pub fn read(input: impl Read + Send + 'static, workers: &Workers) -> Result<Self, Error> {
		let bytes = read_all(gyrase_io::decode(Box::new(input), workers)?)?;
		let mut entries = Vec::new();
		for (number, line) in bytes.split(|&byte| byte == b'\n').enumerate() {
			let line = line.strip_suffix(b"\r").unwrap_or(line);
			if line.is_empty() {
				continue;
			}
			entries.push(CraiEntry::parse(line).ok_or(Error::InvalidCrai { line: number + 1 })?);
		}
		Ok(Self { entries })
	}

	pub fn write<W: Write + Send + 'static>(&self, out: W, workers: &Workers) -> Result<W, Error> {
		let mut writer = BgzfWriter::new(out, Level::DEFAULT, workers);
		for entry in &self.entries {
			writeln!(
				writer,
				"{}\t{}\t{}\t{}\t{}\t{}",
				entry.reference.code(),
				entry.start,
				entry.span,
				entry.container_offset,
				entry.slice_offset,
				entry.slice_size,
			)?;
		}
		Ok(writer.finish()?.inner)
	}
}
