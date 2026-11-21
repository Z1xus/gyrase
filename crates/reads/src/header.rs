mod bam;
mod fields;
mod reference;
mod text;

use std::collections::HashMap;

use gyrase_genome::{Contig, ContigId, Contigs};

pub use fields::{Field, Fields};
pub use reference::Reference;

pub(crate) use bam::ByteSource;

use crate::error::{Error, HeaderError};
use reference::{MAX_REFERENCE_LEN, md5_hex};
use text::{ParsedLine, parse_line, write_fields};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Line {
	Header(Fields),
	Reference(ContigId),
	ReadGroup(Fields),
	Program(Fields),
	Comment(String),
	Other(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
	Unknown,
	Unsorted,
	QueryName,
	Coordinate,
}

impl SortOrder {
	const fn as_str(self) -> &'static str {
		match self {
			Self::Unknown => "unknown",
			Self::Unsorted => "unsorted",
			Self::QueryName => "queryname",
			Self::Coordinate => "coordinate",
		}
	}
}

const DEFAULT_VERSION: &str = "1.6";

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Header {
	lines: Vec<Line>,
	references: Vec<Reference>,
	by_name: HashMap<String, ContigId>,
}

impl Header {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn parse(text: &[u8]) -> Result<Self, HeaderError> {
		let text = std::str::from_utf8(text).map_err(|_| HeaderError::InvalidUtf8)?;
		let mut header = Self::new();
		for (index, line) in text.trim_end_matches('\n').lines().enumerate() {
			let number = index + 1;
			match parse_line(line, number)? {
				ParsedLine::Header(fields) => header.lines.push(Line::Header(fields)),
				ParsedLine::Reference(fields) => {
					header.push_reference(Reference::from_fields(fields, number)?)?;
				}
				ParsedLine::ReadGroup(fields) => header.lines.push(Line::ReadGroup(fields)),
				ParsedLine::Program(fields) => header.lines.push(Line::Program(fields)),
				ParsedLine::Comment(comment) => {
					header.lines.push(Line::Comment(comment.to_owned()));
				}
				ParsedLine::Other(raw) => header.lines.push(Line::Other(raw.to_owned())),
			}
		}
		Ok(header)
	}

	pub fn from_contigs(contigs: &Contigs) -> Result<Self, HeaderError> {
		let mut header = Self::new();
		for (_, contig) in contigs.iter() {
			let length = u32::try_from(contig.length)
				.ok()
				.filter(|&length| length <= MAX_REFERENCE_LEN)
				.ok_or_else(|| HeaderError::ReferenceTooLong(contig.name.clone()))?;
			let mut reference = Reference::new(contig.name.clone(), length)?;
			if let Some(checksum) = &contig.checksum {
				reference = reference.with_field(*b"M5", md5_hex(checksum));
			}
			header.push_reference(reference)?;
		}
		Ok(header)
	}

	pub fn push_header_line(&mut self, fields: Fields) {
		self.lines.push(Line::Header(fields));
	}

	pub fn set_sort_order(&mut self, order: SortOrder) {
		let existing = self.lines.iter_mut().find_map(|line| match line {
			Line::Header(fields) => Some(fields),
			_ => None,
		});
		match existing {
			Some(fields) => {
				fields.set(*b"SO", order.as_str());
				fields.remove(*b"GO");
				fields.remove(*b"SS");
			}
			None => self.lines.insert(
				0,
				Line::Header(
					Fields::new()
						.with(*b"VN", DEFAULT_VERSION)
						.with(*b"SO", order.as_str()),
				),
			),
		}
	}

	pub fn push_reference(&mut self, reference: Reference) -> Result<ContigId, HeaderError> {
		let id = u32::try_from(self.references.len())
			.map(ContigId::new)
			.map_err(|_| HeaderError::TooLarge { limit: usize::MAX })?;
		if self.by_name.contains_key(reference.name()) {
			return Err(HeaderError::DuplicateReference(reference.name().to_owned()));
		}
		self.by_name.insert(reference.name().to_owned(), id);
		self.references.push(reference);
		self.lines.push(Line::Reference(id));
		Ok(id)
	}

	pub fn push_read_group(&mut self, fields: Fields) {
		self.lines.push(Line::ReadGroup(fields));
	}

	pub fn push_program(&mut self, fields: Fields) {
		self.lines.push(Line::Program(fields));
	}

	pub fn push_comment(&mut self, comment: impl Into<String>) {
		self.lines.push(Line::Comment(comment.into()));
	}

	pub fn lines(&self) -> &[Line] {
		&self.lines
	}

	pub fn references(&self) -> &[Reference] {
		&self.references
	}

	pub fn reference(&self, id: ContigId) -> Option<&Reference> {
		self.references.get(id.index())
	}

	pub fn contig_id(&self, name: &str) -> Option<ContigId> {
		self.by_name.get(name).copied()
	}

	pub fn text(&self) -> String {
		let mut out = String::new();
		for line in &self.lines {
			match line {
				Line::Header(fields) => write_fields(&mut out, "HD", fields),
				Line::Reference(id) => {
					if let Some(reference) = self.reference(*id) {
						write_fields(&mut out, "SQ", reference.fields());
					}
				}
				Line::ReadGroup(fields) => write_fields(&mut out, "RG", fields),
				Line::Program(fields) => write_fields(&mut out, "PG", fields),
				Line::Comment(comment) => {
					out.push_str("@CO\t");
					out.push_str(comment);
					out.push('\n');
				}
				Line::Other(raw) => {
					out.push_str(raw);
					out.push('\n');
				}
			}
		}
		out
	}

	pub fn contigs(&self) -> Result<Contigs, gyrase_genome::Error> {
		Contigs::new(self.references.iter().map(|reference| Contig {
			name: reference.name().to_owned(),
			length: u64::from(reference.length()),
			checksum: reference.md5(),
		}))
	}

	pub fn write_bam(&self, out: &mut Vec<u8>) -> Result<(), HeaderError> {
		bam::write(self, out)
	}

	pub(crate) fn read_bam(source: &mut impl ByteSource) -> Result<Self, Error> {
		bam::read(source)
	}

	fn insert_references(&mut self, references: Vec<Reference>) -> Result<(), HeaderError> {
		let after_header = usize::from(matches!(self.lines.first(), Some(Line::Header(_))));
		let mut inserted = Vec::with_capacity(references.len());
		for reference in references {
			let id = ContigId::new(
				u32::try_from(self.references.len())
					.map_err(|_| HeaderError::TooLarge { limit: usize::MAX })?,
			);
			if self
				.by_name
				.insert(reference.name().to_owned(), id)
				.is_some()
			{
				return Err(HeaderError::DuplicateReference(reference.name().to_owned()));
			}
			self.references.push(reference);
			inserted.push(Line::Reference(id));
		}
		self.lines.splice(after_header..after_header, inserted);
		Ok(())
	}
}
