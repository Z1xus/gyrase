use memchr::{memchr, memchr_iter};

use crate::Error;
use crate::header::Header;
use crate::record::Record;
use crate::site::{Alleles, Site};
use crate::value::Loaded;

use super::column::{kind_of, load_column};
use super::text::parse_position;

const FIXED_COLUMNS: usize = 8;
const CHROM: usize = 0;
const POS: usize = 1;
const ID: usize = 2;
const REF: usize = 3;
const ALT: usize = 4;
const QUAL: usize = 5;
const FILTER: usize = 6;
const INFO: usize = 7;
const FORMAT: usize = 8;

#[derive(Debug, Clone, Copy)]
pub struct VcfLine<'a> {
	text: &'a [u8],
	tabs: [usize; FORMAT + 1],
	tab_count: usize,
}

impl<'a> VcfLine<'a> {
	pub fn new(text: &'a [u8]) -> Result<Self, Error> {
		let mut tabs = [0; FORMAT + 1];
		let mut tab_count = 0;
		for (slot, tab) in tabs.iter_mut().zip(memchr_iter(b'\t', text)) {
			*slot = tab;
			tab_count += 1;
		}
		if tab_count < FIXED_COLUMNS - 1 {
			return Err(Error::Vcf("record has fewer than eight columns"));
		}
		Ok(Self {
			text,
			tabs,
			tab_count,
		})
	}

	pub fn text(&self) -> &'a [u8] {
		self.text
	}

	fn column(&self, index: usize) -> &'a [u8] {
		let start = match index.checked_sub(1) {
			Some(previous) => self
				.tabs
				.get(previous)
				.map_or(self.text.len(), |tab| tab + 1),
			None => 0,
		};
		let end = if index < self.tab_count {
			self.tabs.get(index).copied().unwrap_or(self.text.len())
		} else {
			self.text.len()
		};
		self.text.get(start..end).unwrap_or_default()
	}

	pub fn chrom(&self) -> &'a [u8] {
		self.column(CHROM)
	}

	pub fn position(&self) -> Result<u32, Error> {
		parse_position(self.column(POS))
	}

	pub fn reference(&self) -> &'a [u8] {
		self.column(REF)
	}

	pub fn alt_text(&self) -> &'a [u8] {
		self.column(ALT)
	}

	pub fn qual_text(&self) -> &'a [u8] {
		self.column(QUAL)
	}

	pub fn filter_text(&self) -> &'a [u8] {
		self.column(FILTER)
	}

	pub fn info_text(&self) -> &'a [u8] {
		self.column(INFO)
	}

	pub fn has_samples(&self) -> bool {
		self.tab_count > FORMAT
	}

	pub fn format_text(&self) -> &'a [u8] {
		if self.has_samples() {
			self.column(FORMAT)
		} else {
			&[]
		}
	}

	pub fn sample_text(&self) -> &'a [u8] {
		match self.tabs.get(FORMAT) {
			Some(&tab) if self.has_samples() => self.text.get(tab + 1..).unwrap_or_default(),
			_ => &[],
		}
	}

	pub fn prefix_before_samples(&self) -> &'a [u8] {
		match self.tabs.get(FORMAT) {
			Some(&tab) if self.has_samples() => self.text.get(..tab).unwrap_or_default(),
			_ => self.text,
		}
	}

	pub fn samples(&self) -> impl Iterator<Item = &'a [u8]> {
		let rest = self.sample_text();
		let present = self.has_samples();
		present
			.then_some(rest)
			.into_iter()
			.flat_map(|rest| rest.split(|byte| *byte == b'\t'))
	}

	pub fn info_value(&self, name: &str) -> Option<InfoText<'a>> {
		let info = self.info_text();
		let wanted = name.as_bytes();
		let mut rest = info;
		loop {
			let end = memchr(b';', rest).unwrap_or(rest.len());
			let entry = rest.get(..end).unwrap_or_default();
			match entry.strip_prefix(wanted) {
				Some([]) => return Some(InfoText::Bare),
				Some([b'=', value @ ..]) => return Some(InfoText::Value(value)),
				_ => {}
			}
			rest = rest.get(end + 1..)?;
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InfoText<'a> {
	Bare,
	Value(&'a [u8]),
}

impl Site for VcfLine<'_> {
	fn contig(&self, header: &Header) -> Result<u32, Error> {
		let name =
			std::str::from_utf8(self.chrom()).map_err(|_| Error::Vcf("contig is not utf-8"))?;
		header
			.contig_id(name)
			.ok_or_else(|| Error::UndefinedContig(name.to_owned()))
	}

	fn start(&self) -> Result<u32, Error> {
		self.position()?
			.checked_sub(1)
			.ok_or(Error::Vcf("position is zero"))
	}

	fn end(&self, header: &Header) -> Result<u32, Error> {
		let start = self.start()?;
		let reference_end =
			start.saturating_add(u32::try_from(self.reference().len()).unwrap_or(u32::MAX));
		let declared = header
			.key("END")
			.and_then(|_| self.info_value("END"))
			.and_then(|value| match value {
				InfoText::Value(text) => super::text::parse_leading_uint(text),
				InfoText::Bare => None,
			});
		Ok(declared.map_or(reference_end, |end| end.max(reference_end)))
	}

	fn id(&self) -> Result<&[u8], Error> {
		Ok(match self.column(ID) {
			[b'.'] => &[],
			other => other,
		})
	}

	fn alleles(&self, out: &mut Alleles) -> Result<(), Error> {
		out.clear();
		out.push(self.reference());
		match self.alt_text() {
			[b'.'] => {}
			alts => {
				for alt in alts.split(|byte| *byte == b',') {
					out.push(alt);
				}
			}
		}
		Ok(())
	}

	fn qual(&self) -> Result<Option<f32>, Error> {
		super::text::parse_qual(self.qual_text())
	}

	fn filters(&self, header: &Header, out: &mut Vec<u32>) -> Result<(), Error> {
		out.clear();
		let text = self.filter_text();
		if text == b"." {
			return Ok(());
		}
		for name in text.split(|byte| *byte == b';') {
			let name = std::str::from_utf8(name).map_err(|_| Error::Vcf("filter is not utf-8"))?;
			out.push(
				header
					.key(name)
					.ok_or_else(|| Error::UndefinedTag(name.to_owned()))?,
			);
		}
		Ok(())
	}

	fn info(&self, header: &Header, key: u32, out: &mut Loaded) -> Result<bool, Error> {
		out.clear();
		let name = header.name(key).ok_or(Error::Vcf("unknown tag id"))?;
		let value_type = header
			.type_of_info(key)
			.ok_or_else(|| Error::UndefinedTag(name.to_owned()))?;
		match self.info_value(name) {
			None => Ok(false),
			Some(InfoText::Bare) => Ok(true),
			Some(InfoText::Value(text)) => {
				let (cells, width) =
					super::text::parse_info_value(value_type, text, &mut out.store);
				out.cells = cells;
				out.width = width;
				Ok(true)
			}
		}
	}

	fn format_keys(&self, header: &Header, out: &mut Vec<u32>) -> Result<(), Error> {
		out.clear();
		if !self.has_samples() || self.format_text() == b"." {
			return Ok(());
		}
		for name in self.format_text().split(|byte| *byte == b':') {
			let name =
				std::str::from_utf8(name).map_err(|_| Error::Vcf("FORMAT key is not utf-8"))?;
			out.push(
				header
					.key(name)
					.ok_or_else(|| Error::UndefinedTag(name.to_owned()))?,
			);
		}
		Ok(())
	}

	fn format(&self, header: &Header, key: u32, out: &mut Loaded) -> Result<bool, Error> {
		out.clear();
		if !self.has_samples() {
			return Ok(false);
		}
		let name = header.name(key).ok_or(Error::Vcf("unknown tag id"))?;
		let Some(index) = self
			.format_text()
			.split(|byte| *byte == b':')
			.position(|candidate| candidate == name.as_bytes())
		else {
			return Ok(false);
		};
		let kind = kind_of(header, key).ok_or_else(|| Error::UndefinedTag(name.to_owned()))?;
		let (cells, width) = load_column(self.sample_text(), index, kind, header, &mut out.store)?;
		out.cells = cells;
		out.width = width;
		Ok(true)
	}

	fn read_record(&self, header: &Header, out: &mut Record) -> Result<(), Error> {
		out.read_vcf(self, header)
	}
}
