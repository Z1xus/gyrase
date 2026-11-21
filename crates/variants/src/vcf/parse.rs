use memchr::memchr;

use crate::Error;
use crate::header::Header;
use crate::record::{Field, Record};
use crate::site::Site;

use super::column::{kind_of, load_column};
use super::line::VcfLine;
use super::text::parse_info_value;

const MAX_FORMAT_FIELDS: usize = 255;

impl Record {
	pub fn read_vcf(&mut self, line: &VcfLine<'_>, header: &Header) -> Result<(), Error> {
		self.clear();
		self.contig = line.contig(header)?;
		self.start = line.start()?;
		self.set_id(line.id()?);
		self.push_allele(line.reference());
		if line.alt_text() != b"." {
			for alt in line.alt_text().split(|byte| *byte == b',') {
				self.push_allele(alt);
			}
		}
		self.qual = line.qual()?;
		line.filters(header, &mut self.filters)?;
		self.read_info(line.info_text(), header)?;
		if header.sample_count() > 0 && line.has_samples() {
			self.read_format(line, header)?;
		}
		Ok(())
	}

	fn read_info(&mut self, text: &[u8], header: &Header) -> Result<(), Error> {
		if text == b"." {
			return Ok(());
		}
		let mut rest = text.strip_suffix(b";").unwrap_or(text);
		loop {
			let end = memchr(b';', rest).unwrap_or(rest.len());
			let entry = rest.get(..end).unwrap_or_default();
			if !entry.is_empty() {
				self.read_info_entry(entry, header)?;
			}
			match rest.get(end + 1..) {
				Some(tail) => rest = tail,
				None => return Ok(()),
			}
		}
	}

	fn read_info_entry(&mut self, entry: &[u8], header: &Header) -> Result<(), Error> {
		let (name, value) = match memchr(b'=', entry) {
			Some(at) => (entry.get(..at).unwrap_or_default(), entry.get(at + 1..)),
			None => (entry, None),
		};
		let name = std::str::from_utf8(name).map_err(|_| Error::Vcf("INFO key is not utf-8"))?;
		let key = header
			.key(name)
			.ok_or_else(|| Error::UndefinedTag(name.to_owned()))?;
		let value_type = header
			.type_of_info(key)
			.ok_or_else(|| Error::UndefinedTag(name.to_owned()))?;
		match value {
			None => self.push_info_flag(key),
			Some(text) => {
				let (cells, width) = parse_info_value(value_type, text, &mut self.store);
				self.info.push(Field { key, width, cells });
			}
		}
		Ok(())
	}

	fn read_format(&mut self, line: &VcfLine<'_>, header: &Header) -> Result<(), Error> {
		let spec = line.format_text();
		if spec == b"." {
			return Ok(());
		}
		let sample_text = line.sample_text();
		let samples = sample_text.split(|byte| *byte == b'\t').count();
		if samples != header.sample_count() {
			return Err(Error::Vcf("sample columns do not match the header"));
		}
		for (index, name) in spec.split(|byte| *byte == b':').enumerate() {
			let name =
				std::str::from_utf8(name).map_err(|_| Error::Vcf("FORMAT key is not utf-8"))?;
			let key = header
				.key(name)
				.ok_or_else(|| Error::UndefinedTag(name.to_owned()))?;
			if self.format.iter().any(|field| field.key == key) {
				continue;
			}
			if self.format.len() == MAX_FORMAT_FIELDS {
				return Err(Error::Limits("FORMAT fields"));
			}
			let kind = kind_of(header, key).ok_or_else(|| Error::UndefinedTag(name.to_owned()))?;
			let (cells, width) = load_column(sample_text, index, kind, header, &mut self.store)?;
			self.format.push(Field { key, width, cells });
		}
		Ok(())
	}
}
