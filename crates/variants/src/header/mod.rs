mod line;

use std::collections::{HashMap, HashSet};

use gyrase_genome::{Contig, Contigs};

use crate::Error;
use crate::value::{Number, ValueType};
pub use line::FieldDefinition;
use line::{Line, LineKind, Parsed, parse_line};

const PASS_LINE: &str = "##FILTER=<ID=PASS,Description=\"All filters passed\">";
const COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
const UNKNOWN_CONTIG_LENGTH: u64 = i32::MAX as u64;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Version {
	pub major: u8,
	pub minor: u8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StringEntry {
	pub name: String,
	pub info: Option<FieldDefinition>,
	pub format: Option<FieldDefinition>,
	pub filter: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigEntry {
	pub name: String,
	pub length: Option<u64>,
}

#[derive(Debug, Clone)]
pub struct Header {
	lines: Vec<Line>,
	samples: Vec<String>,
	strings: Vec<Option<StringEntry>>,
	string_ids: HashMap<String, u32>,
	contigs: Vec<ContigEntry>,
	contig_ids: HashMap<String, u32>,
	version: Version,
}

impl Header {
	pub fn parse(text: &[u8]) -> Result<Self, Error> {
		let text = std::str::from_utf8(text).map_err(|_| Error::Header("header is not utf-8"))?;
		let mut meta = Vec::new();
		let mut samples = None;
		for line in text.trim_end_matches('\0').lines() {
			if line.starts_with("#CHROM") {
				samples = Some(parse_samples(line)?);
				break;
			}
			if !line.is_empty() {
				meta.push(line);
			}
		}
		Self::build(&meta, samples.ok_or(Error::Header("missing #CHROM line"))?)
	}

	pub fn from_lines(meta: &[&str], samples: Vec<String>) -> Result<Self, Error> {
		Self::build(meta, samples)
	}

	fn build(meta: &[&str], samples: Vec<String>) -> Result<Self, Error> {
		let mut header = Self {
			lines: Vec::new(),
			samples,
			strings: Vec::new(),
			string_ids: HashMap::new(),
			contigs: Vec::new(),
			contig_ids: HashMap::new(),
			version: Version { major: 4, minor: 2 },
		};
		let mut seen = HashSet::new();
		let first_is_format = meta
			.first()
			.is_some_and(|line| line.starts_with("##fileformat="));
		if !first_is_format {
			return Err(Error::Header("first line is not ##fileformat"));
		}
		let mut parsed = meta
			.iter()
			.map(|text| parse_line(text.trim_end()))
			.collect::<Result<Vec<_>, _>>()?;
		let has_pass = parsed.iter().any(|line| {
			line.line.kind == LineKind::Filter && line.line.id.as_deref() == Some("PASS")
		});
		let first_slot_claimed = parsed.iter().any(|line| {
			matches!(
				line.line.kind,
				LineKind::Info | LineKind::Format | LineKind::Filter
			) && line.explicit_index == Some(0)
		});
		if !has_pass {
			let pass = parse_line(PASS_LINE)?;
			if first_slot_claimed {
				parsed.push(pass);
			} else {
				parsed.insert(1, pass);
			}
		}
		for line in parsed {
			header.add_line(line, &mut seen)?;
		}
		Ok(header)
	}

	fn add_line(&mut self, parsed: Parsed, seen: &mut HashSet<String>) -> Result<(), Error> {
		let mut line = parsed.line;
		let identity = match &line.id {
			Some(id) => format!("{:?}:{id}", line.kind),
			None => line.text.clone(),
		};
		if !seen.insert(identity) {
			return Ok(());
		}
		match (line.kind, line.id.clone()) {
			(LineKind::FileFormat, _) => {
				self.version = parse_version(&line.text).unwrap_or(self.version);
			}
			(LineKind::Contig, Some(id)) => {
				let index =
					self.register_contig(&id, parsed.contig_length, parsed.explicit_index)?;
				line.dictionary_id = Some(index);
			}
			(LineKind::Info | LineKind::Format | LineKind::Filter, Some(id)) => {
				let index = self.register_string(&id, parsed.explicit_index)?;
				let entry = self
					.strings
					.get_mut(index as usize)
					.and_then(Option::as_mut)
					.ok_or(Error::Header("dictionary entry missing"))?;
				match line.kind {
					LineKind::Info => entry.info = parsed.definition,
					LineKind::Format => entry.format = parsed.definition,
					_ => entry.filter = true,
				}
				line.dictionary_id = Some(index);
			}
			_ => {}
		}
		self.lines.push(line);
		Ok(())
	}

	fn register_contig(
		&mut self,
		name: &str,
		length: Option<u64>,
		explicit: Option<u32>,
	) -> Result<u32, Error> {
		let index = u32::try_from(self.contigs.len()).map_err(|_| Error::Limits("contigs"))?;
		if explicit.is_some_and(|explicit| explicit != index) {
			return Err(Error::Header("contig IDX values are not sequential"));
		}
		self.contig_ids.insert(name.to_owned(), index);
		self.contigs.push(ContigEntry {
			name: name.to_owned(),
			length,
		});
		Ok(index)
	}

	fn register_string(&mut self, name: &str, explicit: Option<u32>) -> Result<u32, Error> {
		if let Some(&existing) = self.string_ids.get(name) {
			return Ok(existing);
		}
		let index = match explicit {
			Some(index) => index,
			None => u32::try_from(self.strings.len()).map_err(|_| Error::Limits("tags"))?,
		};
		let slot = index as usize;
		if slot >= self.strings.len() {
			let wanted = slot.checked_add(1).ok_or(Error::Limits("tags"))?;
			if wanted > MAX_DICTIONARY {
				return Err(Error::Limits("tags"));
			}
			self.strings.resize(wanted, None);
		}
		let entry = self
			.strings
			.get_mut(slot)
			.ok_or(Error::Header("dictionary entry missing"))?;
		if entry.is_some() {
			return Err(Error::Header("duplicate IDX"));
		}
		*entry = Some(StringEntry {
			name: name.to_owned(),
			info: None,
			format: None,
			filter: false,
		});
		self.string_ids.insert(name.to_owned(), index);
		Ok(index)
	}

	pub fn version(&self) -> Version {
		self.version
	}

	pub fn samples(&self) -> &[String] {
		&self.samples
	}

	pub fn sample_count(&self) -> usize {
		self.samples.len()
	}

	pub fn sample_index(&self, name: &str) -> Option<usize> {
		self.samples.iter().position(|sample| sample == name)
	}

	pub fn key(&self, name: &str) -> Option<u32> {
		self.string_ids.get(name).copied()
	}

	pub fn key_count(&self) -> usize {
		self.strings.len()
	}

	pub fn entry(&self, key: u32) -> Option<&StringEntry> {
		self.strings.get(key as usize)?.as_ref()
	}

	pub fn name(&self, key: u32) -> Option<&str> {
		self.entry(key).map(|entry| entry.name.as_str())
	}

	pub fn info_definition(&self, key: u32) -> Option<&FieldDefinition> {
		self.entry(key)?.info.as_ref()
	}

	pub fn format_definition(&self, key: u32) -> Option<&FieldDefinition> {
		self.entry(key)?.format.as_ref()
	}

	pub fn contig_id(&self, name: &str) -> Option<u32> {
		self.contig_ids.get(name).copied()
	}

	pub fn contig(&self, id: u32) -> Option<&ContigEntry> {
		self.contigs.get(id as usize)
	}

	pub fn contig_count(&self) -> usize {
		self.contigs.len()
	}

	pub fn contig_names(&self) -> impl Iterator<Item = &str> {
		self.contigs.iter().map(|contig| contig.name.as_str())
	}

	pub fn meta_lines(&self) -> impl Iterator<Item = &str> {
		self.lines.iter().map(|line| line.text.as_str())
	}

	pub fn contigs(&self) -> Result<Contigs, Error> {
		Ok(Contigs::new(self.contigs.iter().map(|contig| Contig {
			name: contig.name.clone(),
			length: contig.length.unwrap_or(UNKNOWN_CONTIG_LENGTH),
			checksum: None,
		}))?)
	}

	pub fn with_samples(&self, keep: &[usize]) -> Result<Self, Error> {
		let mut header = self.clone();
		header.samples = keep
			.iter()
			.map(|&index| self.samples.get(index).cloned())
			.collect::<Option<_>>()
			.ok_or(Error::Header("sample index out of range"))?;
		Ok(header)
	}

	pub fn without_formats(&self) -> Self {
		let mut header = self.clone();
		header
			.lines
			.retain(|line| !line.text.starts_with("##FORMAT="));
		header
	}

	pub fn with_lines(&self, extra: &[&str]) -> Result<Self, Error> {
		let mut meta: Vec<&str> = self.meta_lines().collect();
		meta.extend_from_slice(extra);
		Self::build(&meta, self.samples.clone())
	}

	fn with_lines_of(&self, other: &Self) -> Result<Self, Error> {
		for (_, entry) in other.keys() {
			let Some(existing) = self.key(&entry.name).and_then(|key| self.entry(key)) else {
				continue;
			};
			if conflicts(existing.info.as_ref(), entry.info.as_ref())
				|| conflicts(existing.format.as_ref(), entry.format.as_ref())
			{
				return Err(Error::Incompatible(
					"a tag has different Number or Type definitions",
				));
			}
		}
		let extra: Vec<&str> = other
			.meta_lines()
			.filter(|line| !line.starts_with("##fileformat="))
			.collect();
		self.with_lines(&extra)
	}

	pub fn merge(&self, other: &Self) -> Result<Self, Error> {
		if self.samples != other.samples {
			return Err(Error::Incompatible("sample names differ"));
		}
		self.with_lines_of(other)
	}

	pub fn join(&self, other: &Self) -> Result<Self, Error> {
		let mut header = self.with_lines_of(other)?;
		let known: HashSet<&str> = self.samples.iter().map(String::as_str).collect();
		if let Some(name) = other
			.samples
			.iter()
			.find(|name| known.contains(name.as_str()))
		{
			return Err(Error::DuplicateSample(name.clone()));
		}
		header.samples.extend_from_slice(&other.samples);
		Ok(header)
	}

	pub fn write_vcf(&self, out: &mut Vec<u8>) {
		for line in &self.lines {
			out.extend_from_slice(line.text.as_bytes());
			out.push(b'\n');
		}
		self.write_columns(out);
	}

	pub fn write_bcf(&self, out: &mut Vec<u8>) {
		for line in &self.lines {
			match (line.dictionary_id, line.text.strip_suffix('>')) {
				(Some(index), Some(body)) => {
					out.extend_from_slice(body.as_bytes());
					out.extend_from_slice(format!(",IDX={index}>").as_bytes());
				}
				_ => out.extend_from_slice(line.text.as_bytes()),
			}
			out.push(b'\n');
		}
		self.write_columns(out);
	}

	fn write_columns(&self, out: &mut Vec<u8>) {
		out.extend_from_slice(COLUMNS.as_bytes());
		if !self.samples.is_empty() {
			out.extend_from_slice(b"\tFORMAT");
			for sample in &self.samples {
				out.push(b'\t');
				out.extend_from_slice(sample.as_bytes());
			}
		}
		out.push(b'\n');
	}

	pub fn number_of(&self, key: u32) -> Option<Number> {
		self.info_definition(key)
			.or_else(|| self.format_definition(key))
			.map(|definition| definition.number)
	}

	pub fn type_of_info_by_name(&self, name: &str) -> Option<ValueType> {
		self.type_of_info(self.key(name)?)
	}

	pub fn type_of_info(&self, key: u32) -> Option<ValueType> {
		self.info_definition(key)
			.map(|definition| definition.value_type)
	}

	pub fn type_of_format(&self, key: u32) -> Option<ValueType> {
		self.format_definition(key)
			.map(|definition| definition.value_type)
	}

	pub fn keys(&self) -> impl Iterator<Item = (u32, &StringEntry)> {
		self.strings
			.iter()
			.zip(0u32..)
			.filter_map(|(entry, key)| Some((key, entry.as_ref()?)))
	}
}

const MAX_DICTIONARY: usize = 1 << 22;

fn conflicts(first: Option<&FieldDefinition>, second: Option<&FieldDefinition>) -> bool {
	first
		.zip(second)
		.is_some_and(|(first, second)| first != second)
}

fn parse_samples(line: &str) -> Result<Vec<String>, Error> {
	let mut columns = line.split('\t');
	let fixed = columns.by_ref().take(8).count();
	if fixed < 8 {
		return Err(Error::Header("#CHROM line has too few columns"));
	}
	match columns.next() {
		None => Ok(Vec::new()),
		Some("FORMAT") => Ok(columns.map(str::to_owned).collect()),
		Some(_) => Err(Error::Header("column nine is not FORMAT")),
	}
}

fn parse_version(line: &str) -> Option<Version> {
	let text = line.strip_prefix("##fileformat=VCFv")?;
	let (major, minor) = text.split_once('.')?;
	Some(Version {
		major: major.parse().ok()?,
		minor: minor.parse().ok()?,
	})
}
