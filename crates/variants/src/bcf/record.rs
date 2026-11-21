use super::typed::{Cursor, Typed};
use crate::Error;
use crate::header::Header;
use crate::record::Record;
use crate::site::{Alleles, Site};
use crate::value::{Cells, Loaded, MISSING_FLOAT_BITS};

const FIXED_LEN: usize = 24;

#[derive(Debug, Clone, Copy)]
pub struct BcfRecord<'a> {
	shared: &'a [u8],
	indiv: &'a [u8],
	allele_count: usize,
	info_count: usize,
	format_count: usize,
	sample_count: usize,
	id_at: usize,
	alleles_at: usize,
	filter_at: usize,
	info_at: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct InfoEntry<'a> {
	pub key: u32,
	pub value: Typed<'a>,
}

#[derive(Debug, Clone, Copy)]
pub struct FormatEntry<'a> {
	pub key: u32,
	pub values: Typed<'a>,
}

impl<'a> BcfRecord<'a> {
	pub fn new(shared: &'a [u8], indiv: &'a [u8]) -> Result<Self, Error> {
		if shared.len() < FIXED_LEN {
			return Err(Error::Bcf("shared part is shorter than the fixed fields"));
		}
		let allele_info = word(shared, 16)?;
		let format_samples = word(shared, 20)?;
		let allele_count = usize::from(u16::try_from(allele_info >> 16).unwrap_or(0));
		let mut record = Self {
			shared,
			indiv,
			allele_count,
			info_count: usize::from(u16::try_from(allele_info & 0xffff).unwrap_or(0)),
			format_count: usize::try_from(format_samples >> 24).unwrap_or(0),
			sample_count: usize::try_from(format_samples & 0x00ff_ffff).unwrap_or(0),
			id_at: FIXED_LEN,
			alleles_at: 0,
			filter_at: 0,
			info_at: 0,
		};
		let mut cursor = Cursor::at(shared, FIXED_LEN);
		cursor.typed()?;
		record.alleles_at = cursor.position();
		for _ in 0..allele_count {
			cursor.typed()?;
		}
		record.filter_at = cursor.position();
		cursor.typed()?;
		record.info_at = cursor.position();
		Ok(record)
	}

	pub fn contig(&self) -> Result<u32, Error> {
		let raw = i32::from_le_bytes(array(self.shared, 0)?);
		u32::try_from(raw).map_err(|_| Error::Bcf("negative contig id"))
	}

	pub fn start(&self) -> Result<u32, Error> {
		let raw = i32::from_le_bytes(array(self.shared, 4)?);
		u32::try_from(raw).map_err(|_| Error::Bcf("negative position"))
	}

	pub fn length(&self) -> Result<u32, Error> {
		let raw = i32::from_le_bytes(array(self.shared, 8)?);
		u32::try_from(raw).map_err(|_| Error::Bcf("negative reference length"))
	}

	pub fn qual(&self) -> Result<Option<f32>, Error> {
		let bits = word(self.shared, 12)?;
		Ok((bits != MISSING_FLOAT_BITS).then(|| f32::from_bits(bits)))
	}

	pub fn allele_count(&self) -> usize {
		self.allele_count
	}

	pub fn info_count(&self) -> usize {
		self.info_count
	}

	pub fn format_count(&self) -> usize {
		self.format_count
	}

	pub fn sample_count(&self) -> usize {
		self.sample_count
	}

	pub fn shared(&self) -> &'a [u8] {
		self.shared
	}

	pub fn indiv(&self) -> &'a [u8] {
		self.indiv
	}

	pub fn id(&self) -> Result<&'a [u8], Error> {
		let typed = Cursor::at(self.shared, self.id_at).typed()?;
		Ok(typed.bytes)
	}

	pub fn alleles(&self) -> Result<impl Iterator<Item = Result<&'a [u8], Error>>, Error> {
		let mut cursor = Cursor::at(self.shared, self.alleles_at);
		let count = self.allele_count;
		Ok((0..count).map(move |_| cursor.typed().map(|typed| typed.bytes)))
	}

	pub fn filters(&self) -> Result<Typed<'a>, Error> {
		Cursor::at(self.shared, self.filter_at).typed()
	}

	pub fn info(&self) -> impl Iterator<Item = Result<InfoEntry<'a>, Error>> {
		let mut cursor = Cursor::at(self.shared, self.info_at);
		(0..self.info_count).map(move |_| {
			let key = cursor
				.typed()?
				.scalar_int()
				.and_then(|key| u32::try_from(key).ok())
				.ok_or(Error::Bcf("invalid INFO key"))?;
			let value = cursor.typed()?;
			Ok(InfoEntry { key, value })
		})
	}

	pub fn format(&self) -> impl Iterator<Item = Result<FormatEntry<'a>, Error>> {
		let mut cursor = Cursor::new(self.indiv);
		let samples = self.sample_count;
		(0..self.format_count).map(move |_| {
			let key = cursor
				.typed()?
				.scalar_int()
				.and_then(|key| u32::try_from(key).ok())
				.ok_or(Error::Bcf("invalid FORMAT key"))?;
			let values = cursor.typed_repeated(samples)?;
			Ok(FormatEntry { key, values })
		})
	}

	pub fn load_info(&self, key: u32, out: &mut Loaded) -> Result<bool, Error> {
		out.clear();
		for entry in self.info() {
			let entry = entry?;
			if entry.key == key {
				out.cells = entry.value.push_to(&mut out.store);
				out.width = entry.value.len;
				return Ok(true);
			}
		}
		Ok(false)
	}

	pub fn load_format(&self, key: u32, out: &mut Loaded) -> Result<bool, Error> {
		out.clear();
		for entry in self.format() {
			let entry = entry?;
			if entry.key == key {
				out.cells = entry.values.push_to(&mut out.store);
				out.width = entry.values.len;
				if matches!(out.cells, Cells::Flag) {
					return Err(Error::Bcf("FORMAT field without a type"));
				}
				return Ok(true);
			}
		}
		Ok(false)
	}
}

fn array(bytes: &[u8], at: usize) -> Result<[u8; 4], Error> {
	bytes
		.get(at..at + 4)
		.and_then(|slice| slice.try_into().ok())
		.ok_or(Error::Bcf("record is truncated"))
}

fn word(bytes: &[u8], at: usize) -> Result<u32, Error> {
	Ok(u32::from_le_bytes(array(bytes, at)?))
}

impl Site for BcfRecord<'_> {
	fn contig(&self, _header: &Header) -> Result<u32, Error> {
		BcfRecord::contig(self)
	}

	fn start(&self) -> Result<u32, Error> {
		BcfRecord::start(self)
	}

	fn end(&self, _header: &Header) -> Result<u32, Error> {
		BcfRecord::start(self)?
			.checked_add(self.length()?)
			.ok_or(Error::Bcf("record end overflows"))
	}

	fn id(&self) -> Result<&[u8], Error> {
		BcfRecord::id(self)
	}

	fn alleles(&self, out: &mut Alleles) -> Result<(), Error> {
		out.clear();
		for allele in BcfRecord::alleles(self)? {
			out.push(allele?);
		}
		Ok(())
	}

	fn qual(&self) -> Result<Option<f32>, Error> {
		BcfRecord::qual(self)
	}

	fn filters(&self, _header: &Header, out: &mut Vec<u32>) -> Result<(), Error> {
		out.clear();
		for key in self.filters()?.ints() {
			out.push(u32::try_from(key).map_err(|_| Error::Bcf("invalid FILTER key"))?);
		}
		Ok(())
	}

	fn info(&self, _header: &Header, key: u32, out: &mut Loaded) -> Result<bool, Error> {
		self.load_info(key, out)
	}

	fn format_keys(&self, _header: &Header, out: &mut Vec<u32>) -> Result<(), Error> {
		out.clear();
		for entry in BcfRecord::format(self) {
			out.push(entry?.key);
		}
		Ok(())
	}

	fn format(&self, _header: &Header, key: u32, out: &mut Loaded) -> Result<bool, Error> {
		self.load_format(key, out)
	}

	fn read_record(&self, _header: &Header, out: &mut Record) -> Result<(), Error> {
		out.read_bcf(self)
	}
}
