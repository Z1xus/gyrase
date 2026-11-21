use crate::Error;
use crate::header::Header;
use crate::record::Record;
use crate::value::Loaded;

#[derive(Debug, Clone, Default)]
pub struct Alleles {
	bytes: Vec<u8>,
	ends: Vec<usize>,
}

impl Alleles {
	pub fn clear(&mut self) {
		self.bytes.clear();
		self.ends.clear();
	}

	pub fn push(&mut self, allele: &[u8]) {
		self.bytes.extend_from_slice(allele);
		self.ends.push(self.bytes.len());
	}

	pub fn len(&self) -> usize {
		self.ends.len()
	}

	pub fn is_empty(&self) -> bool {
		self.ends.is_empty()
	}

	pub fn get(&self, index: usize) -> Option<&[u8]> {
		let end = *self.ends.get(index)?;
		let start = match index.checked_sub(1) {
			Some(previous) => *self.ends.get(previous)?,
			None => 0,
		};
		self.bytes.get(start..end)
	}

	pub fn iter(&self) -> impl Iterator<Item = &[u8]> {
		(0..self.len()).filter_map(|index| self.get(index))
	}
}

pub trait Site {
	fn contig(&self, header: &Header) -> Result<u32, Error>;
	fn start(&self) -> Result<u32, Error>;
	fn end(&self, header: &Header) -> Result<u32, Error>;
	fn id(&self) -> Result<&[u8], Error>;
	fn alleles(&self, out: &mut Alleles) -> Result<(), Error>;
	fn qual(&self) -> Result<Option<f32>, Error>;
	fn filters(&self, header: &Header, out: &mut Vec<u32>) -> Result<(), Error>;
	fn info(&self, header: &Header, key: u32, out: &mut Loaded) -> Result<bool, Error>;
	fn format_keys(&self, header: &Header, out: &mut Vec<u32>) -> Result<(), Error>;
	fn format(&self, header: &Header, key: u32, out: &mut Loaded) -> Result<bool, Error>;
	fn read_record(&self, header: &Header, out: &mut Record) -> Result<(), Error>;
}
