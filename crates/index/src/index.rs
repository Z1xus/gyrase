use std::io::{Read, Write};

use gyrase_io::{BgzfWriter, Level, Workers};

use crate::decode::decode;
use crate::encode::encode;
use crate::wire::read_all;
use crate::{Error, Geometry, RefIndex, TabixHeader};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
	Bai,
	Csi,
	Tbi,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Index {
	geometry: Geometry,
	references: Vec<RefIndex>,
	unplaced: u64,
	aux: Vec<u8>,
}

impl Index {
	pub(crate) fn assemble(
		geometry: Geometry,
		references: Vec<RefIndex>,
		unplaced: u64,
		aux: Vec<u8>,
	) -> Self {
		Self {
			geometry,
			references,
			unplaced,
			aux,
		}
	}

	pub fn with_tabix(mut self, header: &TabixHeader) -> Self {
		self.aux = header.to_bytes();
		self
	}

	pub fn geometry(&self) -> Geometry {
		self.geometry
	}

	pub fn references(&self) -> &[RefIndex] {
		&self.references
	}

	pub fn unplaced(&self) -> u64 {
		self.unplaced
	}

	pub fn aux(&self) -> &[u8] {
		&self.aux
	}

	pub fn tabix(&self) -> Result<TabixHeader, Error> {
		TabixHeader::parse(&self.aux)
	}

	pub fn read(input: impl Read + Send + 'static, workers: &Workers) -> Result<Self, Error> {
		decode(&read_all(gyrase_io::decode(Box::new(input), workers)?)?)
	}

	pub fn write<W: Write + Send + 'static>(
		&self,
		format: Format,
		mut out: W,
		workers: &Workers,
	) -> Result<W, Error> {
		let bytes = encode(self, format)?;
		match format {
			Format::Bai => {
				out.write_all(&bytes)?;
				Ok(out)
			}
			Format::Csi | Format::Tbi => {
				let mut writer = BgzfWriter::new(out, Level::DEFAULT, workers);
				writer.write_all(&bytes)?;
				Ok(writer.finish()?.inner)
			}
		}
	}
}
