use std::io::{BufRead, Read};

use gyrase_io::{Decoded, Workers, decode};

use super::bcf::BcfReader;
use super::vcf::VcfReader;
use crate::Error;
use crate::header::Header;

const BCF_MAGIC: &[u8] = b"BCF";

#[derive(Debug)]
pub enum Reader<R = Decoded> {
	Vcf(VcfReader<R>),
	Bcf(BcfReader<R>),
}

impl Reader<Decoded> {
	pub fn open(source: Box<dyn Read + Send>, workers: &Workers) -> Result<Self, Error> {
		let mut decoded = decode(source, workers)?;
		if decoded.fill_buf()?.starts_with(BCF_MAGIC) {
			BcfReader::new(decoded).map(Self::Bcf)
		} else {
			VcfReader::new(decoded).map(Self::Vcf)
		}
	}
}

impl<R: BufRead> Reader<R> {
	pub fn header(&self) -> &Header {
		match self {
			Self::Vcf(reader) => reader.header(),
			Self::Bcf(reader) => reader.header(),
		}
	}
}
