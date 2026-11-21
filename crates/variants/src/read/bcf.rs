use std::io::{self, BufRead, Read};

use crate::Error;
use crate::bcf::BcfRecord;
use crate::header::Header;

const MAJOR: u8 = 2;
const MAX_MINOR: u8 = 2;
const MAX_PART_LEN: u32 = 1 << 28;

#[derive(Debug)]
pub struct BcfReader<R> {
	source: R,
	header: Header,
	scratch: Vec<u8>,
}

#[derive(Debug, Default)]
pub struct RecordBatch {
	bytes: Vec<u8>,
	parts: Vec<(usize, usize)>,
}

impl RecordBatch {
	pub fn clear(&mut self) {
		self.bytes.clear();
		self.parts.clear();
	}

	pub fn len(&self) -> usize {
		self.parts.len()
	}

	pub fn is_empty(&self) -> bool {
		self.parts.is_empty()
	}

	pub fn byte_len(&self) -> usize {
		self.bytes.len()
	}

	pub fn records(&self) -> impl Iterator<Item = Result<BcfRecord<'_>, Error>> {
		let mut start = 0;
		self.parts.iter().map(move |&(shared_len, indiv_len)| {
			let shared_end = start + shared_len;
			let indiv_end = shared_end + indiv_len;
			let shared = self.bytes.get(start..shared_end).unwrap_or_default();
			let indiv = self.bytes.get(shared_end..indiv_end).unwrap_or_default();
			start = indiv_end;
			BcfRecord::new(shared, indiv)
		})
	}
}

impl<R: BufRead> BcfReader<R> {
	pub fn new(mut source: R) -> Result<Self, Error> {
		let mut magic = [0; 5];
		source.read_exact(&mut magic).map_err(|_| Error::NotBcf)?;
		let [b'B', b'C', b'F', major, minor] = magic else {
			return Err(Error::NotBcf);
		};
		if major != MAJOR || minor > MAX_MINOR {
			return Err(Error::BcfVersion { major, minor });
		}
		let text_len = read_length(&mut source)?.ok_or(Error::Header("truncated BCF header"))?;
		let mut text = Vec::new();
		(&mut source)
			.take(u64::from(text_len))
			.read_to_end(&mut text)?;
		if text.len() != text_len as usize {
			return Err(Error::Header("truncated BCF header"));
		}
		let header = Header::parse(&text)?;
		Ok(Self::with_header(source, header))
	}

	pub fn with_header(source: R, header: Header) -> Self {
		Self {
			source,
			header,
			scratch: Vec::new(),
		}
	}

	pub fn header(&self) -> &Header {
		&self.header
	}

	pub fn next_record(&mut self) -> Result<Option<BcfRecord<'_>>, Error> {
		self.scratch.clear();
		let Some((shared_len, _)) = read_record(&mut self.source, &mut self.scratch)? else {
			return Ok(None);
		};
		let (shared, indiv) = self.scratch.split_at(shared_len);
		BcfRecord::new(shared, indiv).map(Some)
	}

	pub fn read_batch(&mut self, batch: &mut RecordBatch, limit: usize) -> Result<bool, Error> {
		batch.clear();
		while batch.bytes.len() < limit {
			let Some(parts) = read_record(&mut self.source, &mut batch.bytes)? else {
				break;
			};
			batch.parts.push(parts);
		}
		Ok(!batch.is_empty())
	}
}

pub(crate) fn read_record(
	source: &mut impl BufRead,
	out: &mut Vec<u8>,
) -> Result<Option<(usize, usize)>, Error> {
	let Some(shared_len) = read_length(source)? else {
		return Ok(None);
	};
	let indiv_len = read_length(source)?.ok_or(Error::Bcf("record is truncated"))?;
	let shared_len = shared_len as usize;
	let indiv_len = indiv_len as usize;
	let start = out.len();
	out.resize(start + shared_len + indiv_len, 0);
	let body = out
		.get_mut(start..)
		.ok_or(Error::Bcf("record is truncated"))?;
	source.read_exact(body).map_err(|error| {
		if error.kind() == io::ErrorKind::UnexpectedEof {
			Error::Bcf("record is truncated")
		} else {
			Error::Io(error)
		}
	})?;
	Ok(Some((shared_len, indiv_len)))
}

fn read_length(source: &mut impl BufRead) -> Result<Option<u32>, Error> {
	if source.fill_buf()?.is_empty() {
		return Ok(None);
	}
	let mut word = [0; 4];
	source.read_exact(&mut word)?;
	let length = u32::from_le_bytes(word);
	if length > MAX_PART_LEN {
		return Err(Error::Limits("record part length"));
	}
	Ok(Some(length))
}
