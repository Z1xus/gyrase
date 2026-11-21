use std::io::BufRead;

use super::header::read_header_text;
use crate::Error;
use crate::header::Header;
use crate::vcf::VcfLine;

#[derive(Debug)]
pub struct VcfReader<R> {
	source: R,
	header: Header,
	line: Vec<u8>,
}

#[derive(Debug, Default)]
pub struct LineBatch {
	bytes: Vec<u8>,
	ends: Vec<usize>,
}

impl LineBatch {
	pub fn clear(&mut self) {
		self.bytes.clear();
		self.ends.clear();
	}

	pub fn len(&self) -> usize {
		self.ends.len()
	}

	pub fn is_empty(&self) -> bool {
		self.ends.is_empty()
	}

	pub fn byte_len(&self) -> usize {
		self.bytes.len()
	}

	pub fn lines(&self) -> impl Iterator<Item = Result<VcfLine<'_>, Error>> {
		let mut start = 0;
		self.ends.iter().map(move |&end| {
			let line = self.bytes.get(start..end).unwrap_or_default();
			start = end;
			VcfLine::new(line)
		})
	}
}

impl<R: BufRead> VcfReader<R> {
	pub fn new(mut source: R) -> Result<Self, Error> {
		let header = Header::parse(&read_header_text(&mut source)?)?;
		Ok(Self::with_header(source, header))
	}

	pub fn with_header(source: R, header: Header) -> Self {
		Self {
			source,
			header,
			line: Vec::new(),
		}
	}

	pub fn header(&self) -> &Header {
		&self.header
	}

	pub fn next_line(&mut self) -> Result<Option<VcfLine<'_>>, Error> {
		self.line.clear();
		if !read_record_line(&mut self.source, &mut self.line)? {
			return Ok(None);
		}
		VcfLine::new(&self.line).map(Some)
	}

	pub fn read_batch(&mut self, batch: &mut LineBatch, limit: usize) -> Result<bool, Error> {
		batch.clear();
		while batch.bytes.len() < limit {
			if !read_record_line(&mut self.source, &mut batch.bytes)? {
				break;
			}
			batch.ends.push(batch.bytes.len());
		}
		Ok(!batch.is_empty())
	}
}

fn read_record_line(source: &mut impl BufRead, out: &mut Vec<u8>) -> Result<bool, Error> {
	loop {
		let start = out.len();
		if source.read_until(b'\n', out)? == 0 {
			return Ok(false);
		}
		while matches!(out.last(), Some(b'\n' | b'\r')) && out.len() > start {
			out.pop();
		}
		if out.len() > start {
			return Ok(true);
		}
	}
}
