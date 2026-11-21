use std::io::BufRead;

use super::parse::Parser;
use crate::error::Error;
use crate::header::Header;
use crate::record::Record;

#[derive(Debug)]
pub struct SamReader<R: BufRead> {
	source: R,
	header: Header,
	line: Vec<u8>,
	number: u64,
	parser: Parser,
}

impl<R: BufRead> SamReader<R> {
	pub fn new(mut source: R) -> Result<Self, Error> {
		let mut text = Vec::new();
		let mut number = 0;
		while source.fill_buf()?.first() == Some(&b'@') {
			source.read_until(b'\n', &mut text)?;
			number += 1;
		}
		Ok(Self {
			source,
			header: Header::parse(&text)?,
			line: Vec::new(),
			number,
			parser: Parser::default(),
		})
	}

	pub fn header(&self) -> &Header {
		&self.header
	}

	pub fn read_record(&mut self, record: &mut Record) -> Result<bool, Error> {
		loop {
			self.line.clear();
			if self.source.read_until(b'\n', &mut self.line)? == 0 {
				return Ok(false);
			}
			self.number += 1;
			while matches!(self.line.last(), Some(b'\n' | b'\r')) {
				self.line.pop();
			}
			if self.line.is_empty() {
				continue;
			}
			self.parser
				.parse(&self.line, self.number, &self.header, record)?;
			return Ok(true);
		}
	}
}
