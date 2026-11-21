use std::io::Write;

use super::format::format_record;
use crate::error::Error;
use crate::header::Header;
use crate::record::RecordRef;

const FLUSH_BYTES: usize = 1 << 16;

#[derive(Debug)]
pub struct SamWriter<W: Write> {
	inner: W,
	header: Header,
	line: Vec<u8>,
}

impl<W: Write> SamWriter<W> {
	pub fn new(mut inner: W, header: &Header) -> Result<Self, Error> {
		inner.write_all(header.text().as_bytes())?;
		Ok(Self {
			inner,
			header: header.clone(),
			line: Vec::with_capacity(2 * FLUSH_BYTES),
		})
	}

	pub fn write(&mut self, record: RecordRef<'_>) -> Result<(), Error> {
		format_record(&mut self.line, &self.header, record)?;
		self.line.push(b'\n');
		if self.line.len() >= FLUSH_BYTES {
			self.drain()?;
		}
		Ok(())
	}

	pub fn finish(mut self) -> Result<W, Error> {
		self.drain()?;
		self.inner.flush()?;
		Ok(self.inner)
	}

	fn drain(&mut self) -> Result<(), Error> {
		self.inner.write_all(&self.line)?;
		self.line.clear();
		Ok(())
	}
}
