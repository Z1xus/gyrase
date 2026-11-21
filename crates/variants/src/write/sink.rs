use std::io::{self, Write};

use gyrase_io::{BgzfWriter, BlockPosition, Finished};

use crate::Error;

pub trait Sink: Write {
	type Mark: Copy;
	type Done;

	fn mark(&self) -> Self::Mark;
	fn end_block(&mut self) -> Result<(), Error>;
	fn finish(self) -> Result<Self::Done, Error>;
}

impl<W: Write + Send + 'static> Sink for BgzfWriter<W> {
	type Mark = BlockPosition;
	type Done = Finished<W>;

	fn mark(&self) -> BlockPosition {
		self.position()
	}

	fn end_block(&mut self) -> Result<(), Error> {
		BgzfWriter::end_block(self).map_err(Error::from)
	}

	fn finish(self) -> Result<Finished<W>, Error> {
		BgzfWriter::finish(self).map_err(Error::from)
	}
}

#[derive(Debug)]
pub struct PlainSink<W> {
	inner: W,
	written: u64,
}

impl<W: Write> PlainSink<W> {
	pub fn new(inner: W) -> Self {
		Self { inner, written: 0 }
	}
}

impl<W: Write> Write for PlainSink<W> {
	fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
		let written = self.inner.write(bytes)?;
		self.written += written as u64;
		Ok(written)
	}

	fn flush(&mut self) -> io::Result<()> {
		self.inner.flush()
	}
}

impl<W: Write> Sink for PlainSink<W> {
	type Mark = u64;
	type Done = W;

	fn mark(&self) -> u64 {
		self.written
	}

	fn end_block(&mut self) -> Result<(), Error> {
		Ok(())
	}

	fn finish(mut self) -> Result<W, Error> {
		self.inner.flush()?;
		Ok(self.inner)
	}
}
