use std::io::Write;

use gyrase_io::{BgzfWriter, BlockPosition, Finished, Level, VirtualOffset, Workers};

use crate::error::{BamError, Error};
use crate::header::Header;
use crate::record::RecordRef;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Written {
	pub start: BlockPosition,
	pub end: BlockPosition,
}

#[derive(Debug)]
pub struct BamWriter<W: Write + Send + 'static> {
	inner: BgzfWriter<W>,
}

impl<W: Write + Send + 'static> BamWriter<W> {
	pub fn new(inner: W, header: &Header, level: Level, workers: &Workers) -> Result<Self, Error> {
		let mut encoded = Vec::new();
		header.write_bam(&mut encoded)?;
		let mut inner = BgzfWriter::new(inner, level, workers);
		inner.write_all(&encoded)?;
		inner.end_block()?;
		Ok(Self { inner })
	}

	pub fn position(&self) -> BlockPosition {
		self.inner.position()
	}

	pub fn resolve(&mut self, position: BlockPosition) -> Option<VirtualOffset> {
		self.inner.resolve(position)
	}

	pub fn write(&mut self, record: RecordRef<'_>) -> Result<Written, Error> {
		let (core, rest) = record.parts();
		let len = u32::try_from(core.len() + rest.len())
			.map_err(|_| BamError::RecordSize { len: u64::MAX })?;
		let start = self.inner.position();
		self.inner.write_all(&len.to_le_bytes())?;
		self.inner.write_all(core)?;
		self.inner.write_all(rest)?;
		Ok(Written {
			start,
			end: self.inner.position(),
		})
	}

	pub fn finish(self) -> Result<Finished<W>, Error> {
		Ok(self.inner.finish()?)
	}
}
