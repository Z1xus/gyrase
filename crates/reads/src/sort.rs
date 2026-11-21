mod buffer;
mod merge;
mod run;

use std::io::Write;
use std::path::PathBuf;

use gyrase_index::{Format, Index};
use gyrase_io::{Level, Workers};

use self::buffer::Buffer;
use self::run::RunFile;
use crate::bam::{IndexedWriter, index_layout};
use crate::error::{Error, SortError};
use crate::header::{Header, SortOrder};
use crate::markdup::{MarkSettings, MarkStats, Marker};
use crate::record::RecordRef;

#[derive(Debug)]
pub struct Sorter<'w> {
	header: Header,
	workers: &'w Workers,
	directory: PathBuf,
	buffer: Buffer,
	runs: Vec<RunFile>,
}

#[derive(Debug)]
pub struct Sorted<W> {
	pub output: W,
	pub index: Index,
	pub format: Format,
}

impl<'w> Sorter<'w> {
	pub fn new(
		header: &Header,
		budget: usize,
		directory: impl Into<PathBuf>,
		workers: &'w Workers,
	) -> Result<Self, Error> {
		let mut header = header.clone();
		header.set_sort_order(SortOrder::Coordinate);
		Ok(Self {
			header,
			workers,
			directory: directory.into(),
			buffer: Buffer::new(budget)?,
			runs: Vec::new(),
		})
	}

	pub fn push(&mut self, record: RecordRef<'_>) -> Result<(), Error> {
		let (core, rest) = record.parts();
		let len = u32::try_from(core.len() + rest.len()).map_err(|_| SortError::RecordTooLarge)?;
		if !self.buffer.fits(len as usize) {
			self.spill()?;
			if !self.buffer.fits(len as usize) {
				return Err(SortError::RecordTooLarge.into());
			}
		}
		self.buffer.push(core, rest, len);
		Ok(())
	}

	pub fn finish<W: Write + Send + 'static>(
		self,
		output: W,
		level: Level,
	) -> Result<Sorted<W>, Error> {
		Ok(self.write(output, level, None)?.0)
	}

	pub fn finish_marking<W: Write + Send + 'static>(
		self,
		output: W,
		level: Level,
		settings: MarkSettings,
	) -> Result<(Sorted<W>, MarkStats), Error> {
		let (sorted, stats) = self.write(output, level, Some(Marker::new(settings)))?;
		Ok((sorted, stats.unwrap_or_default()))
	}

	fn write<W: Write + Send + 'static>(
		mut self,
		output: W,
		level: Level,
		mut marker: Option<Marker>,
	) -> Result<(Sorted<W>, Option<MarkStats>), Error> {
		self.buffer.sort(self.workers);
		let (format, geometry) = index_layout(&self.header)?;
		let mut writer = IndexedWriter::new(output, &self.header, geometry, level, self.workers)?;
		let mut emit = |bytes: &[u8]| writer.write(RecordRef::new(bytes)?);
		merge::merge(
			&self.buffer,
			&self.runs,
			self.workers,
			|bytes| match marker.as_mut() {
				Some(marker) => marker.push(bytes, &mut emit),
				None => emit(bytes),
			},
		)?;
		let stats = marker.map(|marker| marker.finish(&mut emit)).transpose()?;
		let (output, index) = writer.finish()?;
		Ok((
			Sorted {
				output,
				index,
				format,
			},
			stats,
		))
	}

	fn spill(&mut self) -> Result<(), Error> {
		if self.buffer.len() == 0 {
			return Ok(());
		}
		self.buffer.sort(self.workers);
		self.runs
			.push(RunFile::write(&self.directory, &self.buffer, self.workers)?);
		self.buffer.clear();
		Ok(())
	}
}
