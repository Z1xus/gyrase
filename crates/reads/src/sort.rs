mod buffer;
mod merge;
mod run;

use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};

use gyrase_index::{Format, Index};
use gyrase_io::{Level, Workers};

use self::buffer::Buffer;
pub use self::buffer::Encoded;
use self::run::RunFile;
use crate::bam::{IndexedWriter, index_layout};
use crate::error::{Error, RecordError, SortError};
use crate::header::{Header, SortOrder};
use crate::markdup::{MarkSettings, MarkStats, Marker, is_duplicate, set_duplicate};
use crate::record::RecordRef;

const CHECK_EVERY: usize = 1 << 16;
const POSITION_LIMIT: u64 = 1 << 31;

#[derive(Debug)]
pub struct Sorter<'w> {
	header: Header,
	workers: &'w Workers,
	cancel: &'w AtomicBool,
	directory: PathBuf,
	buffer: Buffer,
	runs: Vec<RunFile>,
}

#[derive(Debug)]
pub struct Records {
	buffer: Buffer,
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
		cancel: &'w AtomicBool,
	) -> Result<Self, Error> {
		let mut header = header.clone();
		header.set_sort_order(SortOrder::Coordinate);
		Ok(Self {
			header,
			workers,
			cancel,
			directory: directory.into(),
			buffer: Buffer::new(budget)?,
			runs: Vec::new(),
		})
	}

	pub fn push(&mut self, record: RecordRef<'_>) -> Result<(), Error> {
		let (core, rest) = record.parts();
		let len = u32::try_from(core.len() + rest.len()).map_err(|_| SortError::RecordTooLarge)?;
		self.reserve(len as usize, 1)?;
		self.buffer.push(record, len);
		Ok(())
	}

	pub fn push_encoded(&mut self, chunk: &Encoded) -> Result<(), Error> {
		self.reserve(chunk.bytes(), chunk.records())?;
		self.buffer.append(chunk);
		Ok(())
	}

	pub fn spilled(&self) -> bool {
		!self.runs.is_empty()
	}

	pub fn finish_records(
		mut self,
		settings: MarkSettings,
	) -> Result<Result<(Records, MarkStats), Self>, Error> {
		if self.spilled() {
			return Ok(Err(self));
		}
		self.buffer.sort(self.workers, self.cancel)?;
		let mut marker = Marker::new(settings, &self.header);
		let mut changed = Vec::new();
		let mut next = 0;
		let buffer = &self.buffer;
		let mut emit = |bytes: &[u8]| {
			if buffer.record_at(next).map(is_duplicate) != Some(is_duplicate(bytes)) {
				changed.push((next, is_duplicate(bytes)));
			}
			next += 1;
			Ok::<(), Error>(())
		};
		for index in 0..buffer.len() {
			if index.is_multiple_of(CHECK_EVERY) {
				check(self.cancel)?;
			}
			marker.push(buffer.record_at(index).unwrap_or_default(), &mut emit)?;
		}
		let stats = marker.finish(&mut emit)?;
		for (index, duplicate) in changed {
			if let Some(bytes) = self.buffer.record_mut(index) {
				set_duplicate(bytes, duplicate);
			}
		}
		Ok(Ok((
			Records {
				buffer: self.buffer,
			},
			stats,
		)))
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
		let (sorted, stats) = self.write(output, level, Some(settings))?;
		Ok((sorted, stats.unwrap_or_default()))
	}

	fn write<W: Write + Send + 'static>(
		mut self,
		output: W,
		level: Level,
		settings: Option<MarkSettings>,
	) -> Result<(Sorted<W>, Option<MarkStats>), Error> {
		self.buffer.sort(self.workers, self.cancel)?;
		merge::collapse(&mut self.runs, &self.directory, self.workers, self.cancel)?;
		let (format, geometry) = index_layout(&self.header)?;
		let mut marker = settings.map(|settings| Marker::new(settings, &self.header));
		let mut writer = IndexedWriter::new(output, &self.header, geometry, level, self.workers)?;
		let mut emit = |bytes: &[u8]| writer.write(RecordRef::new(bytes)?);
		merge::merge(
			Some(&self.buffer),
			&self.runs,
			self.workers,
			self.cancel,
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

	fn reserve(&mut self, bytes: usize, records: usize) -> Result<(), Error> {
		if !self.buffer.fits(bytes, records) {
			self.spill()?;
			if !self.buffer.fits(bytes, records) {
				return Err(SortError::RecordTooLarge.into());
			}
		}
		Ok(())
	}

	fn spill(&mut self) -> Result<(), Error> {
		if self.buffer.len() == 0 {
			return Ok(());
		}
		self.buffer.sort(self.workers, self.cancel)?;
		self.runs.push(RunFile::write(
			&self.directory,
			&self.buffer,
			self.workers,
			self.cancel,
		)?);
		self.buffer.clear();
		Ok(())
	}
}

impl Records {
	pub fn len(&self) -> usize {
		self.buffer.len()
	}

	pub fn is_empty(&self) -> bool {
		self.buffer.len() == 0
	}

	pub fn count(&self, contig: u32, start: u64, end: u64) -> usize {
		let (first, last) = self.candidates(contig, start, end);
		last - first
	}

	pub fn overlapping(
		&self,
		contig: u32,
		start: u64,
		end: u64,
	) -> impl Iterator<Item = Result<RecordRef<'_>, RecordError>> {
		let (first, last) = self.candidates(contig, start, end);
		(first..last)
			.filter_map(|index| self.buffer.record_at(index))
			.map(RecordRef::new)
			.filter(move |record| {
				record.as_ref().map_or(true, |record| {
					record.reference_end().is_some_and(|stop| stop > start)
				})
			})
	}

	fn candidates(&self, contig: u32, start: u64, end: u64) -> (usize, usize) {
		let key = |position: u64| u64::from(contig) << 32 | position.min(POSITION_LIMIT) << 1;
		(
			self.buffer
				.lower_bound(key(start.saturating_sub(self.buffer.longest()))),
			self.buffer.lower_bound(key(end)),
		)
	}
}

fn check(cancel: &AtomicBool) -> Result<(), Error> {
	if cancel.load(Ordering::Relaxed) {
		Err(Error::Cancelled)
	} else {
		Ok(())
	}
}
