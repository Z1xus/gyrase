use std::fs::{self, File};
use std::io::{self, BufRead, Read, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};

use gyrase_io::{BgzfReader, BgzfWriter, Input, Level, Readahead, Workers};

use super::buffer::{Buffer, key};
use crate::error::{Error, SortError};

static SERIAL: AtomicU64 = AtomicU64::new(0);

#[derive(Debug)]
pub(super) struct RunFile {
	path: PathBuf,
}

impl Drop for RunFile {
	fn drop(&mut self) {
		let _ = fs::remove_file(&self.path);
	}
}

impl RunFile {
	pub(super) fn write(
		directory: &Path,
		buffer: &Buffer,
		workers: &Workers,
	) -> Result<Self, Error> {
		let path = directory.join(format!(
			"gyrase-sort-{}-{}.bgzf",
			process::id(),
			SERIAL.fetch_add(1, Ordering::Relaxed)
		));
		let file = File::options().write(true).create_new(true).open(&path)?;
		let run = Self { path };
		let mut writer = BgzfWriter::new(file, Level::FASTEST, workers);
		for index in 0..buffer.len() {
			let record = buffer.record_at(index).unwrap_or_default();
			let len = u32::try_from(record.len()).unwrap_or(u32::MAX);
			writer.write_all(&len.to_le_bytes())?;
			writer.write_all(record)?;
		}
		writer.finish()?;
		Ok(run)
	}

	pub(super) fn open(&self, workers: &Workers) -> Result<RunReader, Error> {
		let file = File::open(&self.path)?;
		Ok(RunReader {
			reader: BgzfReader::with_readahead(
				Input::Stream(Box::new(file)),
				workers,
				Readahead::COMPACT,
			),
			record: Vec::new(),
			key: 0,
		})
	}
}

#[derive(Debug)]
pub(super) struct RunReader {
	reader: BgzfReader,
	record: Vec<u8>,
	key: u64,
}

impl RunReader {
	pub(super) fn key(&self) -> u64 {
		self.key
	}

	pub(super) fn record(&self) -> &[u8] {
		&self.record
	}

	pub(super) fn advance(&mut self) -> Result<bool, Error> {
		if self.reader.fill_buf()?.is_empty() {
			return Ok(false);
		}
		let mut word = [0; 4];
		self.read(&mut word)?;
		self.record.resize(u32::from_le_bytes(word) as usize, 0);
		let mut record = std::mem::take(&mut self.record);
		let read = self.read(&mut record);
		self.record = record;
		read?;
		self.key = key(&self.record);
		Ok(true)
	}

	fn read(&mut self, bytes: &mut [u8]) -> Result<(), Error> {
		self.reader.read_exact(bytes).map_err(|error| {
			if error.kind() == io::ErrorKind::UnexpectedEof {
				SortError::TruncatedRun.into()
			} else {
				Error::from(error)
			}
		})
	}
}
