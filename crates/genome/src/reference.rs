use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use md5::{Digest, Md5};

use crate::{Contig, Contigs, Error, Interval};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaiEntry {
	pub name: String,
	pub length: u64,
	pub offset: u64,
	pub line_bases: u64,
	pub line_width: u64,
}

#[derive(Debug)]
pub struct Reference {
	file: File,
	file_len: u64,
	entries: Vec<FaiEntry>,
	contigs: Contigs,
}

const CHECKSUM_CHUNK: u64 = 1 << 20;

impl Reference {
	pub fn open(path: &Path) -> Result<Self, Error> {
		Self::load(path, false)
	}

	pub fn open_indexed(path: &Path) -> Result<Self, Error> {
		Self::load(path, true)
	}

	fn load(path: &Path, require_index: bool) -> Result<Self, Error> {
		let file = File::open(path)?;
		let index_path = index_path(path);
		let entries = if let Ok(index) = File::open(&index_path) {
			read_index(BufReader::new(index))?
		} else {
			let entries = build_index(BufReader::with_capacity(1 << 20, File::open(path)?))?;
			let written = write_index(&index_path, &entries);
			if require_index {
				written?;
			}
			entries
		};
		let contigs = Contigs::new(entries.iter().map(|entry| Contig {
			name: entry.name.clone(),
			length: entry.length,
			checksum: None,
		}))?;
		let file_len = file.metadata()?.len();
		Ok(Self {
			file,
			file_len,
			entries,
			contigs,
		})
	}

	pub fn contigs(&self) -> &Contigs {
		&self.contigs
	}

	pub fn entries(&self) -> &[FaiEntry] {
		&self.entries
	}

	pub fn fetch(&self, interval: &Interval, bases: &mut Vec<u8>) -> Result<(), Error> {
		self.fetch_cased(interval, bases)?;
		bases.make_ascii_uppercase();
		Ok(())
	}

	pub fn fetch_cased(&self, interval: &Interval, bases: &mut Vec<u8>) -> Result<(), Error> {
		bases.clear();
		let entry = self
			.entries
			.get(interval.contig.index())
			.ok_or_else(|| Error::UnknownContig(interval.contig.get().to_string()))?;
		if interval.end > entry.length || interval.start > interval.end {
			return Err(Error::OutOfBounds {
				contig: entry.name.clone(),
				start: interval.start,
				end: interval.end,
			});
		}
		if interval.is_empty() {
			return Ok(());
		}
		let invalid = |offset| Error::InvalidFasta { offset };
		let first = entry
			.byte_offset(interval.start)
			.ok_or_else(|| invalid(entry.offset))?;
		let last = entry
			.byte_offset(interval.end - 1)
			.and_then(|last| last.checked_add(1))
			.filter(|&last| last >= first && last <= self.file_len)
			.ok_or_else(|| invalid(first))?;
		let len = usize::try_from(last - first).map_err(|_| invalid(first))?;
		bases.resize(len, 0);
		read_exact_at(&self.file, bases, first)?;
		bases.retain(|byte| !byte.is_ascii_whitespace());
		if bases.len() as u64 != interval.len() {
			return Err(Error::InvalidFasta { offset: first });
		}
		Ok(())
	}

	pub fn checksum(&self, interval: &Interval) -> Result<[u8; 16], Error> {
		let mut hasher = Md5::new();
		let mut bases = Vec::new();
		let mut start = interval.start;
		while start < interval.end {
			let end = (start + CHECKSUM_CHUNK).min(interval.end);
			self.fetch(
				&Interval {
					contig: interval.contig,
					start,
					end,
				},
				&mut bases,
			)?;
			hasher.update(&bases);
			start = end;
		}
		Ok(hasher.finalize().into())
	}
}

impl FaiEntry {
	fn byte_offset(&self, position: u64) -> Option<u64> {
		(position / self.line_bases)
			.checked_mul(self.line_width)?
			.checked_add(position % self.line_bases)?
			.checked_add(self.offset)
	}
}

fn index_path(path: &Path) -> PathBuf {
	let mut name = path.as_os_str().to_owned();
	name.push(".fai");
	PathBuf::from(name)
}

fn read_index(reader: impl BufRead) -> Result<Vec<FaiEntry>, Error> {
	reader
		.lines()
		.enumerate()
		.filter(|(_, line)| line.as_ref().map_or(true, |line| !line.is_empty()))
		.map(|(number, line)| {
			let line = line?;
			let invalid = || Error::InvalidIndex { line: number + 1 };
			let mut fields = line.split('\t');
			let name = fields.next().ok_or_else(invalid)?.to_owned();
			let mut number_field = || {
				fields
					.next()
					.and_then(|field| field.parse::<u64>().ok())
					.ok_or_else(invalid)
			};
			let entry = FaiEntry {
				name,
				length: number_field()?,
				offset: number_field()?,
				line_bases: number_field()?,
				line_width: number_field()?,
			};
			if entry.line_bases == 0 || entry.line_width < entry.line_bases {
				return Err(invalid());
			}
			Ok(entry)
		})
		.collect()
}

fn write_index(path: &Path, entries: &[FaiEntry]) -> Result<(), Error> {
	let mut temporary = path.as_os_str().to_owned();
	temporary.push(format!(".{}.tmp", std::process::id()));
	let temporary = PathBuf::from(temporary);
	let written = write_entries(&temporary, entries).and_then(|()| fs::rename(&temporary, path));
	if written.is_err() {
		let _ = fs::remove_file(&temporary);
	}
	Ok(written?)
}

fn write_entries(path: &Path, entries: &[FaiEntry]) -> std::io::Result<()> {
	let mut writer = BufWriter::new(File::create(path)?);
	for entry in entries {
		writeln!(
			writer,
			"{}\t{}\t{}\t{}\t{}",
			entry.name, entry.length, entry.offset, entry.line_bases, entry.line_width
		)?;
	}
	writer.flush()
}

struct Scan {
	entry: FaiEntry,
	short_line_seen: bool,
}

fn build_index(mut reader: impl BufRead) -> Result<Vec<FaiEntry>, Error> {
	let mut entries = Vec::new();
	let mut current: Option<Scan> = None;
	let mut offset = 0u64;
	let mut line = Vec::new();
	loop {
		line.clear();
		let read = reader.read_until(b'\n', &mut line)? as u64;
		if read == 0 {
			break;
		}
		let line_start = offset;
		offset += read;
		if let Some(header) = line.strip_prefix(b">") {
			entries.extend(current.take().map(|scan| scan.entry));
			let name = header
				.split(u8::is_ascii_whitespace)
				.next()
				.unwrap_or_default();
			let name = String::from_utf8(name.to_vec())
				.map_err(|_| Error::InvalidFasta { offset: line_start })?;
			current = Some(Scan {
				entry: FaiEntry {
					name,
					length: 0,
					offset,
					line_bases: 0,
					line_width: 0,
				},
				short_line_seen: false,
			});
			continue;
		}
		let scan = current
			.as_mut()
			.ok_or(Error::InvalidFasta { offset: line_start })?;
		let bases = line
			.iter()
			.take_while(|byte| !byte.is_ascii_whitespace())
			.count();
		if !matches!(line.get(bases..), Some(b"" | b"\n" | b"\r\n")) {
			return Err(Error::InvalidFasta { offset: line_start });
		}
		let bases = bases as u64;
		if bases == 0 {
			if scan.entry.line_bases == 0 {
				scan.entry.offset = offset;
			} else {
				scan.short_line_seen = true;
			}
			continue;
		}
		let terminated = line.ends_with(b"\n");
		if scan.entry.line_bases == 0 {
			scan.entry.line_bases = bases;
			scan.entry.line_width = read;
		} else if scan.short_line_seen
			|| bases > scan.entry.line_bases
			|| (terminated && read - bases != scan.entry.line_width - scan.entry.line_bases)
		{
			return Err(Error::InvalidFasta { offset: line_start });
		}
		if bases < scan.entry.line_bases {
			scan.short_line_seen = true;
		}
		scan.entry.length += bases;
	}
	entries.extend(current.map(|scan| scan.entry));
	for entry in &mut entries {
		if entry.line_bases == 0 {
			entry.line_bases = 1;
			entry.line_width = 1;
		}
	}
	Ok(entries)
}

#[cfg(unix)]
fn read_exact_at(file: &File, buffer: &mut [u8], offset: u64) -> std::io::Result<()> {
	std::os::unix::fs::FileExt::read_exact_at(file, buffer, offset)
}

#[cfg(windows)]
fn read_exact_at(file: &File, mut buffer: &mut [u8], mut offset: u64) -> std::io::Result<()> {
	while !buffer.is_empty() {
		match std::os::windows::fs::FileExt::seek_read(file, buffer, offset) {
			Ok(0) => return Err(std::io::ErrorKind::UnexpectedEof.into()),
			Ok(read) => {
				buffer = &mut buffer[read..];
				offset += read as u64;
			}
			Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
			Err(error) => return Err(error),
		}
	}
	Ok(())
}
