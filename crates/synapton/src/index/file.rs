use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::Path;

use gyrase_genome::{Contig, Contigs};
use gyrase_io::Workers;
use rayon::prelude::*;

use super::Index;
use crate::Error;
use crate::config::SeedConfig;
use crate::hash::{checksum_chunk, checksum_words};
use crate::num::u32_of;
use crate::reference::{ContigSpan, PackedReference};

const MAGIC: [u8; 8] = *b"SYNAPTON";
const VERSION: u32 = 2;
const HEADER_LEN: usize = 128;
const HEADER_WORDS: usize = HEADER_LEN / 8;
const READ_CHUNK_BYTES: usize = 4 << 20;
const WRITE_CHUNK_WORDS: usize = 1 << 17;
const MAX_CONTIG_NAME: usize = 1 << 16;

pub(super) enum Reuse {
	Loaded(Box<Index>),
	Stale(PackedReference),
}

struct Header {
	config: SeedConfig,
	bucket_bits: u32,
	filter_cutoff: u32,
	reference_checksum: u64,
	reference_length: u64,
	meta_len: u64,
	word_count: u64,
	bucket_count: u64,
	entry_count: u64,
	meta_checksum: u64,
	words_checksum: u64,
	buckets_checksum: u64,
	entries_checksum: u64,
}

impl Header {
	fn encode(&self) -> [u8; HEADER_LEN] {
		let mut bytes = [0u8; HEADER_LEN];
		let mut cursor = 0;
		let mut put = |value: &[u8]| {
			if let Some(target) = bytes.get_mut(cursor..cursor + value.len()) {
				target.copy_from_slice(value);
			}
			cursor += value.len();
		};
		put(&MAGIC);
		put(&VERSION.to_le_bytes());
		for field in [
			self.config.k,
			self.config.s,
			self.config.w_min,
			self.config.w_max,
		] {
			put(&u32_of(field).to_le_bytes());
		}
		put(&self.bucket_bits.to_le_bytes());
		put(&self.filter_cutoff.to_le_bytes());
		put(&0u32.to_le_bytes());
		for field in [
			self.reference_checksum,
			self.reference_length,
			self.meta_len,
			self.word_count,
			self.bucket_count,
			self.entry_count,
			self.meta_checksum,
			self.words_checksum,
			self.buckets_checksum,
			self.entries_checksum,
		] {
			put(&field.to_le_bytes());
		}
		let words = header_words(&bytes);
		let checksum = checksum_chunk(words.get(..HEADER_WORDS - 1).unwrap_or_default());
		if let Some(target) = bytes.get_mut(HEADER_LEN - 8..) {
			target.copy_from_slice(&checksum.to_le_bytes());
		}
		bytes
	}

	fn decode(bytes: &[u8; HEADER_LEN]) -> Result<Self, Error> {
		if bytes.get(..8) != Some(&MAGIC) {
			return Err(Error::NotAnIndex);
		}
		let mut reader = ByteReader::new(bytes);
		reader.take(8)?;
		let version = reader.u32()?;
		if version != VERSION {
			return Err(Error::IndexVersion {
				found: version,
				expected: VERSION,
			});
		}
		let words = header_words(bytes);
		let stored = words.last().copied().unwrap_or(0);
		if checksum_chunk(words.get(..HEADER_WORDS - 1).unwrap_or_default()) != stored {
			return Err(Error::IndexChecksum("header"));
		}
		let k = reader.u32()? as usize;
		let s = reader.u32()? as usize;
		let w_min = reader.u32()? as usize;
		let w_max = reader.u32()? as usize;
		let config = SeedConfig::new(k, s, w_min, w_max)?;
		let bucket_bits = reader.u32()?;
		let filter_cutoff = reader.u32()?;
		reader.u32()?;
		Ok(Self {
			config,
			bucket_bits,
			filter_cutoff,
			reference_checksum: reader.u64()?,
			reference_length: reader.u64()?,
			meta_len: reader.u64()?,
			word_count: reader.u64()?,
			bucket_count: reader.u64()?,
			entry_count: reader.u64()?,
			meta_checksum: reader.u64()?,
			words_checksum: reader.u64()?,
			buckets_checksum: reader.u64()?,
			entries_checksum: reader.u64()?,
		})
	}

	fn words_offset(&self) -> u64 {
		HEADER_LEN as u64 + self.meta_len
	}

	fn buckets_offset(&self) -> u64 {
		self.words_offset() + self.word_count * 8
	}

	fn entries_offset(&self) -> u64 {
		self.buckets_offset() + self.bucket_count * 4
	}

	fn file_len(&self) -> u64 {
		self.entries_offset() + self.entry_count * 8
	}

	fn validate_sizes(&self, file_len: u64) -> Result<(), Error> {
		let sane = self.bucket_count == (1u64 << self.bucket_bits) + 1
			&& self.word_count < (1 << 40)
			&& self.bucket_count < (1 << 40)
			&& self.entry_count < (1 << 40)
			&& self.meta_len < (1 << 40);
		if !sane {
			return Err(Error::CorruptIndex("section sizes"));
		}
		if self.file_len() != file_len {
			return Err(Error::CorruptIndex("file length"));
		}
		Ok(())
	}
}

fn header_words(bytes: &[u8; HEADER_LEN]) -> Vec<u64> {
	bytes
		.as_chunks::<8>()
		.0
		.iter()
		.map(|chunk| u64::from_le_bytes(*chunk))
		.collect()
}

struct ByteReader<'a> {
	bytes: &'a [u8],
}

impl<'a> ByteReader<'a> {
	fn new(bytes: &'a [u8]) -> Self {
		Self { bytes }
	}

	fn take(&mut self, count: usize) -> Result<&'a [u8], Error> {
		if count > self.bytes.len() {
			return Err(Error::CorruptIndex("metadata"));
		}
		let (head, tail) = self.bytes.split_at(count);
		self.bytes = tail;
		Ok(head)
	}

	fn u32(&mut self) -> Result<u32, Error> {
		let bytes = self.take(4)?;
		Ok(u32::from_le_bytes(bytes.try_into().unwrap_or([0; 4])))
	}

	fn u64(&mut self) -> Result<u64, Error> {
		let bytes = self.take(8)?;
		Ok(u64::from_le_bytes(bytes.try_into().unwrap_or([0; 8])))
	}

	fn remaining(&self) -> usize {
		self.bytes.len()
	}
}

fn encode_meta(reference: &PackedReference) -> Vec<u8> {
	let mut meta = Vec::new();
	meta.extend_from_slice(&(reference.spans().len() as u64).to_le_bytes());
	for (span, (_, contig)) in reference.spans().iter().zip(reference.contigs().iter()) {
		meta.extend_from_slice(&u32_of(contig.name.len()).to_le_bytes());
		meta.extend_from_slice(contig.name.as_bytes());
		meta.extend_from_slice(&span.length.to_le_bytes());
		meta.extend_from_slice(&span.start.to_le_bytes());
	}
	meta.extend_from_slice(&(reference.ambiguous().len() as u64).to_le_bytes());
	for run in reference.ambiguous() {
		meta.extend_from_slice(&run.start.to_le_bytes());
		meta.extend_from_slice(&run.end.to_le_bytes());
	}
	meta
}

struct Meta {
	spans: Vec<ContigSpan>,
	contigs: Contigs,
	ambiguous: Vec<std::ops::Range<u64>>,
}

fn decode_meta(bytes: &[u8]) -> Result<Meta, Error> {
	let mut reader = ByteReader::new(bytes);
	let contig_count = reader.u64()?;
	if contig_count > reader.remaining() as u64 {
		return Err(Error::CorruptIndex("contig count"));
	}
	let mut spans = Vec::new();
	let mut contigs = Vec::new();
	for _ in 0..contig_count {
		let name_len = reader.u32()? as usize;
		if name_len > MAX_CONTIG_NAME {
			return Err(Error::CorruptIndex("contig name"));
		}
		let name = String::from_utf8(reader.take(name_len)?.to_vec())
			.map_err(|_| Error::CorruptIndex("contig name"))?;
		let length = reader.u64()?;
		let start = reader.u64()?;
		spans.push(ContigSpan { start, length });
		contigs.push(Contig {
			name,
			length,
			checksum: None,
		});
	}
	let run_count = reader.u64()?;
	if run_count > (reader.remaining() / 16) as u64 {
		return Err(Error::CorruptIndex("ambiguous runs"));
	}
	let mut ambiguous = Vec::new();
	for _ in 0..run_count {
		let start = reader.u64()?;
		let end = reader.u64()?;
		ambiguous.push(start..end);
	}
	Ok(Meta {
		spans,
		contigs: Contigs::new(contigs)?,
		ambiguous,
	})
}

fn checksum_bytes(bytes: &[u8]) -> u64 {
	let words: Vec<u64> = bytes
		.chunks(8)
		.map(|chunk| {
			let mut padded = [0u8; 8];
			if let Some(target) = padded.get_mut(..chunk.len()) {
				target.copy_from_slice(chunk);
			}
			u64::from_le_bytes(padded)
		})
		.collect();
	checksum_chunk(&words)
}

fn checksum_u32s(values: &[u32]) -> u64 {
	let chunks: Vec<u64> = values
		.par_chunks(1 << 23)
		.map(|chunk| {
			chunk.iter().fold(chunk.len() as u64, |hash, &value| {
				(hash.rotate_left(23) ^ u64::from(value)).wrapping_mul(0xbf58_476d_1ce4_e5b9)
			})
		})
		.collect();
	checksum_chunk(&chunks)
}

pub(super) fn save(index: &Index, path: &Path, workers: &Workers) -> Result<(), Error> {
	let meta = encode_meta(&index.reference);
	let (words_checksum, buckets_checksum, entries_checksum) = workers.install(|| {
		(
			checksum_words(index.reference.words()),
			checksum_u32s(&index.buckets),
			checksum_words(&index.entries),
		)
	});
	let header = Header {
		config: index.config,
		bucket_bits: index.bucket_bits,
		filter_cutoff: index.filter_cutoff,
		reference_checksum: index.reference.checksum(),
		reference_length: index.reference.length(),
		meta_len: meta.len() as u64,
		word_count: index.reference.words().len() as u64,
		bucket_count: index.buckets.len() as u64,
		entry_count: index.entries.len() as u64,
		meta_checksum: checksum_bytes(&meta),
		words_checksum,
		buckets_checksum,
		entries_checksum,
	};
	let mut partial = OsString::from(path);
	partial.push(".part");
	let written = (|| {
		let mut writer = BufWriter::with_capacity(8 << 20, File::create(&partial)?);
		writer.write_all(&header.encode())?;
		writer.write_all(&meta)?;
		write_u64s(&mut writer, index.reference.words())?;
		write_u32s(&mut writer, &index.buckets)?;
		write_u64s(&mut writer, &index.entries)?;
		writer
			.into_inner()
			.map_err(io::IntoInnerError::into_error)?
			.sync_all()
	})();
	if let Err(error) = written {
		let _ = fs::remove_file(&partial);
		return Err(error.into());
	}
	fs::rename(&partial, path)?;
	Ok(())
}

fn write_u64s(writer: &mut impl Write, values: &[u64]) -> io::Result<()> {
	let mut bytes = Vec::with_capacity(WRITE_CHUNK_WORDS * 8);
	for chunk in values.chunks(WRITE_CHUNK_WORDS) {
		bytes.clear();
		bytes.extend(chunk.iter().flat_map(|value| value.to_le_bytes()));
		writer.write_all(&bytes)?;
	}
	Ok(())
}

fn write_u32s(writer: &mut impl Write, values: &[u32]) -> io::Result<()> {
	let mut bytes = Vec::with_capacity(WRITE_CHUNK_WORDS * 4);
	for chunk in values.chunks(WRITE_CHUNK_WORDS) {
		bytes.clear();
		bytes.extend(chunk.iter().flat_map(|value| value.to_le_bytes()));
		writer.write_all(&bytes)?;
	}
	Ok(())
}

fn open(path: &Path) -> Result<(File, Header), Error> {
	let file = File::open(path)?;
	let file_len = file.metadata()?.len();
	if file_len < HEADER_LEN as u64 {
		return Err(Error::NotAnIndex);
	}
	let mut bytes = [0u8; HEADER_LEN];
	read_exact_at(&file, &mut bytes, 0)?;
	let header = Header::decode(&bytes)?;
	header.validate_sizes(file_len)?;
	Ok((file, header))
}

fn read_meta(file: &File, header: &Header) -> Result<Meta, Error> {
	let mut meta = vec![0u8; usize::try_from(header.meta_len)?];
	read_exact_at(file, &mut meta, HEADER_LEN as u64)?;
	if checksum_bytes(&meta) != header.meta_checksum {
		return Err(Error::IndexChecksum("metadata"));
	}
	decode_meta(&meta)
}

fn read_tables(
	file: &File,
	header: &Header,
	workers: &Workers,
) -> Result<(Vec<u32>, Vec<u64>), Error> {
	let buckets = read_u32s(
		file,
		header.buckets_offset(),
		usize::try_from(header.bucket_count)?,
		workers,
	)?;
	let entries = read_u64s(
		file,
		header.entries_offset(),
		usize::try_from(header.entry_count)?,
		workers,
	)?;
	let (buckets_ok, entries_ok) = workers.install(|| {
		(
			checksum_u32s(&buckets) == header.buckets_checksum,
			checksum_words(&entries) == header.entries_checksum,
		)
	});
	if !buckets_ok {
		return Err(Error::IndexChecksum("buckets"));
	}
	if !entries_ok {
		return Err(Error::IndexChecksum("entries"));
	}
	Ok((buckets, entries))
}

pub(super) fn load(path: &Path, workers: &Workers) -> Result<Index, Error> {
	let (file, header) = open(path)?;
	let meta = read_meta(&file, &header)?;
	let words = read_u64s(
		&file,
		header.words_offset(),
		usize::try_from(header.word_count)?,
		workers,
	)?;
	if workers.install(|| checksum_words(&words)) != header.words_checksum {
		return Err(Error::IndexChecksum("reference"));
	}
	let reference = PackedReference::from_parts(
		words,
		header.reference_length,
		meta.ambiguous,
		meta.spans,
		meta.contigs,
	)?;
	if reference.checksum() != header.reference_checksum {
		return Err(Error::IndexChecksum("reference identity"));
	}
	let (buckets, entries) = read_tables(&file, &header, workers)?;
	Ok(Index {
		config: header.config,
		bucket_bits: header.bucket_bits,
		filter_cutoff: header.filter_cutoff,
		buckets,
		entries,
		reference,
	})
}

pub(super) fn reuse(
	path: &Path,
	reference: PackedReference,
	config: SeedConfig,
	workers: &Workers,
) -> Reuse {
	let Ok((file, header)) = open(path) else {
		return Reuse::Stale(reference);
	};
	if header.reference_checksum != reference.checksum() || header.config != config {
		return Reuse::Stale(reference);
	}
	let Ok((buckets, entries)) = read_tables(&file, &header, workers) else {
		return Reuse::Stale(reference);
	};
	Reuse::Loaded(Box::new(Index {
		config,
		bucket_bits: header.bucket_bits,
		filter_cutoff: header.filter_cutoff,
		buckets,
		entries,
		reference,
	}))
}

fn read_u64s(file: &File, offset: u64, count: usize, workers: &Workers) -> Result<Vec<u64>, Error> {
	let mut values = vec![0u64; count];
	let per_chunk = READ_CHUNK_BYTES / 8;
	workers.install(|| {
		values
			.par_chunks_mut(per_chunk)
			.enumerate()
			.try_for_each_init(Vec::new, |bytes, (index, chunk)| {
				bytes.resize(chunk.len() * 8, 0);
				read_exact_at(file, bytes, offset + (index * per_chunk * 8) as u64)?;
				for (value, raw) in chunk.iter_mut().zip(bytes.as_chunks::<8>().0) {
					*value = u64::from_le_bytes(*raw);
				}
				Ok::<(), io::Error>(())
			})
	})?;
	Ok(values)
}

fn read_u32s(file: &File, offset: u64, count: usize, workers: &Workers) -> Result<Vec<u32>, Error> {
	let mut values = vec![0u32; count];
	let per_chunk = READ_CHUNK_BYTES / 4;
	workers.install(|| {
		values
			.par_chunks_mut(per_chunk)
			.enumerate()
			.try_for_each_init(Vec::new, |bytes, (index, chunk)| {
				bytes.resize(chunk.len() * 4, 0);
				read_exact_at(file, bytes, offset + (index * per_chunk * 4) as u64)?;
				for (value, raw) in chunk.iter_mut().zip(bytes.as_chunks::<4>().0) {
					*value = u32::from_le_bytes(*raw);
				}
				Ok::<(), io::Error>(())
			})
	})?;
	Ok(values)
}

#[cfg(unix)]
fn read_exact_at(file: &File, buffer: &mut [u8], offset: u64) -> io::Result<()> {
	std::os::unix::fs::FileExt::read_exact_at(file, buffer, offset)
}

#[cfg(windows)]
fn read_exact_at(file: &File, mut buffer: &mut [u8], mut offset: u64) -> io::Result<()> {
	while !buffer.is_empty() {
		match std::os::windows::fs::FileExt::seek_read(file, buffer, offset) {
			Ok(0) => return Err(io::ErrorKind::UnexpectedEof.into()),
			Ok(read) => {
				buffer = buffer.split_at_mut(read).1;
				offset += read as u64;
			}
			Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
			Err(error) => return Err(error),
		}
	}
	Ok(())
}
