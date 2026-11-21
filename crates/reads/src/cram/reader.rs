use std::fmt;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::sync::Arc;

use gyrase_genome::{Interval, Reference};
use gyrase_index::Crai;
use gyrase_io::Workers;
use rayon::prelude::*;

use crate::cram::Error;
use crate::cram::cursor::Cursor;
use crate::cram::decode::{Environment, decode_slice};
use crate::cram::emit::Decoded;
use crate::cram::slice;
use crate::cram::structure::{Block, BlockKind, ContainerHeader, Version};
use crate::{Header, RecordError, RecordRef};

const MAX_HEADER_TEXT: usize = 1 << 28;
const CONTAINERS_PER_WORKER: usize = 2;

#[derive(Debug, Default)]
pub struct CramBatch {
	slices: Vec<Decoded>,
	totals: Vec<usize>,
}

impl CramBatch {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn len(&self) -> usize {
		self.totals.last().copied().unwrap_or(0)
	}

	pub fn is_empty(&self) -> bool {
		self.len() == 0
	}

	pub fn at(&self, index: usize) -> Result<RecordRef<'_>, RecordError> {
		let slice = self.totals.partition_point(|total| *total <= index);
		let before = slice
			.checked_sub(1)
			.and_then(|previous| self.totals.get(previous))
			.copied()
			.unwrap_or(0);
		self.slices
			.get(slice)
			.ok_or(RecordError::Truncated)?
			.record(index - before)
	}

	pub fn iter(&self) -> impl Iterator<Item = Result<RecordRef<'_>, RecordError>> + Clone {
		(0..self.len()).map(|index| self.at(index))
	}

	fn reset(&mut self) {
		self.slices.clear();
		self.totals.clear();
	}

	fn push(&mut self, slice: Decoded) {
		if slice.len() > 0 {
			self.totals.push(self.len() + slice.len());
			self.slices.push(slice);
		}
	}
}

enum Input {
	Stream(Box<dyn Read + Send>),
	Offsets {
		file: File,
		offsets: Vec<u64>,
		next: usize,
	},
}

impl fmt::Debug for Input {
	fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
		match self {
			Self::Stream(_) => formatter.write_str("Stream"),
			Self::Offsets { offsets, next, .. } => formatter
				.debug_struct("Offsets")
				.field("offsets", offsets)
				.field("next", next)
				.finish_non_exhaustive(),
		}
	}
}

impl Input {
	fn next_container(&mut self) -> Result<Option<(ContainerHeader, Vec<u8>)>, Error> {
		let container = match self {
			Self::Stream(source) => ContainerHeader::read(source)?,
			Self::Offsets {
				file,
				offsets,
				next,
			} => {
				let Some(offset) = offsets.get(*next) else {
					return Ok(None);
				};
				*next += 1;
				file.seek(SeekFrom::Start(*offset))?;
				ContainerHeader::read(file)?
			}
		};
		Ok(container.filter(|(header, _)| !header.is_eof()))
	}
}

#[derive(Debug)]
pub struct CramReader {
	header: Header,
	environment: Arc<Environment>,
	input: Input,
	region: Option<Interval>,
	workers: Workers,
	current: CramBatch,
	cursor: usize,
}

impl CramReader {
	pub fn open(
		mut source: Box<dyn Read + Send>,
		reference: Option<Arc<Reference>>,
		name: &str,
		workers: &Workers,
	) -> Result<Self, Error> {
		let header = read_header(&mut source)?;
		Ok(Self::assemble(
			header,
			reference,
			name,
			Input::Stream(source),
			None,
			workers,
		))
	}

	pub fn query(
		mut file: File,
		reference: Option<Arc<Reference>>,
		name: &str,
		crai: &Crai,
		region: Interval,
		workers: &Workers,
	) -> Result<Self, Error> {
		let header = read_header(&mut file)?;
		let mut offsets: Vec<u64> = crai
			.overlapping(region)
			.map(|entry| entry.container_offset)
			.collect();
		offsets.sort_unstable();
		offsets.dedup();
		let input = Input::Offsets {
			file,
			offsets,
			next: 0,
		};
		Ok(Self::assemble(
			header,
			reference,
			name,
			input,
			Some(region),
			workers,
		))
	}

	fn assemble(
		header: Header,
		reference: Option<Arc<Reference>>,
		name: &str,
		input: Input,
		region: Option<Interval>,
		workers: &Workers,
	) -> Self {
		Self {
			environment: Arc::new(Environment::new(&header, reference, name)),
			header,
			input,
			region,
			workers: workers.clone(),
			current: CramBatch::new(),
			cursor: 0,
		}
	}

	pub fn header(&self) -> &Header {
		&self.header
	}

	pub fn read_batch(&mut self, batch: &mut CramBatch) -> Result<bool, Error> {
		let wanted = self.workers.threads().max(1) * CONTAINERS_PER_WORKER;
		let mut containers = Vec::with_capacity(wanted);
		while containers.len() < wanted {
			match self.input.next_container()? {
				Some(container) => containers.push(container),
				None => break,
			}
		}
		if containers.is_empty() {
			return Ok(false);
		}
		let environment = &self.environment;
		let region = self.region;
		let decoded: Vec<Result<Vec<Decoded>, Error>> = self.workers.install(|| {
			containers
				.par_iter()
				.map(|(header, body)| decode_container(environment, header, body, region))
				.collect()
		});
		batch.reset();
		for container in decoded {
			for slice in container? {
				batch.push(slice);
			}
		}
		Ok(true)
	}

	pub fn next_record(&mut self) -> Result<Option<RecordRef<'_>>, Error> {
		while self.cursor >= self.current.len() {
			let mut batch = std::mem::take(&mut self.current);
			let more = self.read_batch(&mut batch)?;
			self.current = batch;
			self.cursor = 0;
			if !more {
				return Ok(None);
			}
		}
		let index = self.cursor;
		self.cursor += 1;
		Ok(Some(self.current.at(index)?))
	}
}

fn read_header(source: &mut impl Read) -> Result<Header, Error> {
	Version::read(source)?;
	let (_, body) = ContainerHeader::read(source)?.ok_or(Error::Truncated)?;
	let block = Block::read(&mut Cursor::new(&body))?;
	if block.kind != BlockKind::FileHeader {
		return Err(Error::Malformed("file header"));
	}
	let raw = block.decode()?;
	let mut cursor = Cursor::new(&raw);
	let length = usize::try_from(cursor.i32_le()?)
		.ok()
		.filter(|length| *length <= MAX_HEADER_TEXT)
		.ok_or(Error::Malformed("file header"))?;
	Ok(Header::parse(cursor.take(length)?)?)
}

fn decode_container(
	environment: &Environment,
	header: &ContainerHeader,
	body: &[u8],
	region: Option<Interval>,
) -> Result<Vec<Decoded>, Error> {
	let parsed = slice::parse(header, body)?;
	let mut decoded = Vec::with_capacity(parsed.slices.len());
	for blocks in &parsed.slices {
		let mut slice = decode_slice(environment, &parsed.compression, blocks)?;
		if let Some(region) = region {
			slice.retain(|record| overlaps(record, region));
		}
		decoded.push(slice);
	}
	Ok(decoded)
}

fn overlaps(record: RecordRef<'_>, region: Interval) -> bool {
	let (Some(contig), Some(start), Some(end)) =
		(record.contig(), record.position(), record.reference_end())
	else {
		return false;
	};
	contig == region.contig && u64::from(start) < region.end && end > region.start
}
