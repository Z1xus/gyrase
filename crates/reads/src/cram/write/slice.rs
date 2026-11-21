use gyrase_genome::ContigId;

use super::header::compression_header;
use super::layout::{BlockKind, ContainerFields, block, put_itf8, put_ltf8, raw_block};
use super::record::RecordEncoder;
use crate::RecordRef;
use crate::cram::Error;
use crate::cram::decode::Environment;

const NO_EMBEDDED_REFERENCE: i32 = -1;
const NO_REFERENCE: i32 = -1;

#[derive(Debug, Default)]
pub(super) struct Batch {
	bytes: Vec<u8>,
	ends: Vec<usize>,
	pub(super) contig: Option<ContigId>,
}

impl Batch {
	pub(super) fn len(&self) -> usize {
		self.ends.len()
	}

	pub(super) fn push(&mut self, record: RecordRef<'_>) {
		let (core, rest) = record.parts();
		self.bytes.extend_from_slice(core);
		self.bytes.extend_from_slice(rest);
		self.ends.push(self.bytes.len());
	}

	fn records(&self) -> impl Iterator<Item = Result<RecordRef<'_>, Error>> {
		let mut start = 0;
		self.ends.iter().map(move |end| {
			let bytes = self.bytes.get(start..*end).unwrap_or_default();
			start = *end;
			Ok(RecordRef::new(bytes)?)
		})
	}
}

#[derive(Debug)]
pub(super) struct Encoded {
	pub(super) fields: ContainerFields,
	pub(super) body: Vec<u8>,
	pub(super) slice_offset: usize,
}

struct Extent {
	start: i32,
	span: i32,
	bases: u64,
}

fn extent(batch: &Batch) -> Result<Extent, Error> {
	let mut start = i64::MAX;
	let mut end = 0i64;
	let mut bases = 0u64;
	for record in batch.records() {
		let record = record?;
		bases += record.sequence().len() as u64;
		let Some(position) = record.position() else {
			continue;
		};
		let first = i64::from(position) + 1;
		let last = i64::try_from(record.reference_end().unwrap_or(0))
			.map_err(|_| Error::IntegerOverflow)?;
		start = start.min(first);
		end = end.max(last.max(first));
	}
	if start == i64::MAX {
		return Ok(Extent {
			start: 0,
			span: 0,
			bases,
		});
	}
	Ok(Extent {
		start: i32::try_from(start).map_err(|_| Error::IntegerOverflow)?,
		span: i32::try_from(end - start + 1).map_err(|_| Error::IntegerOverflow)?,
		bases,
	})
}

pub(super) fn encode(
	environment: &Environment,
	batch: &Batch,
	counter: u64,
) -> Result<Encoded, Error> {
	let extent = extent(batch)?;
	let reference = match batch.contig {
		Some(contig) => i32::try_from(contig.get()).map_err(|_| Error::IntegerOverflow)?,
		None => NO_REFERENCE,
	};
	let mut checksum = [0u8; 16];
	let window = match batch.contig {
		Some(contig) => {
			let start = u64::try_from(extent.start - 1).map_err(|_| Error::IntegerOverflow)?;
			let end = start + u64::try_from(extent.span).map_err(|_| Error::IntegerOverflow)?;
			checksum = environment
				.references
				.checksum(contig.index(), start, end)?;
			Some(environment.references.window(contig.index(), start, end)?)
		}
		None => None,
	};
	let mut encoder = RecordEncoder::new(environment, window.as_ref(), extent.start);
	for record in batch.records() {
		encoder.encode(record?)?;
	}
	let dictionary = encoder.dictionary;
	let header = compression_header(&dictionary.lines, &dictionary.keys);
	let mut body = Vec::new();
	block(&mut body, BlockKind::CompressionHeader, 0, &header);
	let slice_offset = body.len();
	let mut externals = Vec::new();
	let mut content_ids = Vec::new();
	for (id, bytes) in encoder.streams.into_blocks() {
		content_ids.push(id);
		block(&mut externals, BlockKind::ExternalData, id, &bytes);
	}
	let mut slice = Vec::new();
	put_itf8(&mut slice, reference);
	put_itf8(&mut slice, extent.start);
	put_itf8(&mut slice, extent.span);
	put_itf8(
		&mut slice,
		i32::try_from(batch.len()).map_err(|_| Error::IntegerOverflow)?,
	);
	put_ltf8(&mut slice, counter);
	let block_count = i32::try_from(content_ids.len() + 1).map_err(|_| Error::IntegerOverflow)?;
	put_itf8(&mut slice, block_count);
	put_itf8(&mut slice, block_count - 1);
	for id in &content_ids {
		put_itf8(&mut slice, *id);
	}
	put_itf8(&mut slice, NO_EMBEDDED_REFERENCE);
	slice.extend_from_slice(&checksum);
	raw_block(&mut body, BlockKind::SliceHeader, 0, &slice);
	raw_block(&mut body, BlockKind::CoreData, 0, &[]);
	body.extend_from_slice(&externals);
	Ok(Encoded {
		fields: ContainerFields {
			reference,
			start: extent.start,
			span: extent.span,
			records: u32::try_from(batch.len()).map_err(|_| Error::IntegerOverflow)?,
			counter,
			bases: extent.bases,
			blocks: block_count + 2,
			landmarks: vec![u32::try_from(slice_offset).map_err(|_| Error::IntegerOverflow)?],
		},
		body,
		slice_offset,
	})
}
