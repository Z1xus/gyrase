mod header;
mod layout;
mod record;
mod series;
mod slice;

use std::io::Write;
use std::sync::Arc;

use gyrase_genome::{ContigId, Reference};
use gyrase_index::{Crai, CraiEntry, CraiReference};
use gyrase_io::Workers;
use rayon::prelude::*;

use self::layout::{BlockKind, ContainerFields, EOF_START, container_header, raw_block};
use self::slice::{Batch, Encoded};
use crate::cram::Error;
use crate::cram::decode::Environment;
use crate::{Header, RecordRef};

const RECORDS_PER_SLICE: usize = 10_000;
const SLICES_PER_WORKER: usize = 2;
const MAGIC: [u8; 4] = *b"CRAM";
const VERSION: [u8; 2] = [3, 0];
const FILE_ID_LEN: usize = 20;
const EMPTY_COMPRESSION_HEADER: [u8; 6] = [1, 0, 1, 0, 1, 0];

#[derive(Debug)]
pub struct Finished<W> {
	pub sink: W,
	pub index: Crai,
}

#[derive(Debug)]
pub struct CramWriter<W: Write> {
	sink: W,
	environment: Arc<Environment>,
	workers: Workers,
	batch: Batch,
	pending: Vec<(Batch, u64)>,
	counter: u64,
	offset: u64,
	entries: Vec<CraiEntry>,
}

impl<W: Write> CramWriter<W> {
	pub fn new(
		mut sink: W,
		header: &Header,
		reference: Arc<Reference>,
		workers: &Workers,
	) -> Result<Self, Error> {
		let mut start = Vec::new();
		start.extend_from_slice(&MAGIC);
		start.extend_from_slice(&VERSION);
		start.extend_from_slice(&[0; FILE_ID_LEN]);
		let text = header.text();
		let mut payload = Vec::new();
		payload.extend_from_slice(
			&i32::try_from(text.len())
				.map_err(|_| Error::IntegerOverflow)?
				.to_le_bytes(),
		);
		payload.extend_from_slice(text.as_bytes());
		let mut body = Vec::new();
		raw_block(&mut body, BlockKind::FileHeader, 0, &payload);
		let fields = ContainerFields {
			reference: 0,
			start: 0,
			span: 0,
			records: 0,
			counter: 0,
			bases: 0,
			blocks: 1,
			landmarks: Vec::new(),
		};
		start.extend_from_slice(&container_header(&fields, body.len()));
		start.extend_from_slice(&body);
		sink.write_all(&start)?;
		Ok(Self {
			sink,
			environment: Arc::new(Environment::new(header, Some(reference), "")),
			workers: workers.clone(),
			batch: Batch::default(),
			pending: Vec::new(),
			counter: 0,
			offset: start.len() as u64,
			entries: Vec::new(),
		})
	}

	pub fn write(&mut self, record: RecordRef<'_>) -> Result<(), Error> {
		let contig = record.contig();
		if self.batch.len() >= RECORDS_PER_SLICE
			|| (self.batch.len() > 0 && self.batch.contig != contig)
		{
			self.close_batch()?;
		}
		self.batch.contig = contig;
		self.batch.push(record);
		Ok(())
	}

	pub fn finish(mut self) -> Result<Finished<W>, Error> {
		self.close_batch()?;
		self.flush()?;
		let mut body = Vec::new();
		raw_block(
			&mut body,
			BlockKind::CompressionHeader,
			0,
			&EMPTY_COMPRESSION_HEADER,
		);
		let fields = ContainerFields {
			reference: -1,
			start: EOF_START,
			span: 0,
			records: 0,
			counter: 0,
			bases: 0,
			blocks: 1,
			landmarks: Vec::new(),
		};
		self.sink
			.write_all(&container_header(&fields, body.len()))?;
		self.sink.write_all(&body)?;
		self.sink.flush()?;
		Ok(Finished {
			sink: self.sink,
			index: Crai::new(self.entries),
		})
	}

	fn close_batch(&mut self) -> Result<(), Error> {
		if self.batch.len() == 0 {
			return Ok(());
		}
		let batch = std::mem::take(&mut self.batch);
		let records = batch.len() as u64;
		self.pending.push((batch, self.counter));
		self.counter += records;
		if self.pending.len() >= self.workers.threads().max(1) * SLICES_PER_WORKER {
			self.flush()?;
		}
		Ok(())
	}

	fn flush(&mut self) -> Result<(), Error> {
		let pending = std::mem::take(&mut self.pending);
		let environment = &self.environment;
		let encoded: Vec<Result<Encoded, Error>> = self.workers.install(|| {
			pending
				.par_iter()
				.map(|(batch, counter)| slice::encode(environment, batch, *counter))
				.collect()
		});
		for container in encoded {
			self.write_container(&container?)?;
		}
		Ok(())
	}

	fn write_container(&mut self, container: &Encoded) -> Result<(), Error> {
		let header = container_header(&container.fields, container.body.len());
		self.sink.write_all(&header)?;
		self.sink.write_all(&container.body)?;
		let reference = u32::try_from(container.fields.reference)
			.map_or(CraiReference::Unmapped, |id| {
				CraiReference::Contig(ContigId::new(id))
			});
		self.entries.push(CraiEntry {
			reference,
			start: u64::try_from(container.fields.start).unwrap_or(0),
			span: u64::try_from(container.fields.span).unwrap_or(0),
			container_offset: self.offset,
			slice_offset: container.slice_offset as u64,
			slice_size: (container.body.len() - container.slice_offset) as u64,
		});
		self.offset += (header.len() + container.body.len()) as u64;
		Ok(())
	}
}
