use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use gyrase_genome::{ContigId, Interval, Reference};
use gyrase_index::{Crai, DEFAULT_MAX_GAP, Index};
use gyrase_io::{Workers, decode};
use gyrase_reads::cram::{CramBatch, CramReader};
use gyrase_reads::{
	BamReader, Header, Line, Query, RecordBatch, RecordError, RecordRef, SamReader,
};

use crate::{Error, Input, Kind};

pub(crate) const EMPTY_COST: u64 = 1;
const INDEX_SUFFIXES: [&str; 2] = ["bai", "csi"];

#[derive(Debug)]
pub struct Alignments {
	path: PathBuf,
	header: Header,
	lookup: Lookup,
	reference: Arc<Reference>,
}

#[derive(Debug)]
enum Lookup {
	Bam(Index),
	Cram(Crai),
}

#[derive(Debug, Default)]
pub(crate) struct Batches {
	bam: Vec<RecordBatch>,
	cram: Vec<CramBatch>,
}

impl Alignments {
	pub fn open(
		input: &Input,
		reference: Arc<Reference>,
		workers: &Workers,
	) -> Result<Self, Error> {
		let path = &input.path;
		let lookup = match input.kind {
			Kind::Bam => Lookup::Bam(Index::read(
				open(&index_path(path, &INDEX_SUFFIXES)?)?,
				workers,
			)?),
			Kind::Cram => Lookup::Cram(Crai::read(open(&index_path(path, &["crai"])?)?, workers)?),
			_ => return Err(Error::NotAlignments(path.clone())),
		};
		Ok(Self {
			path: path.clone(),
			header: header(input, Some(Arc::clone(&reference)), workers)?,
			lookup,
			reference,
		})
	}

	pub fn header(&self) -> &Header {
		&self.header
	}

	pub fn reference(&self) -> &Reference {
		&self.reference
	}

	pub fn sample(&self) -> String {
		self.header
			.lines()
			.iter()
			.find_map(|line| match line {
				Line::ReadGroup(fields) => fields.get(*b"SM"),
				_ => None,
			})
			.map_or_else(
				|| {
					self.path
						.file_stem()
						.map_or_else(String::new, |stem| stem.to_string_lossy().into_owned())
				},
				str::to_owned,
			)
	}

	pub(crate) fn fetch(&self, interval: Interval, workers: &Workers) -> Result<Batches, Error> {
		let file = open(&self.path)?;
		let mut batches = Batches::default();
		match &self.lookup {
			Lookup::Bam(index) => {
				let chunks = index.plan(&[interval], DEFAULT_MAX_GAP).chunks().to_vec();
				let mut reader = BamReader::query(
					file,
					self.header.clone(),
					Query::new(chunks, vec![interval]),
					workers,
				)?;
				loop {
					let mut batch = RecordBatch::new();
					if !reader.read_batch(&mut batch)? {
						break;
					}
					batches.bam.push(batch);
				}
			}
			Lookup::Cram(crai) => {
				let mut reader = CramReader::query(
					file,
					Some(Arc::clone(&self.reference)),
					&name(&self.path),
					crai,
					interval,
					workers,
				)?;
				loop {
					let mut batch = CramBatch::new();
					if !reader.read_batch(&mut batch)? {
						break;
					}
					batches.cram.push(batch);
				}
			}
		}
		Ok(batches)
	}

	pub(crate) fn cost(&self, interval: Interval) -> u64 {
		let bytes: u64 = match &self.lookup {
			Lookup::Bam(index) => index
				.plan(&[interval], DEFAULT_MAX_GAP)
				.chunks()
				.iter()
				.map(|chunk| chunk.end.block().saturating_sub(chunk.begin.block()))
				.sum(),
			Lookup::Cram(crai) => crai
				.overlapping(interval)
				.map(|entry| entry.slice_size)
				.sum(),
		};
		bytes + EMPTY_COST
	}

	pub(crate) fn reference_contig(&self, contig: ContigId) -> Option<ContigId> {
		let name = &self.header.references().get(contig.index())?.name();
		self.reference.contigs().resolve(name)
	}
}

impl Batches {
	pub(crate) fn records(&self) -> impl Iterator<Item = Result<RecordRef<'_>, RecordError>> {
		let bam = self.bam.iter().flat_map(|batch| {
			batch
				.iter()
				.map(|located| located.map(|located| located.record))
		});
		let cram = self.cram.iter().flat_map(CramBatch::iter);
		bam.chain(cram)
	}
}

pub fn header(
	input: &Input,
	reference: Option<Arc<Reference>>,
	workers: &Workers,
) -> Result<Header, Error> {
	let file = open(&input.path)?;
	match input.kind {
		Kind::Bam => Ok(BamReader::read_header(&file, workers)?),
		Kind::Cram => Ok(
			CramReader::open(Box::new(file), reference, &name(&input.path), workers)?
				.header()
				.clone(),
		),
		Kind::Sam => {
			Ok(SamReader::new(
				decode(Box::new(file), workers).map_err(|source| Error::Read {
					path: input.path.clone(),
					source,
				})?,
			)?
			.header()
			.clone())
		}
		_ => Err(Error::NotAlignments(input.path.clone())),
	}
}

pub fn indexed(input: &Input) -> bool {
	match input.kind {
		Kind::Cram => index_path(&input.path, &["crai"]).is_ok(),
		_ => index_path(&input.path, &INDEX_SUFFIXES).is_ok(),
	}
}

fn open(path: &Path) -> Result<File, Error> {
	File::open(path).map_err(|source| Error::Open {
		path: path.to_owned(),
		source,
	})
}

fn index_path(path: &Path, suffixes: &[&str]) -> Result<PathBuf, Error> {
	suffixes
		.iter()
		.flat_map(|suffix| {
			let mut appended = path.as_os_str().to_owned();
			appended.push(".");
			appended.push(suffix);
			[PathBuf::from(appended), path.with_extension(suffix)]
		})
		.find(|candidate| candidate.is_file())
		.ok_or_else(|| Error::MissingIndex(path.to_owned()))
}

fn name(path: &Path) -> String {
	path.to_string_lossy().into_owned()
}
