use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::sync::Arc;

use gyrase_genome::{Interval, Reference};
use gyrase_index::{Crai, DEFAULT_MAX_GAP, Index};
use gyrase_io::{Workers, decode, peek};
use gyrase_reads::cram::{CramBatch, CramReader};
use gyrase_reads::{
	BamReader, Flagstat, Header, Query, Record, RecordBatch, RecordRef, SamReader, flagstat,
};

use crate::failure::{Failure, file};
use crate::shared::open;

const BAM_MAGIC: &[u8] = b"BAM\x01";
const CRAM_MAGIC: &[u8] = b"CRAM";
const MAGIC_LENGTH: u64 = 4;
const BUFFER: usize = 1 << 20;

enum Source {
	Sam(SamReader<Box<dyn BufRead>>),
	Bam(BamReader),
	Cram(CramReader),
	BamRegions(BamRegions),
	CramRegions(CramRegions),
}

struct BamRegions {
	path: String,
	header: Header,
	index: Index,
	regions: Vec<Interval>,
	multi: bool,
	workers: Workers,
}

struct CramRegions {
	path: String,
	multi: bool,
	crai: Crai,
	regions: Vec<Interval>,
	reference: Option<Arc<Reference>>,
	workers: Workers,
}

pub(super) struct Input {
	pub(super) header: Header,
	source: Source,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Kind {
	Sam,
	Bam,
	Cram,
	Unknown,
}

pub(super) fn kind(raw: impl Read) -> Result<Kind, Failure> {
	let (_, head) = peek(raw, MAGIC_LENGTH)?;
	Ok(if head.starts_with(BAM_MAGIC) {
		Kind::Bam
	} else if head.starts_with(CRAM_MAGIC) {
		Kind::Cram
	} else if head.starts_with(b"@") {
		Kind::Sam
	} else {
		Kind::Unknown
	})
}

pub(super) fn reference(path: &str) -> Result<Arc<Reference>, Failure> {
	Ok(Arc::new(Reference::open(Path::new(path))?))
}

pub(super) enum Need {
	Sequences,
	Positions,
}

impl Input {
	pub(super) fn open(
		path: &str,
		reference: Option<&Arc<Reference>>,
		need: Need,
		workers: &Workers,
	) -> Result<Self, Failure> {
		let mut buffered = BufReader::with_capacity(BUFFER, open(path)?);
		match kind(buffered.fill_buf()?)? {
			Kind::Bam => {
				let reader = BamReader::open(Box::new(buffered), workers)?;
				Ok(Self {
					header: reader.header().clone(),
					source: Source::Bam(reader),
				})
			}
			Kind::Cram => {
				let reader = match (reference, need) {
					(None, Need::Positions) => {
						CramReader::open_unbacked(Box::new(buffered), workers)?
					}
					_ => CramReader::open(Box::new(buffered), reference.cloned(), path, workers)?,
				};
				Ok(Self {
					header: reader.header().clone(),
					source: Source::Cram(reader),
				})
			}
			Kind::Sam | Kind::Unknown => {
				let text = decode(Box::new(buffered), workers)?;
				let reader = SamReader::new(Box::new(text) as Box<dyn BufRead>)?;
				Ok(Self {
					header: reader.header().clone(),
					source: Source::Sam(reader),
				})
			}
		}
	}

	pub(super) fn query(
		path: &str,
		regions: &[String],
		multi: bool,
		reference: Option<&Arc<Reference>>,
		workers: &Workers,
	) -> Result<Self, Failure> {
		if path == "-" {
			return Err(unindexed());
		}
		let source = File::open(path).map_err(|error| file(path, "reading", &error))?;
		match kind(source)? {
			Kind::Bam => bam_query(path, regions, multi, workers),
			Kind::Cram => cram_query(path, regions, multi, reference, workers),
			Kind::Sam | Kind::Unknown => Err(unindexed()),
		}
	}

	pub(super) fn is_bam(&self) -> bool {
		matches!(self.source, Source::Bam(_))
	}

	pub(super) fn flagstat(self, workers: &Workers) -> Result<Flagstat, Failure> {
		let mut stats = Flagstat::default();
		match self.source {
			Source::Bam(mut reader) => return Ok(flagstat(&mut reader, workers)?),
			source => visit(source, |record| {
				stats.add(record);
				Ok(())
			})?,
		}
		Ok(stats)
	}

	pub(super) fn each(
		self,
		visit_record: impl FnMut(RecordRef<'_>) -> Result<(), Failure>,
	) -> Result<(), Failure> {
		visit(self.source, visit_record)
	}
}

fn visit(
	source: Source,
	mut visit: impl FnMut(RecordRef<'_>) -> Result<(), Failure>,
) -> Result<(), Failure> {
	match source {
		Source::Sam(mut reader) => {
			let mut record = Record::new();
			while reader.read_record(&mut record)? {
				visit(record.view())?;
			}
		}
		Source::Bam(reader) => visit_bam(reader, &mut visit)?,
		Source::Cram(reader) => visit_cram(reader, &mut visit)?,
		Source::BamRegions(query) => {
			let groups = if query.multi {
				vec![query.regions.clone()]
			} else {
				query.regions.iter().map(|region| vec![*region]).collect()
			};
			for group in groups {
				let plan = query.index.plan(&group, DEFAULT_MAX_GAP);
				let query_set = Query::new(plan.chunks().to_vec(), group);
				let reader = BamReader::query(
					File::open(&query.path)?,
					query.header.clone(),
					query_set,
					&query.workers,
				)?;
				visit_bam(reader, &mut visit)?;
			}
		}
		Source::CramRegions(query) => {
			for (position, region) in query.regions.iter().enumerate() {
				let earlier = if query.multi {
					query.regions.get(..position).unwrap_or_default()
				} else {
					&[]
				};
				let reader = CramReader::query(
					File::open(&query.path)?,
					query.reference.clone(),
					&query.path,
					&query.crai,
					*region,
					&query.workers,
				)?;
				visit_cram(reader, &mut |record| {
					if overlaps_any(record, earlier) {
						return Ok(());
					}
					visit(record)
				})?;
			}
		}
	}
	Ok(())
}

fn overlaps_any(record: RecordRef<'_>, regions: &[Interval]) -> bool {
	let (Some(contig), Some(start), Some(end)) =
		(record.contig(), record.position(), record.reference_end())
	else {
		return false;
	};
	regions.iter().any(|region| {
		region.contig == contig && u64::from(start) < region.end && end > region.start
	})
}

fn visit_bam(
	mut reader: BamReader,
	visit: &mut impl FnMut(RecordRef<'_>) -> Result<(), Failure>,
) -> Result<(), Failure> {
	let mut batch = RecordBatch::new();
	while reader.read_batch(&mut batch)? {
		for located in batch.iter() {
			visit(located?.record)?;
		}
	}
	Ok(())
}

fn visit_cram(
	mut reader: CramReader,
	visit: &mut impl FnMut(RecordRef<'_>) -> Result<(), Failure>,
) -> Result<(), Failure> {
	let mut batch = CramBatch::new();
	while reader.read_batch(&mut batch)? {
		for record in batch.iter() {
			visit(record?)?;
		}
	}
	Ok(())
}

fn unindexed() -> Failure {
	Failure::new("random alignment retrieval only works for indexed BAM or CRAM files")
}

fn intervals(header: &Header, regions: &[String]) -> Result<Vec<Interval>, Failure> {
	let contigs = header.contigs()?;
	regions
		.iter()
		.map(|region| {
			Interval::parse(region, &contigs).map_err(|_| {
				Failure::new(format!(
					"region \"{region}\" specifies an unknown reference name"
				))
			})
		})
		.collect()
}

pub(super) fn sibling(path: &str, suffixes: &[&str]) -> Option<String> {
	let stem = path.rsplit_once('.').map(|(stem, _)| stem);
	suffixes.iter().find_map(|suffix| {
		[
			Some(format!("{path}{suffix}")),
			stem.map(|stem| format!("{stem}{suffix}")),
		]
		.into_iter()
		.flatten()
		.find(|candidate| Path::new(candidate).exists())
	})
}

fn bam_query(
	path: &str,
	regions: &[String],
	multi: bool,
	workers: &Workers,
) -> Result<Input, Failure> {
	let index_path = sibling(path, &[".bai", ".csi"]).ok_or_else(unindexed)?;
	let header = BamReader::read_header(&File::open(path)?, workers)?;
	let regions = intervals(&header, regions)?;
	Ok(Input {
		source: Source::BamRegions(BamRegions {
			path: path.to_owned(),
			header: header.clone(),
			index: Index::read(File::open(index_path)?, workers)?,
			regions,
			multi,
			workers: workers.clone(),
		}),
		header,
	})
}

fn cram_query(
	path: &str,
	regions: &[String],
	multi: bool,
	reference: Option<&Arc<Reference>>,
	workers: &Workers,
) -> Result<Input, Failure> {
	let crai_path = sibling(path, &[".crai"]).ok_or_else(unindexed)?;
	let header = CramReader::open(Box::new(File::open(path)?), None, path, workers)?
		.header()
		.clone();
	let mut regions = intervals(&header, regions)?;
	if multi {
		regions.sort_unstable();
		regions = merged(regions);
	}
	Ok(Input {
		source: Source::CramRegions(CramRegions {
			path: path.to_owned(),
			multi,
			crai: Crai::read(File::open(crai_path)?, workers)?,
			regions,
			reference: reference.cloned(),
			workers: workers.clone(),
		}),
		header,
	})
}

fn merged(sorted: Vec<Interval>) -> Vec<Interval> {
	let mut merged: Vec<Interval> = Vec::with_capacity(sorted.len());
	for interval in sorted {
		match merged.last_mut() {
			Some(last) if last.contig == interval.contig && interval.start <= last.end => {
				last.end = last.end.max(interval.end);
			}
			_ => merged.push(interval),
		}
	}
	merged
}
