use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use gyrase_genome::Reference;
use gyrase_index::{Format, Index};
use gyrase_io::{Level, Workers, decode};
use gyrase_reads::cram::{CramBatch, CramReader, index_cram};
use gyrase_reads::{
	BamReader, Header, Line, MarkSettings, MarkStats, Record, RecordBatch, SamReader, Sorter,
	index_bam, index_layout,
};

use crate::counted::Counted;
use crate::partial::{Partial, suffixed};
use crate::{Cancel, Error, Estimate, Input, Kind, Progress, Stage};

const MIN_SORT_MEMORY: u64 = 256 << 20;
const MAX_SORT_MEMORY: u64 = 4 << 30;
const SORT_MEMORY_SHARE: u64 = 4;
const COORDINATE: &str = "coordinate";
const RECORDS_PER_REPORT: u32 = 1 << 16;

pub(crate) fn coordinate_sorted(header: &Header) -> bool {
	header.lines().iter().any(|line| match line {
		Line::Header(fields) => fields.get(*b"SO") == Some(COORDINATE),
		_ => false,
	})
}

pub(crate) fn sort_budget(memory: u64) -> usize {
	usize::try_from((memory / SORT_MEMORY_SHARE).clamp(MIN_SORT_MEMORY, MAX_SORT_MEMORY))
		.unwrap_or(usize::MAX)
}

pub(crate) fn index(
	input: &Input,
	workers: &Workers,
	cancel: &Cancel,
	report: &mut impl FnMut(Progress),
) -> Result<(), Error> {
	let mut estimate = Estimate::new(Stage::Index, input.size, Instant::now());
	report(estimate.progress());
	let advance = |position| {
		cancel.check()?;
		report(estimate.advance(position, Instant::now()));
		Ok::<_, Error>(())
	};
	if input.kind == Kind::Cram {
		let crai = index_cram(BufReader::new(open(&input.path)?), advance)?;
		let (partial, file) = Partial::create(&suffixed(&input.path, ".crai"))?;
		crai.write(file, workers)?;
		return partial.commit();
	}
	let mut reader = BamReader::open(Box::new(open(&input.path)?), workers)?;
	let (format, geometry) = index_layout(reader.header())?;
	let index = index_bam(&mut reader, geometry, advance)?;
	write_index(&index, format, &input.path, workers)
}

pub(crate) struct Sorting<'a> {
	pub reference: &'a Arc<Reference>,
	pub budget: usize,
	pub workers: &'a Workers,
	pub cancel: &'a Cancel,
}

impl Sorting<'_> {
	pub(crate) fn run(
		self,
		input: &Input,
		output: &Path,
		report: &mut impl FnMut(Progress),
	) -> Result<(), Error> {
		let read = Arc::new(AtomicU64::new(0));
		let source = Box::new(Counted::new(open(&input.path)?, &read));
		let total = input.size.saturating_mul(2);
		let mut estimate = Estimate::new(Stage::Sort, total, Instant::now());
		report(estimate.progress());
		let mut advance = || {
			self.cancel.check()?;
			report(estimate.advance(read.load(Ordering::Relaxed), Instant::now()));
			Ok::<_, Error>(())
		};
		match input.kind {
			Kind::Cram => {
				let name = input.path.to_string_lossy();
				let reference = Some(Arc::clone(self.reference));
				let mut reader = CramReader::open(source, reference, &name, self.workers)?;
				let header = reader.header().clone();
				self.write(&header, output, |sorter| {
					let mut batch = CramBatch::new();
					while reader.read_batch(&mut batch)? {
						for record in batch.iter() {
							sorter.push(record?)?;
						}
						advance()?;
					}
					Ok(())
				})?;
			}
			Kind::Sam => {
				let mut reader =
					SamReader::new(decode(source, self.workers).map_err(|source| {
						Error::Read {
							path: input.path.clone(),
							source,
						}
					})?)?;
				let header = reader.header().clone();
				self.write(&header, output, |sorter| {
					let mut record = Record::new();
					let mut count = 0_u32;
					while reader.read_record(&mut record)? {
						sorter.push(record.view())?;
						count = count.wrapping_add(1);
						if count.is_multiple_of(RECORDS_PER_REPORT) {
							advance()?;
						}
					}
					Ok(())
				})?;
			}
			_ => {
				let mut reader = BamReader::open(source, self.workers)?;
				let header = reader.header().clone();
				self.write(&header, output, |sorter| {
					let mut batch = RecordBatch::new();
					while reader.read_batch(&mut batch)? {
						for located in batch.iter() {
							sorter.push(located.map_err(gyrase_reads::Error::from)?.record)?;
						}
						advance()?;
					}
					Ok(())
				})?;
			}
		}
		report(estimate.advance(total, Instant::now()));
		Ok(())
	}

	fn write(
		&self,
		header: &Header,
		output: &Path,
		feed: impl FnOnce(&mut Sorter<'_>) -> Result<(), Error>,
	) -> Result<(), Error> {
		write_sorted(
			header,
			output,
			self.budget,
			None,
			self.workers,
			self.cancel,
			feed,
		)?;
		Ok(())
	}
}

pub(crate) fn write_sorted(
	header: &Header,
	output: &Path,
	budget: usize,
	marking: Option<MarkSettings>,
	workers: &Workers,
	cancel: &Cancel,
	feed: impl FnOnce(&mut Sorter<'_>) -> Result<(), Error>,
) -> Result<Option<MarkStats>, Error> {
	let directory = output
		.parent()
		.filter(|parent| !parent.as_os_str().is_empty())
		.unwrap_or(Path::new("."));
	let mut sorter = Sorter::new(header, budget, directory, workers, cancel.flag())?;
	let loaded = feed(&mut sorter);
	cancel.check()?;
	loaded?;
	let (partial, sink) = Partial::create(output)?;
	let finished = match marking {
		Some(settings) => sorter
			.finish_marking(sink, Level::DEFAULT, settings)
			.map(|(written, stats)| (written, Some(stats))),
		None => sorter
			.finish(sink, Level::DEFAULT)
			.map(|written| (written, None)),
	};
	cancel.check()?;
	let (written, stats) = finished?;
	drop(written.output);
	partial.commit()?;
	write_index(&written.index, written.format, output, workers)?;
	Ok(stats)
}

fn write_index(index: &Index, format: Format, bam: &Path, workers: &Workers) -> Result<(), Error> {
	let path = suffixed(
		bam,
		match format {
			Format::Csi => ".csi",
			Format::Bai | Format::Tbi => ".bai",
		},
	);
	let (partial, file) = Partial::create(&path)?;
	index.write(format, file, workers)?;
	partial.commit()
}

fn open(path: &Path) -> Result<File, Error> {
	File::open(path).map_err(|source| Error::Open {
		path: path.to_owned(),
		source,
	})
}
