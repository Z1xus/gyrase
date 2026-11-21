use std::ffi::OsString;
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::time::Instant;

use gyrase_index::{Format, Index};
use gyrase_io::{Level, Workers};
use gyrase_reads::{
	BamReader, Header, Line, MarkSettings, MarkStats, RecordBatch, Sorter, index_bam, index_layout,
};

use crate::{Cancel, Error, Estimate, Input, Progress, Stage};

const MIN_SORT_MEMORY: u64 = 256 << 20;
const MAX_SORT_MEMORY: u64 = 4 << 30;
const SORT_MEMORY_SHARE: u64 = 4;
const COORDINATE: &str = "coordinate";

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
	let mut reader = BamReader::open(Box::new(open(&input.path)?), workers)?;
	let (format, geometry) = index_layout(reader.header())?;
	let mut estimate = Estimate::new(Stage::Index, input.size, Instant::now());
	report(estimate.progress());
	let index = index_bam(&mut reader, geometry, |position| {
		cancel.check()?;
		report(estimate.advance(position, Instant::now()));
		Ok::<_, Error>(())
	})?;
	write_index(&index, format, &input.path, workers)
}

pub(crate) fn sort(
	input: &Input,
	output: &Path,
	budget: usize,
	workers: &Workers,
	cancel: &Cancel,
	report: &mut impl FnMut(Progress),
) -> Result<(), Error> {
	let mut reader = BamReader::open(Box::new(open(&input.path)?), workers)?;
	let header = reader.header().clone();
	let total = input.size.saturating_mul(2);
	let mut estimate = Estimate::new(Stage::Sort, total, Instant::now());
	report(estimate.progress());
	write_sorted(&header, output, budget, None, workers, cancel, |sorter| {
		let mut batch = RecordBatch::new();
		while reader.read_batch(&mut batch)? {
			let mut position = 0;
			for located in batch.iter() {
				let located = located.map_err(gyrase_reads::Error::from)?;
				position = located.end.block();
				sorter.push(located.record)?;
			}
			cancel.check()?;
			report(estimate.advance(position, Instant::now()));
		}
		Ok(())
	})?;
	report(estimate.advance(total, Instant::now()));
	Ok(())
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
	let mut sorter = Sorter::new(header, budget, directory, workers)?;
	feed(&mut sorter)?;
	cancel.check()?;
	let partial = suffixed(output, ".part");
	let file = create(&partial)?;
	let (written, stats) = match marking {
		Some(settings) => {
			let (written, stats) = sorter.finish_marking(file, Level::DEFAULT, settings)?;
			(written, Some(stats))
		}
		None => (sorter.finish(file, Level::DEFAULT)?, None),
	};
	rename(&partial, output)?;
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
	let partial = suffixed(&path, ".part");
	index.write(format, create(&partial)?, workers)?;
	rename(&partial, &path)
}

pub(crate) fn suffixed(path: &Path, suffix: &str) -> PathBuf {
	let mut name = OsString::from(path);
	name.push(suffix);
	PathBuf::from(name)
}

fn open(path: &Path) -> Result<File, Error> {
	File::open(path).map_err(|source| Error::Open {
		path: path.to_owned(),
		source,
	})
}

pub(crate) fn create(path: &Path) -> Result<File, Error> {
	File::create(path).map_err(|source| Error::Open {
		path: path.to_owned(),
		source,
	})
}

fn rename(from: &Path, to: &Path) -> Result<(), Error> {
	fs::rename(from, to).map_err(|source| Error::Open {
		path: to.to_owned(),
		source,
	})
}
