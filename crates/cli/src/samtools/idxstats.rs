use std::fmt::Write as _;
use std::fs::File;
use std::io::Write;

use gyrase_index::Index;
use gyrase_reads::{Header, Idxstats};
use lexopt::Parser;
use lexopt::prelude::*;

use super::input::{Input, sibling};
use crate::failure::Failure;
use crate::shared::{create, workers};

pub const USAGE: &str = "\
Usage: samtools idxstats [-@ INT] <in.bam>

Options:
  -@, --threads INT  number of additional threads
";

pub fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut threads = None;
	let mut path = None;
	while let Some(arg) = parser.next()? {
		match arg {
			Short('@') | Long("threads") => threads = Some(parser.value()?.parse()?),
			Long("help") => return Err(Failure::Help),
			Value(value) if path.is_none() => path = Some(value.string()?),
			Value(_) => return Err(Failure::usage("too many arguments")),
			_ => return Err(arg.unexpected().into()),
		}
	}
	let path = path.ok_or_else(|| Failure::usage("no input file"))?;
	let workers = workers(threads)?;
	let input = Input::open(&path, None, &workers)?;
	let header = input.header.clone();
	let text = match sibling(&path, &[".bai", ".csi"]) {
		Some(index_path) if input.is_bam() => {
			from_index(&header, &Index::read(File::open(index_path)?, &workers)?)
		}
		_ => {
			let mut stats = Idxstats::new(header.references().len());
			input.each(|record| Ok(stats.add(record)?))?;
			stats.render(&header)
		}
	};
	let mut out = create(None)?;
	out.write_all(text.as_bytes())?;
	out.flush()?;
	Ok(())
}

fn from_index(header: &Header, index: &Index) -> String {
	let mut text = String::new();
	for (position, reference) in header.references().iter().enumerate() {
		let (mapped, unmapped) = index
			.references()
			.get(position)
			.and_then(|entry| entry.metadata(index.geometry()))
			.map_or((0, 0), |metadata| (metadata.mapped, metadata.unmapped));
		let _ = writeln!(
			text,
			"{}\t{}\t{mapped}\t{unmapped}",
			reference.name(),
			reference.length()
		);
	}
	let _ = writeln!(text, "*\t0\t0\t{}", index.unplaced());
	text
}
