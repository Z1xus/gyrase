use std::io::Write;
use std::path::Path;

use gyrase_genome::{Interval, Reference};
use lexopt::Parser;
use lexopt::prelude::*;

use crate::failure::Failure;
use crate::shared::Target;

pub(super) const USAGE: &str = "\
Usage: samtools faidx <file.fa> [region ...]

Options:
  -o, --output FILE  write FASTA to FILE
  -n, --length INT   FASTA line length, 0 for no wrapping [as in the input]
";

pub(super) fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut output = None;
	let mut width = None;
	let mut positional = Vec::new();
	while let Some(arg) = parser.next()? {
		match arg {
			Short('o') | Long("output") => output = Some(parser.value()?.string()?),
			Short('n') | Long("length") => width = Some(parser.value()?.parse::<usize>()?),
			Long("help") => return Err(Failure::Help),
			Value(value) => positional.push(value.string()?),
			_ => return Err(arg.unexpected().into()),
		}
	}
	let mut positional = positional.into_iter();
	let path = positional
		.next()
		.ok_or_else(|| Failure::usage("no input file"))?;
	let regions: Vec<String> = positional.collect();
	let target = Target::new(output.as_deref(), &[&path])?;
	let reference = Reference::open_indexed(Path::new(&path)).map_err(|error| {
		Failure::new(format!("could not build the index of \"{path}\": {error}"))
	})?;
	if regions.is_empty() {
		return Ok(());
	}
	let mut out = target.open()?;
	let mut bases = Vec::new();
	for region in &regions {
		let interval = Interval::parse(region, reference.contigs())
			.map_err(|_| Failure::new(format!("failed to fetch sequence in {region}")))?;
		reference.fetch_cased(&interval, &mut bases)?;
		writeln!(out, ">{region}")?;
		let line_length = match width {
			Some(0) => usize::MAX,
			Some(width) => width,
			None => reference
				.entries()
				.get(interval.contig.index())
				.and_then(|entry| usize::try_from(entry.line_bases).ok())
				.unwrap_or(1)
				.max(1),
		};
		for line in bases.chunks(line_length) {
			out.write_all(line)?;
			out.write_all(b"\n")?;
		}
	}
	target.commit(out)
}
