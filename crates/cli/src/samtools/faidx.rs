use std::io::Write;
use std::path::Path;

use gyrase_genome::{Interval, Reference};
use lexopt::Parser;
use lexopt::prelude::*;

use crate::failure::Failure;
use crate::shared::create;

pub const USAGE: &str = "\
Usage: samtools faidx <file.fa> [region ...]

Options:
  -o, --output FILE  write FASTA to FILE
  -n, --length INT   FASTA line length [as in the input]
";

pub fn run(mut parser: Parser) -> Result<(), Failure> {
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
	let reference = Reference::open(Path::new(&path))?;
	if regions.is_empty() {
		return Ok(());
	}
	let mut out = create(output.as_deref())?;
	let mut bases = Vec::new();
	for region in &regions {
		let interval = Interval::parse(region, reference.contigs())
			.map_err(|_| Failure::new(format!("failed to fetch sequence in {region}")))?;
		reference.fetch_cased(&interval, &mut bases)?;
		writeln!(out, ">{region}")?;
		let line_length = width
			.or_else(|| {
				let entry = reference.entries().get(interval.contig.index())?;
				usize::try_from(entry.line_bases).ok()
			})
			.unwrap_or(1)
			.max(1);
		for line in bases.chunks(line_length) {
			out.write_all(line)?;
			out.write_all(b"\n")?;
		}
	}
	out.flush()?;
	Ok(())
}
