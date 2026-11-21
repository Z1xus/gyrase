use std::io::Write;

use lexopt::Parser;
use lexopt::prelude::*;

use super::input::{Input, Need, reference};
use crate::failure::Failure;
use crate::shared::{stdout, workers};

pub(super) const USAGE: &str = "\
Usage: samtools flagstat [-@ INT] [-T FILE] <in.bam>|<in.sam>|<in.cram>

Options:
  -T, --reference FILE  reference FASTA, needed for CRAM
  -@, --threads INT     number of additional threads
";

pub(super) fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut threads = None;
	let mut reference_path: Option<String> = None;
	let mut path = None;
	while let Some(arg) = parser.next()? {
		match arg {
			Short('@') | Long("threads") => threads = Some(parser.value()?.parse()?),
			Short('T') | Long("reference") => reference_path = Some(parser.value()?.string()?),
			Long("help") => return Err(Failure::Help),
			Value(value) if path.is_none() => path = Some(value.string()?),
			Value(_) => return Err(Failure::usage("too many arguments")),
			_ => return Err(arg.unexpected().into()),
		}
	}
	let path = path.ok_or_else(|| Failure::usage("no input file"))?;
	let workers = workers(threads)?;
	let reference = reference_path.as_deref().map(reference).transpose()?;
	let stats =
		Input::open(&path, reference.as_ref(), Need::Positions, &workers)?.flagstat(&workers)?;
	let mut out = stdout();
	out.write_all(stats.render().as_bytes())?;
	out.flush()?;
	Ok(())
}
