use std::io::Write;

use lexopt::Parser;
use lexopt::prelude::*;

use super::input::Input;
use crate::failure::Failure;
use crate::shared::{create, workers};

pub const USAGE: &str = "\
Usage: samtools flagstat [-@ INT] <in.bam>|<in.sam>|<in.cram>

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
	let stats = Input::open(&path, None, &workers)?.flagstat(&workers)?;
	let mut out = create(None)?;
	out.write_all(stats.render().as_bytes())?;
	out.flush()?;
	Ok(())
}
