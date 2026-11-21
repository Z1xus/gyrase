use gyrase_variants::{Concat, ConcatOptions, Fidelity, Reader};
use lexopt::Parser;
use lexopt::prelude::*;

use super::output::{Collector, Common};
use crate::failure::{Failure, file};
use crate::shared::open;

pub(super) const USAGE: &str = "\
Usage: bcftools concat [options] <a.vcf.gz> [<b.vcf.gz> ...]

Options:
  -f, --file-list FILE       read the input file names from FILE
  -o, --output FILE          write output to FILE [standard output]
  -O, --output-type b|u|z|v  b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: VCF
  -W, --write-index[=FMT]    write an index with the output, csi or tbi
      --threads INT          number of additional threads
";

pub(super) fn list(path: &str) -> Result<Vec<String>, Failure> {
	let text = std::fs::read_to_string(path).map_err(|error| file(path, "reading", &error))?;
	Ok(text
		.lines()
		.map(str::trim)
		.filter(|line| !line.is_empty())
		.map(str::to_owned)
		.collect())
}

pub(super) fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut common = Common::new();
	let mut paths = Vec::new();
	while let Some(arg) = parser.next()? {
		match arg {
			Short('O') | Long("output-type") => common.output_type(&mut parser)?,
			Short('o') | Long("output") => common.output(&mut parser)?,
			Short('W') | Long("write-index") => common.write_index(&mut parser)?,
			Long("threads") => common.thread_count(&mut parser)?,
			Short('f') | Long("file-list") => paths.extend(list(&parser.value()?.string()?)?),
			Long("help") => return Err(Failure::Help),
			Value(value) => paths.push(value.string()?),
			_ => return Err(arg.unexpected().into()),
		}
	}
	if paths.is_empty() {
		return Err(Failure::usage("no input files"));
	}
	if paths.iter().any(|path| path == "-") {
		return Err(Failure::new("standard input cannot be an input"));
	}
	let workers = common.workers()?;
	let destination = &common.destination;
	let headers = paths
		.iter()
		.map(|path| Ok(Reader::open(open(path)?, &workers)?.header().clone()))
		.collect::<Result<Vec<_>, Failure>>()?;
	let concat = Concat::new(
		&headers.iter().collect::<Vec<_>>(),
		ConcatOptions {
			fidelity: Fidelity::Canonical,
			format: destination.kind.format(),
		},
	)?;
	let inputs: Vec<&str> = paths.iter().map(String::as_str).collect();
	let target = common.target(&inputs)?;
	let mut collector = Collector::new();
	let mut writer = destination.writer(concat.header(), concat.format(), &workers, &target)?;
	let readers = paths
		.iter()
		.map(|path| Ok::<_, Failure>(Reader::open(open(path)?, &workers)?));
	concat.run(readers, &mut writer, &workers, collector.mark())?;
	destination.finish(writer, &collector, concat.header(), &workers, &target)
}
