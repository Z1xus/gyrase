use gyrase_variants::{Collapse, Merge, MergeOptions, Reader};
use lexopt::Parser;
use lexopt::prelude::*;

use super::concat::list;
use super::output::{Collector, Common};
use crate::failure::Failure;
use crate::shared::open;

pub const USAGE: &str = "\
Usage: bcftools merge [options] <a.vcf.gz> <b.vcf.gz> [...]

Options:
  -l, --file-list FILE       read the input file names from FILE
  -m, --merge TYPE           allow multiallelic records for none, snps, indels, both, all or exact [both]
  -o, --output FILE          write output to FILE [standard output]
  -O, --output-type b|u|z|v  b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: VCF
  -W, --write-index[=FMT]    write an index with the output, csi or tbi
      --threads INT          number of additional threads
";

fn collapse(text: &str) -> Result<Collapse, Failure> {
	match text {
		"none" => Ok(Collapse::None),
		"snps" => Ok(Collapse::Snps),
		"indels" => Ok(Collapse::Indels),
		"both" => Ok(Collapse::Both),
		"all" | "any" => Ok(Collapse::All),
		"exact" => Ok(Collapse::Exact),
		other => Err(Failure::usage(format!(
			"the argument \"{other}\" of -m is not recognised"
		))),
	}
}

pub fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut common = Common::new();
	let mut paths = Vec::new();
	let mut merging = Collapse::default();
	while let Some(arg) = parser.next()? {
		match arg {
			Short('O') | Long("output-type") => common.output_type(&mut parser)?,
			Short('o') | Long("output") => common.output(&mut parser)?,
			Short('W') | Long("write-index") => common.write_index(&mut parser)?,
			Long("threads") => common.thread_count(&mut parser)?,
			Short('m') | Long("merge") => merging = collapse(&parser.value()?.string()?)?,
			Short('l') | Long("file-list") => paths.extend(list(&parser.value()?.string()?)?),
			Long("help") => return Err(Failure::Help),
			Value(value) => paths.push(value.string()?),
			_ => return Err(arg.unexpected().into()),
		}
	}
	if paths.is_empty() {
		return Err(Failure::usage("no input files"));
	}
	let workers = common.workers()?;
	let destination = &common.destination;
	let readers = paths
		.iter()
		.map(|path| Ok(Reader::open(open(path)?, &workers)?))
		.collect::<Result<Vec<_>, Failure>>()?;
	let headers: Vec<_> = readers.iter().map(Reader::header).collect();
	let merge = Merge::new(
		&headers,
		MergeOptions {
			format: destination.kind.format,
			collapse: merging,
		},
	)?;
	let mut collector = Collector::new();
	let mut writer = destination.writer(merge.header(), merge.format(), &workers)?;
	merge.run(readers, &mut writer, collector.mark())?;
	destination.finish(writer, &collector, merge.header(), &workers)
}
