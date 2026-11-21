use std::io::Write;

use gyrase_variants::{Condition, Query, QueryOptions};
use lexopt::Parser;
use lexopt::prelude::*;

use super::input::{Selection, each, set_condition};
use crate::failure::Failure;
use crate::shared::{Target, workers};

pub(super) const USAGE: &str = "\
Usage: bcftools query [options] <in.vcf.gz>|<in.bcf>

Options:
  -e, --exclude EXPR       exclude sites for which EXPR is true
  -f, --format STRING      output format, for example '%CHROM\\t%POS\\n'
  -i, --include EXPR       include only sites for which EXPR is true
  -l, --list-samples       list the sample names and exit
  -o, --output FILE        write output to FILE [standard output]
  -r, --regions REG[,...]  restrict to regions, needs an index
  -R, --regions-file FILE  restrict to regions listed in FILE, needs an index
  -s, --samples LIST       keep the comma-separated samples
  -t, --targets REG[,...]  restrict to positions in regions, no index needed
  -T, --targets-file FILE  restrict to positions in regions listed in FILE
      --threads INT        number of additional threads
";

pub(super) fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut selection = Selection::default();
	let mut options = QueryOptions::default();
	let mut format = None;
	let mut output = None;
	let mut threads = None;
	let mut list = false;
	let mut path = None;
	while let Some(arg) = parser.next()? {
		match arg {
			Short('f') | Long("format") => format = Some(parser.value()?.string()?),
			Short('o') | Long("output") => output = Some(parser.value()?.string()?),
			Long("threads") => threads = Some(parser.value()?.parse()?),
			Short('l') | Long("list-samples") => list = true,
			Short('r') | Long("regions") => selection.regions(&parser.value()?.string()?),
			Short('R') | Long("regions-file") => {
				selection.regions_file(&parser.value()?.string()?)?;
			}
			Short('t') | Long("targets") => selection.targets(&parser.value()?.string()?),
			Short('T') | Long("targets-file") => {
				selection.targets_file(&parser.value()?.string()?)?;
			}
			Short('i') | Long("include") => {
				set_condition(
					&mut options.condition,
					Condition::Include(parser.value()?.string()?),
				)?;
			}
			Short('e') | Long("exclude") => {
				set_condition(
					&mut options.condition,
					Condition::Exclude(parser.value()?.string()?),
				)?;
			}
			Short('s') | Long("samples") => {
				let list = parser.value()?.string()?;
				options.samples = Some(list.split(',').map(str::to_owned).collect());
			}
			Long("help") => return Err(Failure::Help),
			Value(value) if path.is_none() => path = Some(value.string()?),
			Value(_) => return Err(Failure::usage("too many arguments")),
			_ => return Err(arg.unexpected().into()),
		}
	}
	let workers = workers(threads)?;
	let path = path.unwrap_or_else(|| String::from("-"));
	let target = Target::new(output.as_deref(), &[&path])?;
	let (source, regions) = selection.open(&path, &workers)?;
	if list {
		let mut out = target.open()?;
		for name in source.header().samples() {
			writeln!(out, "{name}")?;
		}
		return target.commit(out);
	}
	let format = format.ok_or_else(|| Failure::usage("the -f option is required"))?;
	options.regions = regions;
	let query = Query::new(source.header(), &format, options)?;
	let mut out = target.open()?;
	each!(source, reader => query.run(reader, &mut out, &workers)?);
	target.commit(out)
}
