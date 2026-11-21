use gyrase_variants::{Condition, Counts, Fidelity, Format, SampleChoice, View, ViewOptions};
use lexopt::Parser;
use lexopt::prelude::*;

use super::input::{Selection, each, set_condition};
use super::output::{Collector, Common};
use crate::failure::{Failure, file};

pub(super) const USAGE: &str = "\
Usage: bcftools view [options] <in.vcf.gz>|<in.bcf>

Options:
  -e, --exclude EXPR         exclude sites for which EXPR is true
  -i, --include EXPR         include only sites for which EXPR is true
  -h, --header-only          print the header only
  -H, --no-header            suppress the header
  -I, --no-update            do not recalculate INFO/AC and INFO/AN after subsetting samples
  -o, --output FILE          write output to FILE [standard output]
  -O, --output-type b|u|z|v  b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: VCF
  -r, --regions REG[,...]    restrict to regions, needs an index
  -R, --regions-file FILE    restrict to regions listed in FILE, needs an index
  -s, --samples [^]LIST      keep or, with ^, drop the comma-separated samples
  -S, --samples-file FILE    keep the samples listed in FILE
  -t, --targets REG[,...]    restrict to positions in regions, no index needed
  -T, --targets-file FILE    restrict to positions in regions listed in FILE
  -W, --write-index[=FMT]    write an index with the output, csi or tbi
      --threads INT          number of additional threads
";

fn samples(list: &str) -> SampleChoice {
	match list.strip_prefix('^') {
		Some(rest) => SampleChoice::Drop(rest.split(',').map(str::to_owned).collect()),
		None => SampleChoice::Keep(list.split(',').map(str::to_owned).collect()),
	}
}

fn sample_file(argument: &str) -> Result<SampleChoice, Failure> {
	let (path, drop) = match argument.strip_prefix('^') {
		Some(path) => (path, true),
		None => (argument, false),
	};
	let text = std::fs::read_to_string(path).map_err(|error| file(path, "reading", &error))?;
	let names = text
		.lines()
		.map(|line| {
			line.split_whitespace()
				.next()
				.unwrap_or_default()
				.to_owned()
		})
		.filter(|name| !name.is_empty())
		.collect();
	Ok(if drop {
		SampleChoice::Drop(names)
	} else {
		SampleChoice::Keep(names)
	})
}

pub(super) fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut common = Common::new();
	let mut selection = Selection::default();
	let mut options = ViewOptions::default();
	let mut condition = None;
	let mut rows = Rows::Both;
	let mut path = None;
	while let Some(arg) = parser.next()? {
		match arg {
			Short('O') | Long("output-type") => common.output_type(&mut parser)?,
			Short('o') | Long("output") => common.output(&mut parser)?,
			Short('W') | Long("write-index") => common.write_index(&mut parser)?,
			Long("threads") => common.thread_count(&mut parser)?,
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
					&mut condition,
					Condition::Include(parser.value()?.string()?),
				)?;
			}
			Short('e') | Long("exclude") => {
				set_condition(
					&mut condition,
					Condition::Exclude(parser.value()?.string()?),
				)?;
			}
			Short('s') | Long("samples") => options.samples = samples(&parser.value()?.string()?),
			Short('S') | Long("samples-file") => {
				options.samples = sample_file(&parser.value()?.string()?)?;
			}
			Short('I') | Long("no-update") => options.counts = Counts::Keep,
			Short('h') | Long("header-only") => rows = rows.choose(Rows::HeaderOnly)?,
			Short('H') | Long("no-header") => rows = rows.choose(Rows::RecordsOnly)?,
			Long("help") => return Err(Failure::Help),
			Value(value) if path.is_none() => path = Some(value.string()?),
			Value(_) => return Err(Failure::usage("too many arguments")),
			_ => return Err(arg.unexpected().into()),
		}
	}
	let workers = common.workers()?;
	let destination = &common.destination;
	if rows == Rows::RecordsOnly && destination.kind.format() == Format::Bcf {
		return Err(Failure::new(
			"BCF output requires a header, -H cannot be used",
		));
	}
	let path = path.unwrap_or_else(|| String::from("-"));
	let target = common.target(&[&path])?;
	let (source, regions) = selection.open(&path, &workers)?;
	options.condition = condition;
	options.regions = regions;
	options.format = destination.kind.format();
	options.fidelity = Fidelity::Canonical;
	let view = View::new(source.header(), options)?;
	let mut collector = Collector::new();
	let mut writer = match rows {
		Rows::RecordsOnly => destination.bare(&workers, &target)?,
		Rows::Both | Rows::HeaderOnly => {
			destination.writer(view.header(), view.format(), &workers, &target)?
		}
	};
	if rows != Rows::HeaderOnly {
		each!(source, reader => view.run(reader, &mut writer, &workers, collector.mark())?);
	}
	destination.finish(writer, &collector, view.header(), &workers, &target)
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Rows {
	Both,
	HeaderOnly,
	RecordsOnly,
}

impl Rows {
	fn choose(self, next: Self) -> Result<Self, Failure> {
		match self {
			Self::Both => Ok(next),
			_ if self == next => Ok(next),
			_ => Err(Failure::usage("-h and -H cannot be combined")),
		}
	}
}
