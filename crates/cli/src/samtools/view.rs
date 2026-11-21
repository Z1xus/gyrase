use std::io::{IsTerminal, Write};

use gyrase_genome::Interval;
use gyrase_io::Level;
use gyrase_reads::{Header, RecordRef};
use gyrase_variants::Regions;
use lexopt::Parser;
use lexopt::prelude::*;

use super::input::{Input, Need, reference};
use super::output::{Format, Output};
use crate::failure::{Failure, file};
use crate::shared::{Output as Sink, Target, workers};

pub(super) const USAGE: &str = "\
Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]

Options:
  -b, --bam               output BAM
  -C, --cram              output CRAM (requires -T)
  -O, --output-fmt FMT    output format: sam, bam or cram
  -u, --uncompressed      uncompressed BAM output
  -1, --fast              fast compression
  -h, --with-header       include the header in SAM output
  -H, --header-only       print the header only
  -c, --count             print only the number of records
  -M, --multi-region      visit overlapping regions once, in file order
  -o, --output FILE       write output to FILE [stdout]
  -T, --reference FILE    reference FASTA
  -q, --min-MQ INT        only records with mapping quality >= INT
  -f, --require-flags INT only records with all of the flags in INT set
  -F, --excl-flags INT    only records with none of the flags in INT set
  -L, --target-file FILE  only records overlapping the regions in a BED file
  -@, --threads INT       number of additional threads
";

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Content {
	Records,
	WithHeader,
	HeaderOnly,
	Count,
}

#[derive(Default)]
struct Filter {
	min_mapq: u8,
	required: u16,
	excluded: u16,
	targets: Option<Regions>,
}

impl Filter {
	fn accepts(&self, record: RecordRef<'_>) -> bool {
		let flags = record.flags().bits();
		if flags & self.required != self.required || flags & self.excluded != 0 {
			return false;
		}
		if record.mapq() < self.min_mapq {
			return false;
		}
		let Some(targets) = &self.targets else {
			return true;
		};
		let (Some(contig), Some(start), Some(end)) =
			(record.contig(), record.position(), record.reference_end())
		else {
			return false;
		};
		targets.overlaps(contig.get(), start, u32::try_from(end).unwrap_or(u32::MAX))
	}
}

pub(super) fn parse_flags(text: &str) -> Result<u16, Failure> {
	let parsed = match text.strip_prefix("0x").or_else(|| text.strip_prefix("0X")) {
		Some(hex) => u16::from_str_radix(hex, 16),
		None => text.parse(),
	};
	parsed.map_err(|_| Failure::usage(format!("invalid flag value \"{text}\"")))
}

fn load_targets(path: &str, header: &Header) -> Result<Regions, Failure> {
	let text = std::fs::read_to_string(path).map_err(|error| file(path, "reading", &error))?;
	let contigs = header.contigs()?;
	let mut intervals = Vec::new();
	for line in text.lines() {
		if line.starts_with(['#']) || line.starts_with("track") || line.starts_with("browser") {
			continue;
		}
		let mut fields = line.split_whitespace();
		let (Some(name), Some(start), Some(end)) = (fields.next(), fields.next(), fields.next())
		else {
			continue;
		};
		let (Some(contig), Ok(start), Ok(end)) =
			(contigs.id(name), start.parse::<u64>(), end.parse::<u64>())
		else {
			continue;
		};
		intervals.push(Interval { contig, start, end });
	}
	Ok(Regions::new(&intervals))
}

fn tally(input: Input, filter: &Filter, mut sink: Sink, target: &Target) -> Result<(), Failure> {
	let mut total = 0u64;
	input.each(|record| {
		total += u64::from(filter.accepts(record));
		Ok(())
	})?;
	writeln!(sink, "{total}")?;
	target.commit(sink)
}

pub(super) fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut format = None;
	let mut level = None;
	let mut content = Content::Records;
	let mut multi = false;
	let mut output = None;
	let mut threads = None;
	let mut reference_path: Option<String> = None;
	let mut target_file: Option<String> = None;
	let mut filter = Filter::default();
	let mut positional = Vec::new();
	while let Some(arg) = parser.next()? {
		match arg {
			Short('b') | Long("bam") => format = Some(Format::Bam),
			Short('C') | Long("cram") => format = Some(Format::Cram),
			Short('O') | Long("output-fmt") => {
				let name = parser.value()?.string()?;
				let name = name.split(',').next().unwrap_or_default();
				format =
					Some(Format::named(name).ok_or_else(|| {
						Failure::usage(format!("unknown output format \"{name}\""))
					})?);
			}
			Short('u') | Long("uncompressed") => {
				level = Some(Level::STORED);
			}
			Short('1') | Long("fast") => level = Some(Level::FASTEST),
			Short('h') | Long("with-header") => content = content.max(Content::WithHeader),
			Short('H') | Long("header-only") => content = content.max(Content::HeaderOnly),
			Short('c') | Long("count") => content = Content::Count,
			Short('M') | Long("multi-region") => multi = true,
			Short('o') | Long("output") => output = Some(parser.value()?.string()?),
			Short('@') | Long("threads") => threads = Some(parser.value()?.parse()?),
			Short('T') | Long("reference") => reference_path = Some(parser.value()?.string()?),
			Short('q') | Long("min-MQ") => filter.min_mapq = parser.value()?.parse()?,
			Short('f') | Long("require-flags") => {
				filter.required = parse_flags(&parser.value()?.string()?)?;
			}
			Short('F') | Long("excl-flags") => {
				filter.excluded = parse_flags(&parser.value()?.string()?)?;
			}
			Short('L') | Long("target-file") => target_file = Some(parser.value()?.string()?),
			Short('S') | Long("no-PG") => {}
			Long("help") => return Err(Failure::Help),
			Value(value) => positional.push(value.string()?),
			_ => return Err(arg.unexpected().into()),
		}
	}
	let mut positional = positional.into_iter();
	let path = match positional.next() {
		Some(path) => path,
		None if std::io::stdin().is_terminal() => return Err(Failure::usage("no input file")),
		None => String::from("-"),
	};
	let regions: Vec<String> = positional.collect();
	let format = format
		.or_else(|| output.as_deref().and_then(Format::of_path))
		.unwrap_or(if level.is_some() {
			Format::Bam
		} else {
			Format::Sam
		});
	let mut inputs = vec![path.as_str()];
	inputs.extend(reference_path.as_deref());
	inputs.extend(target_file.as_deref());
	let target = Target::new(output.as_deref(), &inputs)?;
	let workers = workers(threads)?;
	let reference = reference_path.as_deref().map(reference).transpose()?;
	let input = if regions.is_empty() {
		Input::open(&path, reference.as_ref(), Need::Sequences, &workers)?
	} else {
		Input::query(&path, &regions, multi, reference.as_ref(), &workers)?
	};
	if let Some(targets) = &target_file {
		filter.targets = Some(load_targets(targets, &input.header)?);
	}
	let sink = target.open()?;
	if content == Content::Count {
		return tally(input, &filter, sink, &target);
	}
	let mut writer = Output::new(
		sink,
		format,
		&input.header,
		content != Content::Records,
		level.unwrap_or(Level::DEFAULT),
		reference.as_ref(),
		&workers,
	)?;
	if content != Content::HeaderOnly {
		input.each(|record| {
			if filter.accepts(record) {
				writer.write(record)?;
			}
			Ok(())
		})?;
	}
	target.commit(writer.finish()?)
}
