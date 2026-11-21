use std::path::Path;

use gyrase_genome::Reference;
use gyrase_variants::{Duplicates, Mismatch, Norm, NormOptions, Splitting};
use lexopt::Parser;
use lexopt::prelude::*;

use super::input::{Selection, each};
use super::output::{Collector, Common};
use crate::failure::Failure;

pub const USAGE: &str = "\
Usage: bcftools norm [options] <in.vcf.gz>|<in.bcf>

Options:
  -c, --check-ref e|w|x      reference mismatch: e error, w warn, x exclude [e]
  -d, --rm-dup TYPE          remove duplicates: snps, indels, both, all or exact
  -f, --fasta-ref FILE       reference FASTA, required to left-align indels
  -m, --multiallelics -TYPE  split multiallelic sites, TYPE is snps, indels, both or any [both]
  -r, --regions REG[,...]    restrict to regions, needs an index
  -R, --regions-file FILE    restrict to regions listed in FILE, needs an index
  -t, --targets REG[,...]    restrict to positions in regions, no index needed
  -T, --targets-file FILE    restrict to positions in regions listed in FILE
  -o, --output FILE          write output to FILE [standard output]
  -O, --output-type b|u|z|v  b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: VCF
  -w, --site-window INT      buffer size for sorting records, in sites [1000]
  -W, --write-index[=FMT]    write an index with the output, csi or tbi
      --threads INT          number of additional threads
";

fn splitting(text: &str) -> Result<Splitting, Failure> {
	match text {
		"-any" => Ok(Splitting::Any),
		"-both" | "-" => Ok(Splitting::Both),
		"-snps" => Ok(Splitting::Snps),
		"-indels" => Ok(Splitting::Indels),
		other if other.starts_with('+') => {
			Err(Failure::new("joining records with -m + is not supported"))
		}
		other => Err(Failure::usage(format!(
			"the argument \"{other}\" of -m is not recognised"
		))),
	}
}

fn duplicates(text: &str) -> Result<Duplicates, Failure> {
	match text {
		"exact" | "none" => Ok(Duplicates::Exact),
		"snps" => Ok(Duplicates::Snps),
		"indels" => Ok(Duplicates::Indels),
		"both" => Ok(Duplicates::Both),
		"any" | "all" => Ok(Duplicates::Any),
		other => Err(Failure::usage(format!(
			"the argument \"{other}\" of -d is not recognised"
		))),
	}
}

fn mismatch(text: &str) -> Result<Mismatch, Failure> {
	match text {
		"e" => Ok(Mismatch::Fail),
		"w" => Ok(Mismatch::Keep),
		"x" => Ok(Mismatch::Skip),
		other => Err(Failure::usage(format!(
			"the argument \"{other}\" of -c is not recognised"
		))),
	}
}

pub fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut common = Common::new();
	let mut options = NormOptions::default();
	let mut selection = Selection::default();
	let mut path = None;
	while let Some(arg) = parser.next()? {
		match arg {
			Short('O') | Long("output-type") => common.output_type(&mut parser)?,
			Short('o') | Long("output") => common.output(&mut parser)?,
			Short('W') | Long("write-index") => common.write_index(&mut parser)?,
			Long("threads") => common.thread_count(&mut parser)?,
			Short('f') | Long("fasta-ref") => {
				let fasta = parser.value()?.string()?;
				options.reference = Some(Reference::open(Path::new(&fasta))?);
			}
			Short('m') | Long("multiallelics") => {
				options.split = Some(splitting(&parser.value()?.string()?)?);
			}
			Short('d') | Long("rm-dup") => {
				options.duplicates = duplicates(&parser.value()?.string()?)?;
			}
			Short('c') | Long("check-ref") => {
				options.mismatch = mismatch(&parser.value()?.string()?)?;
			}
			Short('r') | Long("regions") => selection.regions(&parser.value()?.string()?),
			Short('R') | Long("regions-file") => {
				selection.regions_file(&parser.value()?.string()?)?;
			}
			Short('t') | Long("targets") => selection.targets(&parser.value()?.string()?),
			Short('T') | Long("targets-file") => {
				selection.targets_file(&parser.value()?.string()?)?;
			}
			Short('w') | Long("site-window") => options.window = parser.value()?.parse()?,
			Long("help") => return Err(Failure::Help),
			Value(value) if path.is_none() => path = Some(value.string()?),
			Value(_) => return Err(Failure::usage("too many arguments")),
			_ => return Err(arg.unexpected().into()),
		}
	}
	let workers = common.workers()?;
	let destination = &common.destination;
	let (source, regions) = selection.open(path.as_deref().unwrap_or("-"), &workers)?;
	options.format = destination.kind.format;
	options.regions = regions;
	let norm = Norm::new(source.header(), options);
	let mut collector = Collector::new();
	let mut writer = destination.writer(norm.header(), norm.format(), &workers)?;
	each!(source, reader => norm.run(reader, &mut writer, &workers, collector.mark())?);
	destination.finish(writer, &collector, norm.header(), &workers)
}
