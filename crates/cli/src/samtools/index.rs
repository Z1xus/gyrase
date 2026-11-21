use std::fs::File;
use std::io::Read;

use gyrase_index::{Format, Geometry};
use gyrase_reads::cram::index_cram;
use gyrase_reads::{BamReader, csi_geometry, index_bam};
use lexopt::Parser;
use lexopt::prelude::*;

use crate::failure::{Failure, file};
use crate::shared::{Target, open, workers};

pub(super) const USAGE: &str = "\
Usage: samtools index [-bc] [-m INT] [-@ INT] <in.bam|in.cram> [out.index]

Options:
  -b, --bai            BAI index [default]
  -c, --csi            CSI index
  -m, --min-shift INT  CSI interval size is 2^INT [14]
  -o, --output FILE    write the index to FILE
  -@, --threads INT    number of additional threads
";

const BAI_MIN_SHIFT: u32 = 14;

pub(super) fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut csi = false;
	let mut min_shift = BAI_MIN_SHIFT;
	let mut output: Option<String> = None;
	let mut threads = None;
	let mut positional = Vec::new();
	while let Some(arg) = parser.next()? {
		match arg {
			Short('b') | Long("bai") => csi = false,
			Short('c') | Long("csi") => csi = true,
			Short('m') | Long("min-shift") => {
				min_shift = parser.value()?.parse()?;
				csi = true;
			}
			Short('o') | Long("output") => output = Some(parser.value()?.string()?),
			Short('@') | Long("threads") => threads = Some(parser.value()?.parse()?),
			Long("help") => return Err(Failure::Help),
			Value(value) => positional.push(value.string()?),
			_ => return Err(arg.unexpected().into()),
		}
	}
	let mut positional = positional.into_iter();
	let path = positional
		.next()
		.ok_or_else(|| Failure::usage("no input file"))?;
	if output.is_none() {
		output = positional.next();
	}
	if positional.next().is_some() {
		return Err(Failure::usage("too many arguments"));
	}
	let workers = workers(threads)?;
	if is_cram(&path)? {
		let output = output.unwrap_or_else(|| format!("{path}.crai"));
		let target = Target::new(Some(&output), &[&path])?;
		let source = File::open(&path).map_err(|error| file(&path, "reading", &error))?;
		let index = index_cram(source, |_| Ok::<_, gyrase_reads::cram::Error>(()))?;
		return target.commit(index.write(target.open()?, &workers)?);
	}
	let output = output.unwrap_or_else(|| format!("{path}.{}", if csi { "csi" } else { "bai" }));
	let target = Target::new(Some(&output), &[&path])?;
	let mut reader = BamReader::open(open(&path)?, &workers)?;
	let geometry = if csi {
		csi_geometry(reader.header(), min_shift)?
	} else {
		Geometry::BAI
	};
	let index = index_bam(&mut reader, geometry, |_| Ok::<_, gyrase_reads::Error>(()))?;
	let format = if csi { Format::Csi } else { Format::Bai };
	target.commit(index.write(format, target.open()?, &workers)?)
}

fn is_cram(path: &str) -> Result<bool, Failure> {
	if path == "-" {
		return Ok(false);
	}
	let mut magic = [0; 4];
	let mut source = File::open(path).map_err(|error| file(path, "reading", &error))?;
	Ok(source.read_exact(&mut magic).is_ok() && &magic == b"CRAM")
}
