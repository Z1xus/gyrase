use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};

use gyrase_io::Workers;
use gyrase_reads::BamReader;
use gyrase_reads::cram::CramReader;
use lexopt::Parser;
use lexopt::prelude::*;

use crate::failure::Failure;
use crate::shared::workers;

pub const USAGE: &str = "\
Usage: samtools quickcheck [options] <input> [...]

Options:
  -v  print the names of failed files to stdout
  -q  suppress warning messages
  -u  unmapped input, do not require references in the header
";

const BGZF_EOF: [u8; 28] = [
	0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0,
	0, 0, 0, 0, 0, 0, 0,
];
const CRAM_EOF: [u8; 38] = [
	0x0f, 0, 0, 0, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0, 0, 0, 0, 0x01, 0, 0x05,
	0xbd, 0xd9, 0x4f, 0, 0x01, 0, 0x06, 0x06, 0x01, 0, 0x01, 0, 0x01, 0, 0xee, 0x63, 0x01, 0x4b,
];

fn has_tail(file: &mut File, tail: &[u8]) -> Result<bool, Failure> {
	let length = i64::try_from(tail.len()).map_err(|_| Failure::new("invalid length"))?;
	let mut found = vec![0; tail.len()];
	Ok(file.seek(SeekFrom::End(-length)).is_ok()
		&& file.read_exact(&mut found).is_ok()
		&& found == tail)
}

fn problem(path: &str, unmapped: bool, workers: &Workers) -> Result<Option<&'static str>, Failure> {
	let Ok(mut file) = File::open(path) else {
		return Ok(Some("could not be opened for reading."));
	};
	let mut head = [0u8; 4];
	let read = file.read(&mut head)?;
	let head = head.get(..read).unwrap_or_default();
	if head.starts_with(&[0x1f, 0x8b]) {
		if !has_tail(&mut file, &BGZF_EOF)? {
			return Ok(Some("was missing EOF block when one should be present."));
		}
		let Ok(header) = BamReader::read_header(&file, workers) else {
			return Ok(Some("could not read the header."));
		};
		if !unmapped && header.references().is_empty() {
			return Ok(Some("had no targets in header."));
		}
	} else if head.starts_with(b"CRAM") {
		file.rewind()?;
		let Ok(reader) = CramReader::open(Box::new(file.try_clone()?), None, path, workers) else {
			return Ok(Some("could not read the header."));
		};
		if !has_tail(&mut file, &CRAM_EOF)? {
			return Ok(Some("was missing EOF block when one should be present."));
		}
		if !unmapped && reader.header().references().is_empty() {
			return Ok(Some("had no targets in header."));
		}
	} else if !head.starts_with(b"@") {
		return Ok(Some("was not identified as sequence data."));
	}
	Ok(None)
}

pub fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut verbose = false;
	let mut quiet = false;
	let mut unmapped = false;
	let mut paths = Vec::new();
	while let Some(arg) = parser.next()? {
		match arg {
			Short('v') => verbose = true,
			Short('q') => quiet = true,
			Short('u') => unmapped = true,
			Long("help") => return Err(Failure::Help),
			Value(value) => paths.push(value.string()?),
			_ => return Err(arg.unexpected().into()),
		}
	}
	if paths.is_empty() {
		return Err(Failure::usage("no input files"));
	}
	let workers = workers(Some(0))?;
	let mut failed = false;
	for path in &paths {
		let Some(message) = problem(path, unmapped, &workers)? else {
			continue;
		};
		failed = true;
		if !quiet {
			let _ = writeln!(std::io::stderr(), "{path} {message}");
		}
		if verbose {
			writeln!(std::io::stdout(), "{path}")?;
		}
	}
	if failed {
		return Err(Failure::Failed);
	}
	Ok(())
}
