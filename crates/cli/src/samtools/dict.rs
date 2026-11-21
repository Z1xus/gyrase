use std::fmt::Write as _;
use std::io::Write;
use std::path::Path;

use gyrase_genome::{Interval, Reference};
use lexopt::Parser;
use lexopt::prelude::*;

use crate::failure::Failure;
use crate::shared::Target;

pub(super) const USAGE: &str = "\
Usage: samtools dict [options] <file.fa>

Options:
  -a, --assembly STR  assembly
  -A, --alias         add an AN tag by adding or removing 'chr'
  -H, --no-header     do not print the @HD line
  -o, --output FILE   write the dictionary to FILE [stdout]
  -s, --species STR   species
  -u, --uri STR       URI [file:///abs/path/to/file.fa]
";

pub(super) fn run(mut parser: Parser) -> Result<(), Failure> {
	let mut assembly = None;
	let mut species = None;
	let mut uri = None;
	let mut output = None;
	let mut alias = false;
	let mut header = true;
	let mut path = None;
	while let Some(arg) = parser.next()? {
		match arg {
			Short('a') | Long("assembly") => assembly = Some(parser.value()?.string()?),
			Short('s') | Long("species") => species = Some(parser.value()?.string()?),
			Short('u') | Long("uri") => uri = Some(parser.value()?.string()?),
			Short('o') | Long("output") => output = Some(parser.value()?.string()?),
			Short('A') | Long("alias" | "alternative-name") => alias = true,
			Short('H') | Long("no-header") => header = false,
			Long("help") => return Err(Failure::Help),
			Value(value) if path.is_none() => path = Some(value.string()?),
			Value(_) => return Err(Failure::usage("too many arguments")),
			_ => return Err(arg.unexpected().into()),
		}
	}
	let path = path.ok_or_else(|| Failure::usage("no input file"))?;
	let target = Target::new(output.as_deref(), &[&path])?;
	let reference = Reference::open(Path::new(&path))?;
	let uri = match uri {
		Some(uri) => uri,
		None => file_uri(Path::new(&path))?,
	};
	let mut text = String::new();
	if header {
		text.push_str("@HD\tVN:1.0\tSO:unsorted\n");
	}
	for (id, contig) in reference.contigs().iter() {
		let checksum = reference.checksum(
			&Interval::whole(id, reference.contigs())
				.ok_or_else(|| Failure::new("invalid contig"))?,
		)?;
		let _ = write!(text, "@SQ\tSN:{}\tLN:{}\tM5:", contig.name, contig.length);
		for byte in checksum {
			let _ = write!(text, "{byte:02x}");
		}
		if alias {
			match contig.name.strip_prefix("chr") {
				Some(short) => {
					let _ = write!(text, "\tAN:{short}");
				}
				None => {
					let _ = write!(text, "\tAN:chr{}", contig.name);
				}
			}
		}
		let _ = write!(text, "\tUR:{uri}");
		if let Some(assembly) = &assembly {
			let _ = write!(text, "\tAS:{assembly}");
		}
		if let Some(species) = &species {
			let _ = write!(text, "\tSP:{species}");
		}
		text.push('\n');
	}
	let mut out = target.open()?;
	out.write_all(text.as_bytes())?;
	target.commit(out)
}

fn file_uri(path: &Path) -> Result<String, Failure> {
	let canonical = std::fs::canonicalize(path)?
		.to_string_lossy()
		.replace('\\', "/");
	let plain = canonical.strip_prefix("//?/").unwrap_or(&canonical);
	Ok(if plain.starts_with('/') {
		format!("file://{plain}")
	} else {
		format!("file:///{plain}")
	})
}
