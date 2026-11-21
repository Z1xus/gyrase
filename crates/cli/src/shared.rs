use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::num::NonZeroUsize;

use gyrase_io::Workers;

use crate::failure::{Failure, file};

pub type Output = BufWriter<Box<dyn Write + Send>>;

const OUTPUT_BUFFER: usize = 1 << 20;

pub fn workers(extra: Option<usize>) -> Result<Workers, Failure> {
	match extra {
		None => Ok(Workers::available()?),
		Some(extra) => {
			let threads = NonZeroUsize::new(extra.saturating_add(1))
				.ok_or_else(|| Failure::new("invalid thread count"))?;
			Ok(Workers::new(threads)?)
		}
	}
}

pub fn open(path: &str) -> Result<Box<dyn Read + Send>, Failure> {
	if path == "-" {
		return Ok(Box::new(std::io::stdin()));
	}
	Ok(Box::new(
		File::open(path).map_err(|error| file(path, "reading", &error))?,
	))
}

pub fn create(path: Option<&str>) -> Result<Output, Failure> {
	let sink: Box<dyn Write + Send> = match path {
		None | Some("-") => Box::new(std::io::stdout()),
		Some(path) => Box::new(File::create(path).map_err(|error| file(path, "writing", &error))?),
	};
	Ok(BufWriter::with_capacity(OUTPUT_BUFFER, sink))
}
