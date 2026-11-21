use std::cell::Cell;
use std::fs::{self, File};
use std::io::{BufWriter, Read, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use gyrase_io::{Partial, Workers};

use crate::failure::{Failure, file};

pub(crate) type Output = BufWriter<Box<dyn Write + Send>>;

const OUTPUT_BUFFER: usize = 1 << 20;

pub(crate) fn workers(extra: Option<usize>) -> Result<Workers, Failure> {
	match extra {
		None => Ok(Workers::available()?),
		Some(extra) => {
			let threads = NonZeroUsize::new(extra.saturating_add(1))
				.ok_or_else(|| Failure::new("invalid thread count"))?;
			Ok(Workers::new(threads)?)
		}
	}
}

pub(crate) fn open(path: &str) -> Result<Box<dyn Read + Send>, Failure> {
	if path == "-" {
		return Ok(Box::new(std::io::stdin()));
	}
	Ok(Box::new(
		File::open(path).map_err(|error| file(path, "reading", &error))?,
	))
}

pub(crate) fn stdout() -> Output {
	BufWriter::with_capacity(OUTPUT_BUFFER, Box::new(std::io::stdout()))
}

pub(crate) struct Target {
	path: Option<String>,
	atomic: bool,
	partial: Cell<Option<Partial>>,
}

impl Target {
	pub(crate) fn new(path: Option<&str>, inputs: &[&str]) -> Result<Self, Failure> {
		let Some(path) = path.filter(|&path| path != "-") else {
			return Ok(Self {
				path: None,
				atomic: false,
				partial: Cell::new(None),
			});
		};
		if let Ok(output) = fs::canonicalize(path)
			&& inputs
				.iter()
				.any(|input| fs::canonicalize(input).is_ok_and(|input| input == output))
		{
			return Err(Failure::new(format!(
				"the output \"{path}\" is also an input"
			)));
		}
		let regular = fs::metadata(path).is_ok_and(|metadata| metadata.is_file());
		let fresh = fs::symlink_metadata(path).is_err();
		Ok(Self {
			path: Some(path.to_owned()),
			atomic: regular || fresh,
			partial: Cell::new(None),
		})
	}

	pub(crate) fn open(&self) -> Result<Output, Failure> {
		let Some(path) = &self.path else {
			return Ok(stdout());
		};
		let sink = if self.atomic {
			let (partial, sink) =
				Partial::create(Path::new(path)).map_err(|error| file(path, "writing", &error))?;
			self.partial.set(Some(partial));
			sink
		} else {
			File::create(path).map_err(|error| file(path, "writing", &error))?
		};
		Ok(BufWriter::with_capacity(OUTPUT_BUFFER, Box::new(sink)))
	}

	pub(crate) fn commit(&self, mut output: Output) -> Result<(), Failure> {
		output.flush()?;
		drop(output);
		if let (Some(partial), Some(path)) = (self.partial.take(), &self.path) {
			partial
				.commit()
				.map_err(|error| file(path, "writing", &error))?;
		}
		Ok(())
	}
}
