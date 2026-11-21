use std::ffi::OsString;
use std::fs::File;
use std::path::{Path, PathBuf};

use crate::Error;

#[derive(Debug)]
pub(crate) struct Partial {
	inner: gyrase_io::Partial,
	target: PathBuf,
}

impl Partial {
	pub(crate) fn create(target: &Path) -> Result<(Self, File), Error> {
		let (inner, file) = gyrase_io::Partial::create(target).map_err(|source| Error::Write {
			path: target.to_owned(),
			source,
		})?;
		let partial = Self {
			inner,
			target: target.to_owned(),
		};
		Ok((partial, file))
	}

	pub(crate) fn path(&self) -> &Path {
		self.inner.path()
	}

	pub(crate) fn commit(self) -> Result<(), Error> {
		self.inner.commit().map_err(|source| Error::Write {
			path: self.target,
			source,
		})
	}
}

pub(crate) fn suffixed(path: &Path, suffix: &str) -> PathBuf {
	let mut name = OsString::from(path);
	name.push(suffix);
	PathBuf::from(name)
}
