use std::ffi::OsString;
use std::fs::{self, File};
use std::path::{Path, PathBuf};

use crate::Error;

#[derive(Debug)]
pub(crate) struct Partial {
	path: PathBuf,
	target: PathBuf,
}

impl Partial {
	pub(crate) fn create(target: &Path) -> Result<(Self, File), Error> {
		let partial = Self {
			path: suffixed(target, &format!(".{}.part", std::process::id())),
			target: target.to_owned(),
		};
		let file = File::create(&partial.path).map_err(|source| Error::Write {
			path: target.to_owned(),
			source,
		})?;
		Ok((partial, file))
	}

	pub(crate) fn path(&self) -> &Path {
		&self.path
	}

	pub(crate) fn commit(self) -> Result<(), Error> {
		fs::rename(&self.path, &self.target).map_err(|source| Error::Write {
			path: self.target.clone(),
			source,
		})
	}
}

impl Drop for Partial {
	fn drop(&mut self) {
		let _ = fs::remove_file(&self.path);
	}
}

pub(crate) fn suffixed(path: &Path, suffix: &str) -> PathBuf {
	let mut name = OsString::from(path);
	name.push(suffix);
	PathBuf::from(name)
}
