use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT: AtomicU64 = AtomicU64::new(0);

#[derive(Debug)]
pub struct Partial {
	path: PathBuf,
	target: PathBuf,
	committed: bool,
}

impl Partial {
	pub fn create(target: &Path) -> io::Result<(Self, File)> {
		loop {
			let mut name = OsString::from(target);
			name.push(format!(
				".{}.{}.part",
				std::process::id(),
				NEXT.fetch_add(1, Ordering::Relaxed)
			));
			let path = PathBuf::from(name);
			match OpenOptions::new().write(true).create_new(true).open(&path) {
				Ok(file) => {
					let partial = Self {
						path,
						target: target.to_owned(),
						committed: false,
					};
					return Ok((partial, file));
				}
				Err(error) if error.kind() == io::ErrorKind::AlreadyExists => {}
				Err(error) => return Err(error),
			}
		}
	}

	pub fn path(&self) -> &Path {
		&self.path
	}

	pub fn commit(mut self) -> io::Result<()> {
		fs::rename(&self.path, &self.target)?;
		self.committed = true;
		Ok(())
	}
}

impl Drop for Partial {
	fn drop(&mut self) {
		if !self.committed {
			let _ = fs::remove_file(&self.path);
		}
	}
}
