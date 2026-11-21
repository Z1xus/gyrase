use std::fs::File;
use std::io::{self, Read};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

pub(crate) struct Counted {
	inner: File,
	read: Arc<AtomicU64>,
}

impl Counted {
	pub(crate) fn new(inner: File, read: &Arc<AtomicU64>) -> Self {
		Self {
			inner,
			read: Arc::clone(read),
		}
	}
}

impl Read for Counted {
	fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
		let read = self.inner.read(buffer)?;
		self.read
			.fetch_add(u64::try_from(read).unwrap_or(0), Ordering::Relaxed);
		Ok(read)
	}
}
