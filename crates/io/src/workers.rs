use std::num::NonZeroUsize;
use std::sync::Arc;

use rayon::{ThreadPool, ThreadPoolBuilder};

use crate::Error;

#[derive(Debug, Clone)]
pub struct Workers {
	pool: Arc<ThreadPool>,
}

impl Workers {
	pub fn new(threads: NonZeroUsize) -> Result<Self, Error> {
		let pool = ThreadPoolBuilder::new()
			.num_threads(threads.get())
			.thread_name(|index| format!("gyrase-{index}"))
			.build()?;
		Ok(Self {
			pool: Arc::new(pool),
		})
	}

	pub fn available() -> Result<Self, Error> {
		Self::new(std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN))
	}

	pub fn threads(&self) -> usize {
		self.pool.current_num_threads()
	}

	pub fn install<R: Send>(&self, work: impl FnOnce() -> R + Send) -> R {
		self.pool.install(work)
	}
}
