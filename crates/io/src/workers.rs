use std::num::NonZeroUsize;
use std::sync::Arc;
use std::thread::{self, JoinHandle};
use std::time::Duration;

use rayon::{ThreadPool, ThreadPoolBuilder};

use crate::Error;

const IDLE_WAIT: Duration = Duration::from_micros(50);

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
		Self::new(thread::available_parallelism().unwrap_or(NonZeroUsize::MIN))
	}

	pub fn threads(&self) -> usize {
		self.pool.current_num_threads()
	}

	pub fn install<R: Send>(&self, work: impl FnOnce() -> R + Send) -> R {
		self.pool.install(work)
	}

	pub fn join<T>(&self, handle: JoinHandle<T>) -> thread::Result<T> {
		if self.pool.current_thread_index().is_some() {
			while !handle.is_finished() {
				if !matches!(rayon::yield_now(), Some(rayon::Yield::Executed)) {
					thread::sleep(IDLE_WAIT);
				}
			}
		}
		handle.join()
	}
}
