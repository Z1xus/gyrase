use std::num::NonZeroUsize;
use std::thread;

use sysinfo::{MemoryRefreshKind, RefreshKind, System};

use crate::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Resources {
	pub threads: NonZeroUsize,
	pub memory: u64,
}

impl Resources {
	pub fn detect() -> Result<Self, Error> {
		let system = System::new_with_specifics(
			RefreshKind::nothing().with_memory(MemoryRefreshKind::nothing().with_ram()),
		)?;
		let available = system.available_memory();
		let memory = system
			.cgroup_limits()
			.map_or(available, |limits| limits.free_memory.min(available));
		Ok(Self {
			threads: thread::available_parallelism().unwrap_or(NonZeroUsize::MIN),
			memory,
		})
	}
}
