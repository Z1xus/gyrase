use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::Error;

#[derive(Debug, Clone, Default)]
pub struct Cancel(Arc<AtomicBool>);

impl Cancel {
	pub fn cancel(&self) {
		self.0.store(true, Ordering::Relaxed);
	}

	pub(crate) fn flag(&self) -> &AtomicBool {
		&self.0
	}

	pub fn check(&self) -> Result<(), Error> {
		if self.0.load(Ordering::Relaxed) {
			Err(Error::Cancelled)
		} else {
			Ok(())
		}
	}
}
