use std::path::PathBuf;

use gyrase_extract::Format;
use gyrase_jobs::{Cancel, Error, Report};

use crate::message;
use crate::prepared::Prepared;

pub struct Finished {
	pub prepared: Prepared,
	pub report: Report,
}

pub enum Screen {
	Pick(Option<String>),
	Ready(Prepared),
	Running(Prepared, Cancel),
	Stopping(Prepared),
	Done(Finished),
	Failed(Prepared, String),
}

impl Default for Screen {
	fn default() -> Self {
		Self::Pick(None)
	}
}

impl Screen {
	pub fn open(self, paths: &[PathBuf]) -> Self {
		match self {
			Self::Running(..) | Self::Stopping(_) => self,
			Self::Ready(prepared) | Self::Failed(prepared, _) => match Prepared::open(paths) {
				Ok(opened) => Self::Ready(opened),
				Err(message) => Self::Failed(prepared, message),
			},
			Self::Pick(_) | Self::Done(_) => match Prepared::open(paths) {
				Ok(opened) => Self::Ready(opened),
				Err(message) => Self::Pick(Some(message)),
			},
		}
	}

	pub fn clear(self) -> Self {
		match self {
			Self::Running(..) | Self::Stopping(_) => self,
			Self::Pick(_) | Self::Ready(_) | Self::Done(_) | Self::Failed(..) => Self::default(),
		}
	}

	pub fn toggle(self, index: usize) -> Self {
		let Some(&format) = Format::ALL.get(index) else {
			return self;
		};
		self.revise(|prepared| prepared.with(prepared.toggled(format), prepared.output.clone()))
	}

	pub fn relocate(self, folder: PathBuf) -> Self {
		self.revise(|prepared| prepared.with(prepared.formats.clone(), folder))
	}

	pub fn start(self, launch: impl FnOnce(&Prepared) -> Cancel) -> Self {
		match self {
			Self::Ready(prepared) | Self::Failed(prepared, _) => {
				let cancel = launch(&prepared);
				Self::Running(prepared, cancel)
			}
			other => other,
		}
	}

	pub fn cancel(self) -> Self {
		match self {
			Self::Running(prepared, cancel) => {
				cancel.cancel();
				Self::Stopping(prepared)
			}
			other => other,
		}
	}

	pub fn finish(self, outcome: Result<Report, Error>) -> Self {
		match self {
			Self::Running(prepared, _) | Self::Stopping(prepared) => match outcome {
				Ok(report) => Self::Done(Finished { prepared, report }),
				Err(Error::Cancelled) => Self::Ready(prepared),
				Err(error) => Self::Failed(prepared, message::text(&error)),
			},
			other => other,
		}
	}

	fn revise(self, change: impl FnOnce(&Prepared) -> Result<Prepared, String>) -> Self {
		match self {
			Self::Ready(prepared) | Self::Failed(prepared, _) => match change(&prepared) {
				Ok(revised) => Self::Ready(revised),
				Err(message) => Self::Failed(prepared, message),
			},
			other => other,
		}
	}
}
