use std::path::PathBuf;

use gyrase_extract::Format;
use gyrase_jobs::{Cancel, Error, Report};

use crate::message;
use crate::prepared::Prepared;

pub struct Finished {
	pub prepared: Prepared,
	pub report: Report,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Exit {
	Stay,
	Quit,
}

pub enum Screen {
	Pick(Option<String>),
	Ready(Prepared),
	Running(Prepared, Cancel),
	Stopping(Prepared, Exit),
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
			Self::Running(..) | Self::Stopping(..) => self,
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
			Self::Running(..) | Self::Stopping(..) => self,
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
		self.stop(Exit::Stay)
	}

	pub fn close(self) -> Self {
		self.stop(Exit::Quit)
	}

	pub fn exit(&self) -> Exit {
		match self {
			Self::Stopping(_, exit) => *exit,
			_ => Exit::Stay,
		}
	}

	fn stop(self, exit: Exit) -> Self {
		match self {
			Self::Running(prepared, cancel) => {
				cancel.cancel();
				Self::Stopping(prepared, exit)
			}
			Self::Stopping(prepared, stopping) => Self::Stopping(prepared, stopping.max(exit)),
			other => other,
		}
	}

	pub fn finish(self, outcome: Result<Report, Error>) -> Self {
		match self {
			Self::Running(prepared, _) | Self::Stopping(prepared, _) => match outcome {
				Ok(report) => Self::Done(Finished { prepared, report }),
				Err(Error::Cancelled) => Self::Ready(prepared.replanned()),
				Err(error) => Self::Failed(prepared.replanned(), message::text(&error)),
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
