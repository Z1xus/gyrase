use std::io;

#[derive(Debug, thiserror::Error)]
pub enum Error {
	#[error(transparent)]
	Io(#[from] io::Error),
	#[error("not a bgzf block at byte {offset}")]
	NotBgzf { offset: u64 },
	#[error("truncated bgzf block at byte {offset}")]
	Truncated { offset: u64 },
	#[error("corrupt bgzf block at byte {offset}")]
	Corrupt { offset: u64 },
	#[error("compression failed")]
	Compression,
	#[error("worker threads could not start")]
	Workers(#[from] rayon::ThreadPoolBuildError),
	#[error("a worker thread stopped")]
	Disconnected,
}

impl From<Error> for io::Error {
	fn from(error: Error) -> Self {
		match error {
			Error::Io(error) => error,
			other => io::Error::new(io::ErrorKind::InvalidData, other),
		}
	}
}
