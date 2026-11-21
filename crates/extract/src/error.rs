use std::io;

#[derive(Debug, thiserror::Error)]
pub enum Error {
	#[error(transparent)]
	Io(#[from] io::Error),
	#[error(transparent)]
	Decode(#[from] gyrase_io::Error),
	#[error("not a gyrase panel")]
	NotPanel,
	#[error("unsupported panel version {0}")]
	UnsupportedVersion(u16),
	#[error("panel is truncated")]
	Truncated,
	#[error("panel is corrupt: {0}")]
	Corrupt(&'static str),
	#[error("panel has more than {0} sites")]
	TooManySites(usize),
	#[error("invalid site name {0:?}")]
	InvalidName(String),
	#[error("invalid position for site {0:?}")]
	InvalidPosition(String),
	#[error("duplicate site {name} at {chromosome}:{position}")]
	DuplicateSite {
		name: String,
		chromosome: &'static str,
		position: u32,
	},
	#[error("expected {expected} calls, got {actual}")]
	CallCount { expected: usize, actual: usize },
	#[error("build is not supported for export")]
	UnsupportedBuild,
	#[error("invalid timestamp")]
	InvalidTimestamp,
}
