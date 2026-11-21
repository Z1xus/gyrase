use std::io;

#[derive(Debug, thiserror::Error)]
pub enum Error {
	#[error(transparent)]
	Io(#[from] io::Error),
	#[error(transparent)]
	Genome(#[from] gyrase_genome::Error),
	#[error(transparent)]
	Workers(#[from] gyrase_io::Error),
	#[error(transparent)]
	Record(#[from] gyrase_reads::RecordError),
	#[error("alignment refers to read {0}, which is not in the given input")]
	ForeignBatch(u32),
	#[error("invalid seeding parameters: {0}")]
	InvalidSeeding(&'static str),
	#[error("reference has {length} bases, the index limit is {limit}")]
	ReferenceTooLong { length: u64, limit: u64 },
	#[error("cancelled")]
	Cancelled,
	#[error("reference has no seedable sequence")]
	EmptyReference,
	#[error("not a synapton index file")]
	NotAnIndex,
	#[error("index file version {found} is not supported, expected {expected}")]
	IndexVersion { found: u32, expected: u32 },
	#[error("index file is truncated or inconsistent: {0}")]
	CorruptIndex(&'static str),
	#[error("index file checksum mismatch in {0}")]
	IndexChecksum(&'static str),
	#[error("read {index} has length {length}, the limit is {limit}")]
	ReadTooLong {
		index: usize,
		length: usize,
		limit: usize,
	},
	#[error("read {index} has {bases} bases but {qualities} qualities")]
	QualityLength {
		index: usize,
		bases: usize,
		qualities: usize,
	},
}

impl From<std::num::TryFromIntError> for Error {
	fn from(_: std::num::TryFromIntError) -> Self {
		Self::CorruptIndex("size overflow")
	}
}
