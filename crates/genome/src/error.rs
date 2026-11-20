use std::io;

#[derive(Debug, thiserror::Error)]
pub enum Error {
	#[error(transparent)]
	Io(#[from] io::Error),
	#[error("unknown contig {0}")]
	UnknownContig(String),
	#[error("duplicate contig {0}")]
	DuplicateContig(String),
	#[error("invalid region {0}")]
	InvalidRegion(String),
	#[error("invalid fasta index line {line}")]
	InvalidIndex { line: usize },
	#[error("invalid fasta at byte {offset}")]
	InvalidFasta { offset: u64 },
	#[error("interval {start}-{end} is outside {contig}")]
	OutOfBounds {
		contig: String,
		start: u64,
		end: u64,
	},
}
