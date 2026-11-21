use std::fmt;
use std::io;

use crate::HeaderError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Codec {
	Gzip,
	Bzip2,
	Lzma,
	Rans4x8,
	RansNx16,
	Arith,
	Fqzcomp,
	Tok3,
}

impl fmt::Display for Codec {
	fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
		formatter.write_str(match self {
			Self::Gzip => "gzip",
			Self::Bzip2 => "bzip2",
			Self::Lzma => "lzma",
			Self::Rans4x8 => "rANS 4x8",
			Self::RansNx16 => "rANS Nx16",
			Self::Arith => "arithmetic",
			Self::Fqzcomp => "fqzcomp",
			Self::Tok3 => "name tokeniser",
		})
	}
}

#[derive(Debug, thiserror::Error)]
pub enum Error {
	#[error(transparent)]
	Io(#[from] io::Error),
	#[error(transparent)]
	Header(#[from] HeaderError),
	#[error(transparent)]
	Reference(#[from] gyrase_genome::Error),
	#[error(transparent)]
	Record(#[from] crate::RecordError),
	#[error(transparent)]
	Index(#[from] gyrase_index::Error),
	#[error("input ended early")]
	Truncated,
	#[error("{0} data is corrupt")]
	Corrupt(Codec),
	#[error("integer is too large")]
	IntegerOverflow,
	#[error("a length of {len} bytes is over the limit")]
	TooLarge { len: u64 },
	#[error("block holds {actual} bytes but its header says {expected}")]
	SizeMismatch { expected: usize, actual: usize },
	#[error("file does not start with the CRAM magic")]
	BadMagic,
	#[error("CRAM version {major}.{minor} is not supported")]
	UnsupportedVersion { major: u8, minor: u8 },
	#[error("checksum of the {0} does not match")]
	ChecksumMismatch(&'static str),
	#[error("malformed {0}")]
	Malformed(&'static str),
	#[error("{0} cannot be written to CRAM")]
	Unsupported(&'static str),
	#[error("a data series that the records need is not defined")]
	MissingSeries,
	#[error("the records need a reference sequence but none is available")]
	ReferenceRequired,
	#[error("reference {0} is not in the reference file")]
	UnknownReference(String),
	#[error("reference checksum of {0} does not match the slice")]
	ReferenceChecksum(String),
}
