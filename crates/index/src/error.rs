use std::io;

use crate::Format;

#[derive(Debug, thiserror::Error)]
pub enum Error {
	#[error(transparent)]
	Io(#[from] io::Error),
	#[error(transparent)]
	Bgzf(#[from] gyrase_io::Error),
	#[error("not an index file")]
	UnknownFormat,
	#[error("index file is truncated")]
	Truncated,
	#[error("index count {count} does not fit the file")]
	InvalidCount { count: u64 },
	#[error("index is larger than {limit} bytes")]
	TooLarge { limit: u64 },
	#[error("duplicate bin {bin}")]
	DuplicateBin { bin: u32 },
	#[error("invalid index geometry: min shift {min_shift}, depth {depth}")]
	InvalidGeometry { min_shift: u32, depth: u32 },
	#[error("invalid tabix header")]
	InvalidTabixHeader,
	#[error("index cannot be written as {format:?}")]
	IncompatibleFormat { format: Format },
	#[error("start {start} follows start {previous} on contig {contig}")]
	Unsorted {
		contig: u32,
		previous: u64,
		start: u64,
	},
	#[error("contig {contig} follows contig {previous}")]
	ContigOrder { contig: u32, previous: u32 },
	#[error("placed record after unplaced records")]
	PlacedAfterUnplaced,
	#[error("invalid interval {start}..{end}")]
	InvalidInterval { start: u64, end: u64 },
	#[error("position {position} is beyond the limit {limit} of this index")]
	PositionOutOfRange { position: u64, limit: u64 },
	#[error("contig {contig} is outside the limit of {limit} references")]
	ReferenceOutOfRange { contig: u32, limit: usize },
	#[error("invalid crai line {line}")]
	InvalidCrai { line: usize },
}
