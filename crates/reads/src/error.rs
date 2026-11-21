use std::io;

#[derive(Debug, thiserror::Error)]
pub enum Error {
	#[error(transparent)]
	Io(#[from] io::Error),
	#[error(transparent)]
	Bgzf(#[from] gyrase_io::Error),
	#[error(transparent)]
	Genome(#[from] gyrase_genome::Error),
	#[error(transparent)]
	Header(#[from] HeaderError),
	#[error(transparent)]
	Record(#[from] RecordError),
	#[error(transparent)]
	Bam(#[from] BamError),
	#[error(transparent)]
	Sam(#[from] SamError),
	#[error(transparent)]
	Fastq(#[from] FastqError),
	#[error(transparent)]
	Index(#[from] gyrase_index::Error),
	#[error(transparent)]
	Sort(#[from] SortError),
	#[error(transparent)]
	Mark(#[from] MarkError),
	#[error("cancelled")]
	Cancelled,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum HeaderError {
	#[error("header text is not valid UTF-8")]
	InvalidUtf8,
	#[error("header line {line} is malformed")]
	MalformedLine { line: usize },
	#[error("header line {line} has no {tag} field")]
	MissingField { line: usize, tag: &'static str },
	#[error("header line {line} has an invalid reference length")]
	InvalidLength { line: usize },
	#[error("duplicate reference {0}")]
	DuplicateReference(String),
	#[error("reference {index} in the text differs from the binary header")]
	ReferenceMismatch { index: usize },
	#[error("header is larger than {limit} bytes")]
	TooLarge { limit: usize },
	#[error("reference {0} is too long for BAM")]
	ReferenceTooLong(String),
	#[error("binary header is truncated")]
	Truncated,
	#[error("reference name is empty or too long")]
	InvalidReferenceName,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum RecordError {
	#[error("record is {len} bytes, shorter than the fixed fields")]
	TooShort { len: usize },
	#[error("read name is empty or not terminated")]
	InvalidName,
	#[error("field lengths exceed the record")]
	Truncated,
	#[error("invalid cigar operation code {code}")]
	InvalidCigarOp { code: u8 },
	#[error("invalid coordinate or reference id")]
	InvalidCoordinate,
	#[error("reference id {contig} is not in the header")]
	ContigOutOfRange { contig: u32 },
	#[error("malformed tag")]
	MalformedTag,
	#[error("read name is longer than 254 bytes or contains a NUL")]
	NameTooLong,
	#[error("cigar operation length is too large")]
	CigarLengthTooLarge,
	#[error("qualities and bases differ in length")]
	LengthMismatch,
	#[error("quality value is out of range")]
	InvalidQuality,
	#[error("sequence is too long")]
	SequenceTooLong,
	#[error("tag array is too long")]
	ArrayTooLong,
	#[error("tag string contains a control character")]
	InvalidTagString,
	#[error("integer tag value is out of range")]
	IntegerOutOfRange,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum BamError {
	#[error("not a BAM file")]
	NotBam,
	#[error("BAM data ends inside a record")]
	TruncatedRecord,
	#[error("record of {len} bytes is outside the allowed size")]
	RecordSize { len: u64 },
	#[error("batch is larger than 4 GiB")]
	BatchTooLarge,
	#[error("query chunk does not start inside the data")]
	MisalignedChunk,
	#[error("records are not on a block boundary")]
	UnexpectedBlock,
	#[error("a written record has no known file offset")]
	UnresolvedOffset,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SamField {
	Name,
	Flag,
	Reference,
	Position,
	MapQ,
	Cigar,
	MateReference,
	MatePosition,
	TemplateLength,
	Sequence,
	Qualities,
	Tag,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SamError {
	#[error("line {line}: too few fields")]
	TooFewFields { line: u64 },
	#[error("line {line}: invalid {field:?}")]
	Invalid { line: u64, field: SamField },
	#[error("line {line}: unknown reference {name}")]
	UnknownReference { line: u64, name: String },
	#[error("line {line}: header line after records")]
	LateHeader { line: u64 },
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum FastqError {
	#[error("record {record}: header line does not start with @")]
	MissingHeader { record: u64 },
	#[error("record {record}: separator line does not start with +")]
	MissingSeparator { record: u64 },
	#[error("record {record}: bases and qualities differ in length")]
	LengthMismatch { record: u64 },
	#[error("record {record}: file ends inside the record")]
	Truncated { record: u64 },
	#[error("record {record}: read names differ between mates")]
	PairNames { record: u64 },
	#[error("the two files have a different number of reads")]
	PairCount,
	#[error("record {record} is larger than the batch limit")]
	RecordTooLarge { record: u64 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum SortError {
	#[error("cannot reserve {bytes} bytes for the sort buffer")]
	Reserve { bytes: usize },
	#[error("a spilled run ends inside a record")]
	TruncatedRun,
	#[error("a record is larger than the sort memory budget")]
	RecordTooLarge,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum MarkError {
	#[error("input is not in coordinate order")]
	NotSorted,
}
