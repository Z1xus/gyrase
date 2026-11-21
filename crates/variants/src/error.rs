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
	Index(#[from] gyrase_index::Error),
	#[error("invalid header: {0}")]
	Header(&'static str),
	#[error("invalid vcf record: {0}")]
	Vcf(&'static str),
	#[error("invalid bcf record: {0}")]
	Bcf(&'static str),
	#[error("index chunk is inconsistent with the file: {0}")]
	Chunk(&'static str),
	#[error("not a bcf file")]
	NotBcf,
	#[error("unsupported bcf version {major}.{minor}")]
	BcfVersion { major: u8, minor: u8 },
	#[error("contig {0} is not defined in the header")]
	UndefinedContig(String),
	#[error("tag {0} is not defined in the header")]
	UndefinedTag(String),
	#[error("sample {0} is not in the header")]
	UndefinedSample(String),
	#[error("invalid expression at byte {position}: {reason}")]
	Expression {
		position: usize,
		reason: &'static str,
	},
	#[error("invalid format string at byte {position}: {reason}")]
	FormatString {
		position: usize,
		reason: &'static str,
	},
	#[error("inputs are not compatible: {0}")]
	Incompatible(&'static str),
	#[error("reference allele does not match the reference at {contig}:{position}")]
	ReferenceMismatch { contig: String, position: u64 },
	#[error("duplicate alleles at {contig}:{position}")]
	DuplicateAllele { contig: String, position: u64 },
	#[error("wrong number of values in {tag} at {contig}:{position}")]
	ValueCount {
		tag: String,
		contig: String,
		position: u64,
	},
	#[error("sample {0} is in more than one input")]
	DuplicateSample(String),
	#[error("records are not sorted")]
	Unsorted,
	#[error("{0} exceeds what bcf can store")]
	Limits(&'static str),
}
