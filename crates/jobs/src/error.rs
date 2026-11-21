use std::io;
use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum Error {
	#[error("{path}: {source}")]
	Open { path: PathBuf, source: io::Error },
	#[error("{path}: {source}")]
	Read {
		path: PathBuf,
		source: gyrase_io::Error,
	},
	#[error("{0}: not a FASTQ, FASTA, SAM, BAM, CRAM, VCF or BCF file")]
	UnknownKind(PathBuf),
	#[error("{0}: no matching mate file")]
	Unpaired(PathBuf),
	#[error("system information: {0}")]
	System(#[from] sysinfo::Error),
	#[error("no data folder")]
	NoDataDir,
	#[error("download failed: {0}")]
	Download(#[from] Box<ureq::Error>),
	#[error("downloaded {} does not match its build", .0.name())]
	Mismatch(gyrase_genome::Build),
	#[error(transparent)]
	Genome(#[from] gyrase_genome::Error),
	#[error("{0}: no BAM or CRAM index next to it")]
	MissingIndex(PathBuf),
	#[error("{0}: not a BAM or CRAM file")]
	NotAlignments(PathBuf),
	#[error("contig missing from the reference")]
	UnknownContig,
	#[error(transparent)]
	Io(#[from] gyrase_io::Error),
	#[error(transparent)]
	Reads(#[from] gyrase_reads::Error),
	#[error(transparent)]
	Cram(#[from] gyrase_reads::cram::Error),
	#[error(transparent)]
	Record(#[from] gyrase_reads::RecordError),
	#[error(transparent)]
	Index(#[from] gyrase_index::Error),
	#[error(transparent)]
	Variants(#[from] gyrase_variants::Error),
	#[error(transparent)]
	Extract(#[from] gyrase_extract::Error),
	#[error("cancelled")]
	Cancelled,
}
