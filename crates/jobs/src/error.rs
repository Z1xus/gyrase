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
	#[error("cancelled")]
	Cancelled,
}
