use std::fs::File;
use std::path::{Path, PathBuf};

use gyrase_io::{Compression, peek};

use crate::Error;

const SNIFF_LEN: u64 = 1 << 16;
const SAM_COLUMNS: usize = 11;
const SAM_HEADER_TAGS: [[u8; 2]; 5] = [*b"HD", *b"SQ", *b"RG", *b"PG", *b"CO"];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Kind {
	Fastq,
	Fasta,
	Sam,
	Bam,
	Cram,
	Vcf,
	Bcf,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Input {
	pub path: PathBuf,
	pub kind: Kind,
	pub compression: Compression,
	pub size: u64,
}

impl Input {
	pub fn open(path: &Path) -> Result<Self, Error> {
		let opened = |source| Error::Open {
			path: path.to_owned(),
			source,
		};
		let file = File::open(path).map_err(opened)?;
		let size = file.metadata().map_err(opened)?.len();
		let (compression, head) = peek(file, SNIFF_LEN).map_err(|source| Error::Read {
			path: path.to_owned(),
			source,
		})?;
		let kind = Kind::of(&head).ok_or_else(|| Error::UnknownKind(path.to_owned()))?;
		Ok(Self {
			path: path.to_owned(),
			kind,
			compression,
			size,
		})
	}
}

impl Kind {
	fn of(head: &[u8]) -> Option<Self> {
		let first_line = head.split(|&byte| byte == b'\n').next().unwrap_or_default();
		match head {
			[b'B', b'A', b'M', 1, ..] => Some(Self::Bam),
			[b'C', b'R', b'A', b'M', ..] => Some(Self::Cram),
			[b'B', b'C', b'F', 2, ..] => Some(Self::Bcf),
			[b'@', first, second, b'\t', ..] if SAM_HEADER_TAGS.contains(&[*first, *second]) => {
				Some(Self::Sam)
			}
			[b'@', ..] => Some(Self::Fastq),
			[b'>', ..] => Some(Self::Fasta),
			_ if head.starts_with(b"##fileformat=VCF") => Some(Self::Vcf),
			_ if first_line.split(|&byte| byte == b'\t').count() >= SAM_COLUMNS => Some(Self::Sam),
			_ => None,
		}
	}
}
