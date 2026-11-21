use std::io::Write;
use std::sync::Arc;

use gyrase_genome::Reference;
use gyrase_io::{Level, Workers};
use gyrase_reads::cram::CramWriter;
use gyrase_reads::{BamWriter, Header, RecordRef, SamWriter};

use crate::failure::Failure;
use crate::shared::Output as Sink;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Format {
	Sam,
	Bam,
	Cram,
}

impl Format {
	pub fn named(name: &str) -> Option<Self> {
		match name.to_ascii_lowercase().as_str() {
			"sam" => Some(Self::Sam),
			"bam" => Some(Self::Bam),
			"cram" => Some(Self::Cram),
			_ => None,
		}
	}

	pub fn of_path(path: &str) -> Self {
		path.rsplit_once('.')
			.and_then(|(_, extension)| Self::named(extension))
			.unwrap_or(Self::Sam)
	}
}

pub enum Output {
	Sam(SamWriter<Sink>),
	Bam(BamWriter<Sink>),
	Cram(CramWriter<Sink>),
}

impl Output {
	pub fn new(
		sink: Sink,
		format: Format,
		header: &Header,
		with_header: bool,
		level: Level,
		reference: Option<&Arc<Reference>>,
		workers: &Workers,
	) -> Result<Self, Failure> {
		Ok(match format {
			Format::Sam if with_header => Self::Sam(SamWriter::new(sink, header)?),
			Format::Sam => Self::Sam(SamWriter::headerless(sink, header)),
			Format::Bam => Self::Bam(BamWriter::new(sink, header, level, workers)?),
			Format::Cram => {
				let reference = reference
					.ok_or_else(|| Failure::new("CRAM output requires a reference, use -T"))?;
				Self::Cram(CramWriter::new(
					sink,
					header,
					Arc::clone(reference),
					workers,
				)?)
			}
		})
	}

	pub fn write(&mut self, record: RecordRef<'_>) -> Result<(), Failure> {
		match self {
			Self::Sam(writer) => writer.write(record)?,
			Self::Bam(writer) => {
				writer.write(record)?;
			}
			Self::Cram(writer) => writer.write(record)?,
		}
		Ok(())
	}

	pub fn finish(self) -> Result<(), Failure> {
		let mut sink = match self {
			Self::Sam(writer) => writer.finish()?,
			Self::Bam(writer) => writer.finish()?.inner,
			Self::Cram(writer) => writer.finish()?.sink,
		};
		sink.flush()?;
		Ok(())
	}
}
