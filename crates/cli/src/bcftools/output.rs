use std::io::{self, Write};

use gyrase_index::{Format as IndexFormat, Index};
use gyrase_io::{BgzfWriter, BlockPosition, Finished, Level, Workers};
use gyrase_variants::{
	Error, Format, Header, IndexCollector, IndexKind, Placed, PlainSink, Sink, Span, Writer,
};
use lexopt::{Parser, ValueExt};

use crate::failure::Failure;
use crate::shared::{Output, Target, workers};

#[derive(Clone, Copy, Default)]
pub(super) enum Kind {
	#[default]
	Vcf,
	VcfGz(Level),
	Bcf(Level),
}

impl Kind {
	fn parse(text: &str) -> Result<Self, Failure> {
		let invalid = || Failure::usage(format!("the output type \"{text}\" is not recognised"));
		let mut chars = text.chars();
		let kind = match chars.next().ok_or_else(invalid)? {
			'v' => Self::Vcf,
			'z' => Self::VcfGz(Level::DEFAULT),
			'b' => Self::Bcf(Level::DEFAULT),
			'u' => Self::Bcf(Level::STORED),
			_ => return Err(invalid()),
		};
		let Some(digit) = chars.next() else {
			return Ok(kind);
		};
		let level = digit
			.to_digit(10)
			.and_then(|digit| u8::try_from(digit).ok())
			.and_then(Level::new)
			.ok_or_else(invalid)?;
		if chars.next().is_some() {
			return Err(invalid());
		}
		Ok(match kind {
			Self::Vcf => Self::Vcf,
			Self::VcfGz(_) => Self::VcfGz(level),
			Self::Bcf(_) => Self::Bcf(level),
		})
	}

	pub(super) fn format(self) -> Format {
		match self {
			Self::Vcf | Self::VcfGz(_) => Format::Vcf,
			Self::Bcf(_) => Format::Bcf,
		}
	}
}

#[derive(Default)]
pub(super) struct Destination {
	pub(super) path: Option<String>,
	pub(super) kind: Kind,
	index: Option<IndexFormat>,
}

pub(super) struct Common {
	pub(super) destination: Destination,
	threads: Option<usize>,
}

impl Common {
	pub(super) fn new() -> Self {
		Self {
			destination: Destination::default(),
			threads: None,
		}
	}

	pub(super) fn output_type(&mut self, parser: &mut Parser) -> Result<(), Failure> {
		self.destination.kind = Kind::parse(&parser.value()?.string()?)?;
		Ok(())
	}

	pub(super) fn output(&mut self, parser: &mut Parser) -> Result<(), Failure> {
		self.destination.path = Some(parser.value()?.string()?).filter(|path| path != "-");
		Ok(())
	}

	pub(super) fn write_index(&mut self, parser: &mut Parser) -> Result<(), Failure> {
		let format = parser.optional_value().map(ValueExt::string).transpose()?;
		let format = format
			.as_deref()
			.map(|text| text.strip_prefix('=').unwrap_or(text));
		self.destination.index = Some(index_format(format)?);
		Ok(())
	}

	pub(super) fn target(&self, inputs: &[&str]) -> Result<Target, Failure> {
		Target::new(self.destination.path.as_deref(), inputs)
	}

	pub(super) fn thread_count(&mut self, parser: &mut Parser) -> Result<(), Failure> {
		self.threads = Some(parser.value()?.parse()?);
		Ok(())
	}

	pub(super) fn workers(&self) -> Result<Workers, Failure> {
		self.destination.check()?;
		workers(self.threads)
	}
}

pub(super) fn index_format(text: Option<&str>) -> Result<IndexFormat, Failure> {
	match text {
		None | Some("csi") => Ok(IndexFormat::Csi),
		Some("tbi") => Ok(IndexFormat::Tbi),
		Some(other) => Err(Failure::usage(format!(
			"the index format \"{other}\" is not recognised"
		))),
	}
}

pub(super) enum OutputSink {
	Plain(PlainSink<Output>),
	Bgzf(BgzfWriter<Output>),
}

pub(super) enum Done {
	Plain(Output),
	Bgzf(Finished<Output>),
}

impl Write for OutputSink {
	fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
		match self {
			Self::Plain(sink) => sink.write(bytes),
			Self::Bgzf(sink) => sink.write(bytes),
		}
	}

	fn flush(&mut self) -> io::Result<()> {
		match self {
			Self::Plain(sink) => sink.flush(),
			Self::Bgzf(sink) => sink.flush(),
		}
	}
}

impl Sink for OutputSink {
	type Mark = Option<BlockPosition>;
	type Done = Done;

	fn mark(&self) -> Self::Mark {
		match self {
			Self::Plain(_) => None,
			Self::Bgzf(sink) => Some(sink.mark()),
		}
	}

	fn end_block(&mut self) -> Result<(), Error> {
		match self {
			Self::Plain(sink) => sink.end_block(),
			Self::Bgzf(sink) => Sink::end_block(sink),
		}
	}

	fn finish(self) -> Result<Done, Error> {
		Ok(match self {
			Self::Plain(sink) => Done::Plain(sink.finish()?),
			Self::Bgzf(sink) => Done::Bgzf(Sink::finish(sink)?),
		})
	}
}

pub(super) struct Collector {
	inner: IndexCollector,
}

impl Collector {
	pub(super) fn new() -> Self {
		Self {
			inner: IndexCollector::new(),
		}
	}

	pub(super) fn mark(&mut self) -> impl FnMut(Span, Placed<Option<BlockPosition>>) + '_ {
		|span, placed| {
			if let (Some(begin), Some(end)) = (placed.begin, placed.end) {
				self.inner.push(span, Placed { begin, end });
			}
		}
	}
}

impl Destination {
	fn check(&self) -> Result<(), Failure> {
		match (&self.index, &self.path) {
			(Some(_), None) => Err(Failure::new(
				"the index cannot be written to standard output",
			)),
			(Some(_), Some(_)) if matches!(self.kind, Kind::Vcf) => Err(Failure::new(
				"an index needs compressed output, use -O z or -O b",
			)),
			(Some(IndexFormat::Tbi), Some(_)) if self.kind.format() == Format::Bcf => {
				Err(Failure::new("a tbi index cannot index BCF output"))
			}
			_ => Ok(()),
		}
	}

	pub(super) fn writer(
		&self,
		header: &Header,
		format: Format,
		workers: &Workers,
		target: &Target,
	) -> Result<Writer<OutputSink>, Failure> {
		Ok(Writer::new(self.sink(workers, target)?, header, format)?)
	}

	pub(super) fn bare(
		&self,
		workers: &Workers,
		target: &Target,
	) -> Result<Writer<OutputSink>, Failure> {
		Ok(Writer::without_header(self.sink(workers, target)?))
	}

	fn sink(&self, workers: &Workers, target: &Target) -> Result<OutputSink, Failure> {
		let out = target.open()?;
		Ok(match self.kind {
			Kind::Vcf => OutputSink::Plain(PlainSink::new(out)),
			Kind::VcfGz(level) | Kind::Bcf(level) => {
				OutputSink::Bgzf(BgzfWriter::new(out, level, workers))
			}
		})
	}

	pub(super) fn finish(
		&self,
		writer: Writer<OutputSink>,
		collector: &Collector,
		header: &Header,
		workers: &Workers,
		target: &Target,
	) -> Result<(), Failure> {
		match writer.finish()? {
			Done::Plain(out) => target.commit(out),
			Done::Bgzf(finished) => {
				let index = match (self.index, &self.path) {
					(Some(format), Some(path)) => {
						let kind = if self.kind.format() == Format::Bcf {
							IndexKind::Bcf
						} else {
							IndexKind::Vcf
						};
						let index = collector.inner.build(&finished, header, kind)?;
						let extension = if format == IndexFormat::Tbi {
							"tbi"
						} else {
							"csi"
						};
						Some((index, format, format!("{path}.{extension}")))
					}
					_ => None,
				};
				target.commit(finished.inner)?;
				match index {
					Some((index, format, path)) => write_index(&index, format, &path, &[], workers),
					None => Ok(()),
				}
			}
		}
	}
}

pub(super) fn write_index(
	index: &Index,
	format: IndexFormat,
	path: &str,
	inputs: &[&str],
	workers: &Workers,
) -> Result<(), Failure> {
	let target = Target::new(Some(path), inputs)?;
	target.commit(index.write(format, target.open()?, workers)?)
}
