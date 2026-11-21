use std::io::{self, Write};
use std::path::Path;

use gyrase_index::{Format as IndexFormat, Index};
use gyrase_io::{BgzfWriter, BlockPosition, Finished, Level, Workers};
use gyrase_variants::{
	Error, Format, Header, IndexCollector, IndexKind, Placed, PlainSink, Sink, Span, Writer,
};
use lexopt::{Parser, ValueExt};

use crate::failure::{Failure, file};
use crate::shared::{Output, create, workers};

#[derive(Clone, Copy)]
pub struct Kind {
	pub format: Format,
	pub compressed: bool,
	pub level: Level,
}

impl Default for Kind {
	fn default() -> Self {
		Self {
			format: Format::Vcf,
			compressed: false,
			level: Level::DEFAULT,
		}
	}
}

impl Kind {
	pub fn parse(text: &str) -> Result<Self, Failure> {
		let invalid = || Failure::usage(format!("the output type \"{text}\" is not recognised"));
		let mut chars = text.chars();
		let mut kind = match chars.next().ok_or_else(invalid)? {
			'v' => Self::default(),
			'z' => Self {
				compressed: true,
				..Self::default()
			},
			'b' => Self {
				format: Format::Bcf,
				compressed: true,
				..Self::default()
			},
			'u' => Self {
				format: Format::Bcf,
				compressed: true,
				level: Level::STORED,
			},
			_ => return Err(invalid()),
		};
		if let Some(digit) = chars.next() {
			let level = digit
				.to_digit(10)
				.and_then(|digit| u8::try_from(digit).ok())
				.and_then(Level::new)
				.ok_or_else(invalid)?;
			if chars.next().is_some() {
				return Err(invalid());
			}
			if kind.compressed {
				kind.level = level;
			}
		}
		Ok(kind)
	}
}

#[derive(Default)]
pub struct Destination {
	pub path: Option<String>,
	pub kind: Kind,
	pub index: Option<IndexFormat>,
}

pub struct Common {
	pub destination: Destination,
	pub threads: Option<usize>,
}

impl Common {
	pub fn new() -> Self {
		Self {
			destination: Destination::default(),
			threads: None,
		}
	}

	pub fn output_type(&mut self, parser: &mut Parser) -> Result<(), Failure> {
		self.destination.kind = Kind::parse(&parser.value()?.string()?)?;
		Ok(())
	}

	pub fn output(&mut self, parser: &mut Parser) -> Result<(), Failure> {
		self.destination.path = Some(parser.value()?.string()?);
		Ok(())
	}

	pub fn write_index(&mut self, parser: &mut Parser) -> Result<(), Failure> {
		let format = parser.optional_value().map(ValueExt::string).transpose()?;
		self.destination.index = Some(index_format(format.as_deref())?);
		Ok(())
	}

	pub fn thread_count(&mut self, parser: &mut Parser) -> Result<(), Failure> {
		self.threads = Some(parser.value()?.parse()?);
		Ok(())
	}

	pub fn workers(&self) -> Result<Workers, Failure> {
		self.destination.check()?;
		workers(self.threads)
	}
}

pub fn index_format(text: Option<&str>) -> Result<IndexFormat, Failure> {
	match text {
		None | Some("csi") => Ok(IndexFormat::Csi),
		Some("tbi") => Ok(IndexFormat::Tbi),
		Some(other) => Err(Failure::usage(format!(
			"the index format \"{other}\" is not recognised"
		))),
	}
}

pub enum Dest {
	Plain(PlainSink<Output>),
	Bgzf(BgzfWriter<Output>),
}

pub enum Done {
	Plain(Output),
	Bgzf(Finished<Output>),
}

impl Write for Dest {
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

impl Sink for Dest {
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

pub struct Collector {
	inner: IndexCollector,
}

impl Collector {
	pub fn new() -> Self {
		Self {
			inner: IndexCollector::new(),
		}
	}

	pub fn mark(&mut self) -> impl FnMut(Span, Placed<Option<BlockPosition>>) + '_ {
		|span, placed| {
			if let (Some(begin), Some(end)) = (placed.begin, placed.end) {
				self.inner.push(span, Placed { begin, end });
			}
		}
	}
}

impl Destination {
	pub fn check(&self) -> Result<(), Failure> {
		match (&self.index, &self.path) {
			(Some(_), None) => Err(Failure::new(
				"the index cannot be written to standard output",
			)),
			(Some(_), Some(_)) if !self.kind.compressed => Err(Failure::new(
				"an index needs compressed output, use -O z or -O b",
			)),
			(Some(IndexFormat::Tbi), Some(_)) if self.kind.format == Format::Bcf => {
				Err(Failure::new("a tbi index cannot index BCF output"))
			}
			_ => Ok(()),
		}
	}

	pub fn writer(
		&self,
		header: &Header,
		format: Format,
		workers: &Workers,
	) -> Result<Writer<Dest>, Failure> {
		Ok(Writer::new(self.sink(workers)?, header, format)?)
	}

	pub fn bare(&self, workers: &Workers) -> Result<Writer<Dest>, Failure> {
		Ok(Writer::without_header(self.sink(workers)?))
	}

	fn sink(&self, workers: &Workers) -> Result<Dest, Failure> {
		let out = create(self.path.as_deref())?;
		Ok(if self.kind.compressed {
			Dest::Bgzf(BgzfWriter::new(out, self.kind.level, workers))
		} else {
			Dest::Plain(PlainSink::new(out))
		})
	}

	pub fn finish(
		&self,
		writer: Writer<Dest>,
		collector: &Collector,
		header: &Header,
		workers: &Workers,
	) -> Result<(), Failure> {
		match writer.finish()? {
			Done::Plain(mut out) => out.flush()?,
			Done::Bgzf(mut finished) => {
				if let (Some(format), Some(path)) = (self.index, &self.path) {
					let kind = if self.kind.format == Format::Bcf {
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
					write_index(&index, format, &format!("{path}.{extension}"), workers)?;
				}
				finished.inner.flush()?;
			}
		}
		Ok(())
	}
}

pub fn write_index(
	index: &Index,
	format: IndexFormat,
	path: &str,
	workers: &Workers,
) -> Result<(), Failure> {
	let out =
		std::fs::File::create(Path::new(path)).map_err(|error| file(path, "writing", &error))?;
	index.write(format, out, workers)?;
	Ok(())
}
