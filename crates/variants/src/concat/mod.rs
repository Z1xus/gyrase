use gyrase_io::Workers;
use std::io::BufRead;

use crate::Error;
use crate::bcf::BcfRecord;
use crate::header::Header;
use crate::par::ordered;
use crate::read::{LineBatch, Reader, RecordBatch};
use crate::record::Record;
use crate::remap::Remap;
use crate::site::Site;
use crate::vcf::VcfLine;
use crate::view::Fidelity;
use crate::write::{Encoded, Format, Formatter, Placed, Sink, Span, Writer};

const BATCH_BYTES: usize = 4 << 20;

#[derive(Debug, Clone, Copy)]
pub struct ConcatOptions {
	pub fidelity: Fidelity,
	pub format: Format,
}

#[derive(Debug)]
pub struct Concat {
	header: Header,
	fidelity: Fidelity,
	format: Format,
}

#[derive(Debug)]
pub struct ConcatInput<'a> {
	concat: &'a Concat,
	remap: Remap,
}

#[derive(Debug)]
pub struct ConcatOutput {
	encoded: Encoded,
	formatter: Formatter,
	record: Record,
}

impl Concat {
	pub fn new(headers: &[&Header], options: ConcatOptions) -> Result<Self, Error> {
		let (first, rest) = headers
			.split_first()
			.ok_or(Error::Incompatible("no inputs"))?;
		let header = rest
			.iter()
			.try_fold((*first).clone(), |merged, header| merged.merge(header))?;
		Ok(Self {
			header,
			fidelity: options.fidelity,
			format: options.format,
		})
	}

	pub fn header(&self) -> &Header {
		&self.header
	}

	pub fn format(&self) -> Format {
		self.format
	}

	pub fn input(&self, header: &Header) -> Result<ConcatInput<'_>, Error> {
		Ok(ConcatInput {
			concat: self,
			remap: Remap::new(header, &self.header)?,
		})
	}

	pub fn output(&self) -> ConcatOutput {
		ConcatOutput {
			encoded: Encoded::default(),
			formatter: Formatter::new(self.format),
			record: Record::new(),
		}
	}

	pub fn run<R: BufRead, S: Sink>(
		&self,
		inputs: Vec<Reader<R>>,
		writer: &mut Writer<S>,
		workers: &Workers,
		mut placed: impl FnMut(Span, Placed<S::Mark>),
	) -> Result<u64, Error> {
		let mut written = 0u64;
		let mut drain = |out: &mut ConcatOutput| {
			written += out.encoded.len() as u64;
			writer.write(&out.encoded, &mut placed)
		};
		for reader in inputs {
			let input = self.input(reader.header())?;
			match reader {
				Reader::Vcf(mut vcf) => ordered(
					workers,
					|| self.output(),
					|batch: &mut LineBatch| vcf.read_batch(batch, BATCH_BYTES),
					|batch, out| input.lines(batch, out),
					&mut drain,
				)?,
				Reader::Bcf(mut bcf) => ordered(
					workers,
					|| self.output(),
					|batch: &mut RecordBatch| bcf.read_batch(batch, BATCH_BYTES),
					|batch, out| input.records(batch, out),
					&mut drain,
				)?,
			}
		}
		Ok(written)
	}
}

impl ConcatInput<'_> {
	fn verbatim(&self, format: Format) -> bool {
		self.concat.fidelity == Fidelity::Verbatim
			&& self.concat.format == format
			&& self.remap.is_identity()
	}

	pub fn line(&self, line: &VcfLine<'_>, out: &mut ConcatOutput) -> Result<(), Error> {
		if self.verbatim(Format::Vcf) {
			let span = Span {
				contig: line.contig(&self.concat.header)?,
				start: line.start()?,
				end: line.end(&self.concat.header)?,
			};
			out.formatter.push_line(line.text(), span, &mut out.encoded);
			return Ok(());
		}
		out.record.read_vcf(line, &self.concat.header)?;
		out.formatter
			.push_record(&self.concat.header, &out.record, &mut out.encoded)
	}

	pub fn bcf(&self, source: &BcfRecord<'_>, out: &mut ConcatOutput) -> Result<(), Error> {
		if self.verbatim(Format::Bcf) {
			let span = Span {
				contig: source.contig()?,
				start: source.start()?,
				end: source.end(&self.concat.header)?,
			};
			return out
				.formatter
				.push_bcf(source.shared(), source.indiv(), span, &mut out.encoded);
		}
		if self.verbatim(Format::Vcf) {
			let span = Span {
				contig: source.contig()?,
				start: source.start()?,
				end: source.end(&self.concat.header)?,
			};
			return out.formatter.push_bcf_as_vcf(
				&self.concat.header,
				source,
				span,
				&mut out.encoded,
			);
		}
		out.record.read_bcf(source)?;
		self.remap.apply(&mut out.record)?;
		out.formatter
			.push_record(&self.concat.header, &out.record, &mut out.encoded)
	}

	pub fn lines(&self, batch: &LineBatch, out: &mut ConcatOutput) -> Result<(), Error> {
		out.encoded.clear();
		for line in batch.lines() {
			self.line(&line?, out)?;
		}
		Ok(())
	}

	pub fn records(&self, batch: &RecordBatch, out: &mut ConcatOutput) -> Result<(), Error> {
		out.encoded.clear();
		for record in batch.records() {
			self.bcf(&record?, out)?;
		}
		Ok(())
	}
}
