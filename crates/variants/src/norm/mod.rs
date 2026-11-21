mod align;
mod kinds;
mod order;
mod split;

use gyrase_genome::{ContigId, Reference};
use gyrase_io::Workers;
use std::io::BufRead;

use self::align::{Aligner, Outcome, Place};
use self::kinds::Kinds;
use self::split::Splitter;
use crate::Error;
use crate::bcf::BcfRecord;
use crate::header::Header;
use crate::par::ordered;
use crate::read::{LineBatch, Reader, RecordBatch};
use crate::record::Record;
use crate::vcf::VcfLine;
use crate::write::{Encoded, Format, Formatter, Placed, Sink, Span, Writer};

pub use self::order::Duplicates;
use self::order::{Dedup, Window};

const BATCH_BYTES: usize = 4 << 20;
const DEFAULT_WINDOW: u32 = 1000;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Mismatch {
	#[default]
	Fail,
	Keep,
	Skip,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Splitting {
	Snps,
	Indels,
	Both,
	Any,
}

#[derive(Debug)]
pub struct NormOptions {
	pub reference: Option<Reference>,
	pub mismatch: Mismatch,
	pub split: Option<Splitting>,
	pub duplicates: Duplicates,
	pub window: u32,
	pub format: Format,
}

impl Default for NormOptions {
	fn default() -> Self {
		Self {
			reference: None,
			mismatch: Mismatch::Fail,
			split: None,
			duplicates: Duplicates::Keep,
			window: DEFAULT_WINDOW,
			format: Format::Vcf,
		}
	}
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NormStats {
	pub total: u64,
	pub split: u64,
	pub realigned: u64,
	pub skipped: u64,
	pub duplicates: u64,
}

impl NormStats {
	fn absorb(&mut self, other: &mut Self) {
		self.total += other.total;
		self.split += other.split;
		self.realigned += other.realigned;
		self.skipped += other.skipped;
		self.duplicates += other.duplicates;
		*other = Self::default();
	}
}

#[derive(Debug, Default)]
struct Pool {
	ready: Vec<Record>,
	spare: Vec<Record>,
}

impl Pool {
	fn take_from(&mut self, record: &mut Record) {
		let mut fresh = self.spare.pop().unwrap_or_default();
		std::mem::swap(&mut fresh, record);
		self.ready.push(fresh);
	}
}

#[derive(Debug, Clone, Copy)]
struct Mark {
	records: usize,
	checked: bool,
}

#[derive(Debug, Default)]
pub struct NormOutput {
	source: Record,
	part: Record,
	kept: Pool,
	marks: Vec<Mark>,
	aligner: Aligner,
	splitter: Splitter,
	stats: NormStats,
}

#[derive(Debug)]
pub struct Norm {
	header: Header,
	reference: Option<Reference>,
	contigs: Vec<Option<ContigId>>,
	mismatch: Mismatch,
	split: Option<Splitting>,
	duplicates: Duplicates,
	window: u32,
	format: Format,
	genotype: Option<u32>,
	end_key: Option<u32>,
	svlen_key: Option<u32>,
}

#[derive(Debug)]
struct Flush {
	window: Window,
	dedup: Dedup,
	formatter: Formatter,
	encoded: Encoded,
	stats: NormStats,
}

impl Norm {
	pub fn new(header: &Header, options: NormOptions) -> Self {
		let contigs = header
			.contig_names()
			.map(|name| {
				options
					.reference
					.as_ref()
					.and_then(|reference| reference.contigs().resolve(name))
			})
			.collect();
		Self {
			header: header.clone(),
			reference: options.reference,
			contigs,
			mismatch: options.mismatch,
			split: options.split,
			duplicates: options.duplicates,
			window: options.window,
			format: options.format,
			genotype: header.key("GT"),
			end_key: header.key("END"),
			svlen_key: header.key("SVLEN"),
		}
	}

	pub fn header(&self) -> &Header {
		&self.header
	}

	pub fn format(&self) -> Format {
		self.format
	}

	pub fn output(&self) -> NormOutput {
		NormOutput::default()
	}

	fn splits(&self, record: &Record) -> bool {
		let Some(split) = self.split else {
			return false;
		};
		if record.allele_count() <= 2 {
			return false;
		}
		let kinds = Kinds::of(record);
		match split {
			Splitting::Snps => kinds.has_snp(),
			Splitting::Indels => kinds.has_indel(),
			Splitting::Both | Splitting::Any => true,
		}
	}

	fn location(&self, record: &Record) -> (String, u64) {
		let name = self
			.header
			.contig(record.contig())
			.map(|contig| contig.name.clone())
			.unwrap_or_default();
		(name, u64::from(record.start()) + 1)
	}

	fn normalize(
		&self,
		aligner: &mut Aligner,
		record: &mut Record,
		stats: &mut NormStats,
	) -> Result<bool, Error> {
		let Some(reference) = &self.reference else {
			return Ok(true);
		};
		let index = usize::try_from(record.contig()).map_err(|_| Error::Limits("contig"))?;
		let contig = self
			.contigs
			.get(index)
			.copied()
			.flatten()
			.ok_or_else(|| Error::UndefinedContig(self.location(record).0))?;
		let place = Place {
			reference,
			contig,
			end_key: self.end_key,
		};
		match aligner.realign(&place, record)? {
			Outcome::Unchanged => Ok(true),
			Outcome::Changed => {
				stats.realigned += 1;
				Ok(true)
			}
			Outcome::Mismatch => match self.mismatch {
				Mismatch::Fail => {
					let (contig, position) = self.location(record);
					Err(Error::ReferenceMismatch { contig, position })
				}
				Mismatch::Keep => Ok(true),
				Mismatch::Skip => {
					stats.skipped += 1;
					Ok(false)
				}
			},
			Outcome::Duplicate => match self.mismatch {
				Mismatch::Fail => {
					let (contig, position) = self.location(record);
					Err(Error::DuplicateAllele { contig, position })
				}
				Mismatch::Keep | Mismatch::Skip => Ok(true),
			},
		}
	}

	fn process(&self, out: &mut NormOutput) -> Result<(), Error> {
		out.stats.total += 1;
		if !self.splits(&out.source) {
			let kept = self.normalize(&mut out.aligner, &mut out.source, &mut out.stats)?;
			if kept {
				out.kept.take_from(&mut out.source);
			}
			out.marks.push(Mark {
				records: usize::from(kept),
				checked: kept,
			});
			return Ok(());
		}
		out.stats.split += 1;
		let before = out.kept.ready.len();
		let mut checked = true;
		for alt in 0..out.source.allele_count().saturating_sub(1) {
			out.splitter
				.split(&self.header, self.genotype, &out.source, alt, &mut out.part)?;
			if !self.normalize(&mut out.aligner, &mut out.part, &mut out.stats)? {
				checked = false;
				break;
			}
			out.kept.take_from(&mut out.part);
		}
		out.marks.push(Mark {
			records: out.kept.ready.len() - before,
			checked,
		});
		Ok(())
	}

	pub fn line(&self, line: &VcfLine<'_>, out: &mut NormOutput) -> Result<(), Error> {
		out.source.read_vcf(line, &self.header)?;
		self.process(out)
	}

	pub fn bcf(&self, source: &BcfRecord<'_>, out: &mut NormOutput) -> Result<(), Error> {
		out.source.read_bcf(source)?;
		self.process(out)
	}

	pub fn lines(&self, batch: &LineBatch, out: &mut NormOutput) -> Result<(), Error> {
		out.marks.clear();
		out.kept.spare.append(&mut out.kept.ready);
		for line in batch.lines() {
			self.line(&line?, out)?;
		}
		Ok(())
	}

	pub fn records(&self, batch: &RecordBatch, out: &mut NormOutput) -> Result<(), Error> {
		out.marks.clear();
		out.kept.spare.append(&mut out.kept.ready);
		for record in batch.records() {
			self.bcf(&record?, out)?;
		}
		Ok(())
	}

	fn flusher(&self) -> Flush {
		Flush {
			window: Window::new(self.window),
			dedup: Dedup::new(self.duplicates, self.svlen_key),
			formatter: Formatter::new(self.format),
			encoded: Encoded::default(),
			stats: NormStats::default(),
		}
	}

	fn absorb(&self, flush: &mut Flush, out: &mut NormOutput) -> Result<(), Error> {
		let mut records = out.kept.ready.drain(..);
		for mark in &out.marks {
			for record in records.by_ref().take(mark.records) {
				flush.window.insert(record);
			}
			if mark.checked && flush.window.spans_positions() {
				let due = flush.window.due();
				self.emit(flush, due)?;
			}
		}
		flush.stats.absorb(&mut out.stats);
		out.kept.spare.append(flush.window.spares());
		Ok(())
	}

	fn emit(&self, flush: &mut Flush, count: usize) -> Result<(), Error> {
		for _ in 0..count {
			let Some(record) = flush.window.pop() else {
				break;
			};
			if flush.dedup.drops(&record) {
				flush.stats.duplicates += 1;
			} else {
				flush
					.formatter
					.push_record(&self.header, &record, &mut flush.encoded)?;
			}
			flush.window.recycle(record);
		}
		Ok(())
	}

	pub fn run<R: BufRead, S: Sink>(
		&self,
		reader: Reader<R>,
		writer: &mut Writer<S>,
		workers: &Workers,
		mut placed: impl FnMut(Span, Placed<S::Mark>),
	) -> Result<NormStats, Error> {
		let mut flush = self.flusher();
		let mut drain = |out: &mut NormOutput| {
			self.absorb(&mut flush, out)?;
			writer.write(&flush.encoded, &mut placed)?;
			flush.encoded.clear();
			Ok(())
		};
		match reader {
			Reader::Vcf(mut vcf) => ordered(
				workers,
				|| self.output(),
				|batch: &mut LineBatch| vcf.read_batch(batch, BATCH_BYTES),
				|batch, out| self.lines(batch, out),
				&mut drain,
			)?,
			Reader::Bcf(mut bcf) => ordered(
				workers,
				|| self.output(),
				|batch: &mut RecordBatch| bcf.read_batch(batch, BATCH_BYTES),
				|batch, out| self.records(batch, out),
				&mut drain,
			)?,
		}
		let remaining = flush.window.len();
		self.emit(&mut flush, remaining)?;
		writer.write(&flush.encoded, &mut placed)?;
		Ok(flush.stats)
	}
}
