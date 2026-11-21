use gyrase_io::Workers;
use std::io::BufRead;

use crate::Error;
use crate::bcf::BcfRecord;
use crate::filter::{Filter, Scratch as FilterScratch};
use crate::header::Header;
use crate::par::ordered;
use crate::read::{LineBatch, Reader, RecordBatch};
use crate::record::{CountKeys, Record};
use crate::region::Regions;
use crate::site::Site;
use crate::vcf::VcfLine;
use crate::write::{Encoded, Format, Formatter, Placed, Sink, Span, Writer};

const BATCH_BYTES: usize = 4 << 20;

#[derive(Debug, Clone)]
pub enum Condition {
	Include(String),
	Exclude(String),
}

#[derive(Debug, Clone, Default)]
pub enum SampleChoice {
	#[default]
	All,
	Keep(Vec<String>),
	Drop(Vec<String>),
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Fidelity {
	#[default]
	Canonical,
	Verbatim,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Counts {
	#[default]
	Recompute,
	Keep,
}

#[derive(Debug, Clone)]
pub struct ViewOptions {
	pub samples: SampleChoice,
	pub counts: Counts,
	pub condition: Option<Condition>,
	pub regions: Option<Regions>,
	pub fidelity: Fidelity,
	pub format: Format,
}

impl Default for ViewOptions {
	fn default() -> Self {
		Self {
			samples: SampleChoice::All,
			counts: Counts::Recompute,
			condition: None,
			regions: None,
			fidelity: Fidelity::Canonical,
			format: Format::Vcf,
		}
	}
}

#[derive(Debug)]
enum Compiled {
	Include(Filter),
	Exclude(Filter),
}

#[derive(Debug)]
pub struct View {
	header: Header,
	keep: Option<Vec<usize>>,
	condition: Option<Compiled>,
	regions: Option<Regions>,
	fidelity: Fidelity,
	format: Format,
	input_header: Header,
	counts: Option<CountKeys>,
}

#[derive(Debug)]
pub struct ViewOutput {
	count_scratch: Vec<i32>,
	encoded: Encoded,
	formatter: Formatter,
	record: Record,
	scratch: FilterScratch,
}

impl View {
	pub fn new(input: &Header, options: ViewOptions) -> Result<Self, Error> {
		let keep = match &options.samples {
			SampleChoice::All => None,
			SampleChoice::Keep(names) => Some(sample_indices(input, names)?),
			SampleChoice::Drop(names) => {
				let dropped = sample_indices(input, names)?;
				Some(
					(0..input.sample_count())
						.filter(|index| !dropped.contains(index))
						.collect(),
				)
			}
		};
		let (header, counts) = match (&keep, options.counts) {
			(Some(indices), Counts::Recompute) => count_header(input.with_samples(indices)?)?,
			(Some(indices), Counts::Keep) => (input.with_samples(indices)?, None),
			(None, _) => (input.clone(), None),
		};
		let header = if header.sample_count() == 0 {
			header.without_formats()
		} else {
			header
		};
		let condition = match &options.condition {
			None => None,
			Some(Condition::Include(text)) => {
				Some(Compiled::Include(Filter::compile(text, input)?))
			}
			Some(Condition::Exclude(text)) => {
				Some(Compiled::Exclude(Filter::compile(text, input)?))
			}
		};
		Ok(Self {
			header,
			keep,
			condition,
			regions: options.regions,
			fidelity: options.fidelity,
			format: options.format,
			input_header: input.clone(),
			counts,
		})
	}

	pub fn header(&self) -> &Header {
		&self.header
	}

	pub fn format(&self) -> Format {
		self.format
	}

	pub fn output(&self) -> ViewOutput {
		ViewOutput {
			count_scratch: Vec::new(),
			encoded: Encoded::default(),
			formatter: Formatter::new(self.format),
			record: Record::new(),
			scratch: FilterScratch::default(),
		}
	}

	fn accepts<S: Site>(&self, site: &S, scratch: &mut FilterScratch) -> Result<bool, Error> {
		if let Some(regions) = &self.regions {
			let contig = site.contig(&self.input_header)?;
			let end = site.end(&self.input_header)?;
			if !regions.overlaps(contig, site.start()?, end) {
				return Ok(false);
			}
		}
		match &self.condition {
			None => Ok(true),
			Some(Compiled::Include(filter)) => filter.test(site, &self.input_header, scratch),
			Some(Compiled::Exclude(filter)) => {
				Ok(!filter.test(site, &self.input_header, scratch)?)
			}
		}
	}

	fn verbatim(&self) -> bool {
		self.fidelity == Fidelity::Verbatim && self.keep.is_none()
	}

	fn finish_record(&self, out: &mut ViewOutput) -> Result<(), Error> {
		let all_dropped = self.keep.as_ref().is_some_and(Vec::is_empty);
		if let Some(keep) = &self.keep
			&& !all_dropped
		{
			out.record.keep_samples(keep)?;
		}
		if let Some(keys) = self.counts {
			out.record.update_counts(keys, &mut out.count_scratch)?;
		}
		if all_dropped {
			out.record.keep_samples(&[])?;
		}
		out.formatter
			.push_record(&self.header, &out.record, &mut out.encoded)
	}

	pub fn line(&self, line: &VcfLine<'_>, out: &mut ViewOutput) -> Result<(), Error> {
		if !self.accepts(line, &mut out.scratch)? {
			return Ok(());
		}
		if self.verbatim() && self.format == Format::Vcf {
			let span = Span {
				contig: line.contig(&self.input_header)?,
				start: line.start()?,
				end: line.end(&self.input_header)?,
			};
			out.formatter.push_line(line.text(), span, &mut out.encoded);
			return Ok(());
		}
		out.record.read_vcf(line, &self.input_header)?;
		self.finish_record(out)
	}

	pub fn bcf(&self, source: &BcfRecord<'_>, out: &mut ViewOutput) -> Result<(), Error> {
		if !self.accepts(source, &mut out.scratch)? {
			return Ok(());
		}
		if self.verbatim() && self.format == Format::Bcf {
			let span = Span {
				contig: source.contig()?,
				start: source.start()?,
				end: source.end(&self.input_header)?,
			};
			return out
				.formatter
				.push_bcf(source.shared(), source.indiv(), span, &mut out.encoded);
		}
		if self.verbatim() && self.format == Format::Vcf {
			let span = Span {
				contig: source.contig()?,
				start: source.start()?,
				end: source.end(&self.input_header)?,
			};
			return out
				.formatter
				.push_bcf_as_vcf(&self.header, source, span, &mut out.encoded);
		}
		out.record.read_bcf(source)?;
		self.finish_record(out)
	}

	pub fn lines(&self, batch: &LineBatch, out: &mut ViewOutput) -> Result<(), Error> {
		out.encoded.clear();
		for line in batch.lines() {
			self.line(&line?, out)?;
		}
		Ok(())
	}

	pub fn records(&self, batch: &RecordBatch, out: &mut ViewOutput) -> Result<(), Error> {
		out.encoded.clear();
		for record in batch.records() {
			self.bcf(&record?, out)?;
		}
		Ok(())
	}

	pub fn run<R: BufRead, S: Sink>(
		&self,
		reader: Reader<R>,
		writer: &mut Writer<S>,
		workers: &Workers,
		mut placed: impl FnMut(Span, Placed<S::Mark>),
	) -> Result<u64, Error> {
		let mut written = 0u64;
		let mut drain = |out: &mut ViewOutput| {
			written += out.encoded.len() as u64;
			writer.write(&out.encoded, &mut placed)
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
		Ok(written)
	}
}

fn sample_indices(header: &Header, names: &[String]) -> Result<Vec<usize>, Error> {
	names
		.iter()
		.map(|name| {
			header
				.sample_index(name)
				.ok_or_else(|| Error::UndefinedSample(name.clone()))
		})
		.collect()
}

const COUNT_LINES: [(&str, &str); 2] = [
	(
		"AC",
		"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">",
	),
	(
		"AN",
		"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
	),
];

fn count_header(header: Header) -> Result<(Header, Option<CountKeys>), Error> {
	let missing: Vec<&str> = COUNT_LINES
		.iter()
		.filter(|(name, _)| header.type_of_info_by_name(name).is_none())
		.map(|(_, line)| *line)
		.collect();
	let header = if missing.is_empty() {
		header
	} else {
		header.with_lines(&missing)?
	};
	let keys = match (header.key("GT"), header.key("AC"), header.key("AN")) {
		(Some(genotype), Some(allele_count), Some(allele_number)) => Some(CountKeys {
			genotype,
			allele_count,
			allele_number,
		}),
		_ => None,
	};
	Ok((header, keys))
}
