mod parse;
mod render;

use std::io::{BufRead, Write};

use gyrase_io::Workers;

use crate::Error;
use crate::filter::{Filter, Scratch as FilterScratch};
use crate::header::Header;
use crate::par::ordered;
use crate::read::{LineBatch, Reader, RecordBatch};
use crate::record::Record;
use crate::region::Regions;
use crate::site::{Alleles, Site};
use crate::value::Loaded;
use crate::view::Condition;

const BATCH_BYTES: usize = 4 << 20;

#[derive(Debug, Clone, Default)]
pub struct QueryOptions {
	pub samples: Option<Vec<String>>,
	pub condition: Option<Condition>,
	pub regions: Option<Regions>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Field {
	Chrom,
	Pos,
	Pos0,
	End,
	End0,
	Id,
	Ref,
	Alt(Option<usize>),
	FirstAlt,
	Qual,
	Filter,
	InfoAll,
	FormatKeys,
	Info { key: u32, subscript: Option<usize> },
	Format { key: u32, subscript: Option<usize> },
	Genotype,
	TargetGenotype,
	Sample,
	Type,
	Line,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Segment {
	Literal(Vec<u8>),
	Field(Field),
	PerSample(Vec<Segment>),
}

#[derive(Debug)]
enum Compiled {
	Include(Filter),
	Exclude(Filter),
}

#[derive(Debug)]
pub struct Query {
	header: Header,
	segments: Vec<Segment>,
	samples: Vec<usize>,
	condition: Option<Compiled>,
	regions: Option<Regions>,
}

#[derive(Debug, Default)]
pub struct QueryOutput {
	bytes: Vec<u8>,
	scratch: FilterScratch,
	loaded: Loaded,
	alleles: Alleles,
	keys: Vec<u32>,
	record: Record,
	line: Vec<u8>,
}

impl Query {
	pub fn new(header: &Header, format: &str, options: QueryOptions) -> Result<Self, Error> {
		let mut samples: Vec<usize> = match &options.samples {
			None => (0..header.sample_count()).collect(),
			Some(names) => names
				.iter()
				.map(|name| {
					header
						.sample_index(name)
						.ok_or_else(|| Error::UndefinedSample(name.clone()))
				})
				.collect::<Result<_, _>>()?,
		};
		samples.sort_unstable();
		samples.dedup();
		let condition = match &options.condition {
			None => None,
			Some(Condition::Include(text)) => {
				Some(Compiled::Include(Filter::compile(text, header)?))
			}
			Some(Condition::Exclude(text)) => {
				Some(Compiled::Exclude(Filter::compile(text, header)?))
			}
		};
		Ok(Self {
			header: header.clone(),
			segments: parse::parse(format, header)?,
			samples,
			condition,
			regions: options.regions,
		})
	}

	pub fn row<S: Site>(&self, site: &S, out: &mut QueryOutput) -> Result<(), Error> {
		if let Some(regions) = &self.regions {
			let end = site.end(&self.header)?;
			if !regions.overlaps(site.contig(&self.header)?, site.start()?, end) {
				return Ok(());
			}
		}
		let rejected = match &self.condition {
			None => false,
			Some(Compiled::Include(filter)) => {
				!filter.test(site, &self.header, &mut out.scratch)?
			}
			Some(Compiled::Exclude(filter)) => filter.test(site, &self.header, &mut out.scratch)?,
		};
		if rejected {
			return Ok(());
		}
		render::render(&self.segments, &self.header, &self.samples, site, out)
	}

	pub fn lines(&self, batch: &LineBatch, out: &mut QueryOutput) -> Result<(), Error> {
		out.bytes.clear();
		for line in batch.lines() {
			self.row(&line?, out)?;
		}
		Ok(())
	}

	pub fn records(&self, batch: &RecordBatch, out: &mut QueryOutput) -> Result<(), Error> {
		out.bytes.clear();
		for record in batch.records() {
			self.row(&record?, out)?;
		}
		Ok(())
	}

	pub fn run<R: BufRead>(
		&self,
		reader: Reader<R>,
		sink: &mut impl Write,
		workers: &Workers,
	) -> Result<(), Error> {
		let drain = |out: &mut QueryOutput| sink.write_all(&out.bytes).map_err(Error::from);
		match reader {
			Reader::Vcf(mut vcf) => ordered(
				workers,
				QueryOutput::default,
				|batch: &mut LineBatch| vcf.read_batch(batch, BATCH_BYTES),
				|batch, out| self.lines(batch, out),
				drain,
			),
			Reader::Bcf(mut bcf) => ordered(
				workers,
				QueryOutput::default,
				|batch: &mut RecordBatch| bcf.read_batch(batch, BATCH_BYTES),
				|batch, out| self.records(batch, out),
				drain,
			),
		}
	}
}
