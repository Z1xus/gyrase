mod combine;
mod group;
mod input;

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::io::BufRead;

use crate::Error;
use crate::header::Header;
use crate::read::Reader;
use crate::record::{CountKeys, Record};
use crate::value::ValueType;
use crate::write::{Encoded, Format, Formatter, Placed, Sink, Span, Writer};
use combine::Combiner;
use group::Groups;
use input::{Input, Place, Taken};

const FLUSH_RECORDS: usize = 4096;
const SUMMED: [&str; 2] = ["DP", "DP4"];
const SUMMED_WITHOUT_SAMPLES: [&str; 2] = ["AC", "AN"];

#[derive(Debug, Clone, Copy)]
pub struct MergeOptions {
	pub format: Format,
}

#[derive(Debug)]
pub(crate) struct Rule {
	key: u32,
	value_type: ValueType,
}

#[derive(Debug)]
pub struct Merge {
	header: Header,
	format: Format,
	rules: Vec<Rule>,
	skipped: Vec<u32>,
	counts: Option<CountKeys>,
	genotype: Option<u32>,
	pass: Option<u32>,
}

impl Merge {
	pub fn new(headers: &[&Header], options: MergeOptions) -> Result<Self, Error> {
		let (first, rest) = headers
			.split_first()
			.ok_or(Error::Incompatible("no inputs"))?;
		let header = rest
			.iter()
			.try_fold((*first).clone(), |joined, next| joined.join(next))?;
		let has_samples = header.sample_count() > 0;
		let mut names: Vec<&str> = SUMMED.to_vec();
		if !has_samples {
			names.extend(SUMMED_WITHOUT_SAMPLES);
		}
		names.sort_unstable();
		let rules = names
			.into_iter()
			.filter_map(|name| {
				let key = header.key(name)?;
				let value_type = header.type_of_info(key)?;
				matches!(value_type, ValueType::Integer | ValueType::Float)
					.then_some(Rule { key, value_type })
			})
			.collect();
		let genotype = header.key("GT");
		let (skipped, counts) = match (has_samples, genotype, header.key("AC"), header.key("AN")) {
			(true, genotype, allele_count, allele_number) => (
				[allele_count, allele_number]
					.into_iter()
					.flatten()
					.collect(),
				match (genotype, allele_count, allele_number) {
					(Some(genotype), Some(allele_count), Some(allele_number)) => Some(CountKeys {
						genotype,
						allele_count,
						allele_number,
					}),
					_ => None,
				},
			),
			(false, ..) => (Vec::new(), None),
		};
		Ok(Self {
			pass: header.key("PASS"),
			header,
			format: options.format,
			rules,
			skipped,
			counts,
			genotype,
		})
	}

	pub fn header(&self) -> &Header {
		&self.header
	}

	pub fn format(&self) -> Format {
		self.format
	}

	pub fn run<R: BufRead, S: Sink>(
		&self,
		readers: Vec<Reader<R>>,
		writer: &mut Writer<S>,
		mut placed: impl FnMut(Span, Placed<S::Mark>),
	) -> Result<u64, Error> {
		let mut inputs = readers
			.into_iter()
			.map(|reader| Input::new(reader, &self.header))
			.collect::<Result<Vec<_>, _>>()?;
		let samples: Vec<usize> = inputs.iter().map(Input::sample_count).collect();
		let mut taken: Vec<Taken> = inputs.iter().map(|_| Taken::default()).collect();
		let mut heap: BinaryHeap<Reverse<(Place, usize)>> = inputs
			.iter()
			.enumerate()
			.filter_map(|(index, input)| Some(Reverse((input.position()?, index))))
			.collect();
		let mut groups = Groups::default();
		let mut combiner = Combiner::new(self, &samples);
		let mut formatter = Formatter::new(self.format);
		let mut encoded = Encoded::default();
		let mut record = Record::new();
		let mut ready: Vec<usize> = Vec::new();
		let mut written = 0u64;
		while let Some(Reverse((at, first))) = heap.pop() {
			ready.clear();
			ready.push(first);
			while heap.peek().is_some_and(|Reverse((next, _))| *next == at) {
				if let Some(Reverse((_, index))) = heap.pop() {
					ready.push(index);
				}
			}
			for &index in &ready {
				let (Some(input), Some(slot)) = (inputs.get_mut(index), taken.get_mut(index))
				else {
					return Err(Error::Limits("inputs"));
				};
				input.take(at, slot)?;
				if let Some(next) = input.position() {
					heap.push(Reverse((next, index)));
				}
			}
			groups.build(&taken, &ready)?;
			for group in groups.ordered() {
				combiner.combine(self, &taken, group, &mut record)?;
				formatter.push_record(&self.header, &record, &mut encoded)?;
			}
			if encoded.len() >= FLUSH_RECORDS {
				written += encoded.len() as u64;
				writer.write(&encoded, &mut placed)?;
				encoded.clear();
			}
		}
		written += encoded.len() as u64;
		writer.write(&encoded, &mut placed)?;
		Ok(written)
	}
}
