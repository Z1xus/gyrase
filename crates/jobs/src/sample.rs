use std::collections::BTreeMap;
use std::path::Path;

use crate::{Error, Input, Kind};

const SEPARATORS: [char; 3] = ['_', '.', '-'];
const EXTENSIONS: [&str; 13] = [
	".gz", ".bgz", ".fastq", ".fq", ".sam", ".bam", ".cram", ".vcf", ".bcf", ".fa", ".fasta",
	".fna", ".txt",
];

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Lane {
	Single(Input),
	Paired(Input, Input),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Source {
	Reads(Vec<Lane>),
	Alignments(Input),
	Variants(Input),
	Reference(Input),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sample {
	pub name: String,
	pub source: Source,
}

impl Sample {
	pub fn path(&self) -> &Path {
		match &self.source {
			Source::Alignments(input) | Source::Variants(input) | Source::Reference(input) => {
				&input.path
			}
			Source::Reads(lanes) => lanes.first().map_or(Path::new(""), |lane| match lane {
				Lane::Single(input) | Lane::Paired(input, _) => &input.path,
			}),
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mate {
	First,
	Second,
}

struct Read {
	sample: String,
	lane: String,
	mate: Option<Mate>,
	input: Input,
}

pub fn samples(inputs: Vec<Input>) -> Result<Vec<Sample>, Error> {
	let mut samples = Vec::new();
	let mut reads = Vec::new();
	for input in inputs {
		let name = stem(&input.path);
		match input.kind {
			Kind::Fastq => reads.push(Read::new(input)),
			Kind::Sam | Kind::Bam | Kind::Cram => {
				samples.push(Sample {
					name,
					source: Source::Alignments(input),
				});
			}
			Kind::Vcf | Kind::Bcf => samples.push(Sample {
				name,
				source: Source::Variants(input),
			}),
			Kind::Fasta => samples.push(Sample {
				name,
				source: Source::Reference(input),
			}),
		}
	}
	let mut lanes: BTreeMap<String, BTreeMap<String, Vec<Read>>> = BTreeMap::new();
	for read in reads {
		lanes
			.entry(read.sample.clone())
			.or_default()
			.entry(read.lane.clone())
			.or_default()
			.push(read);
	}
	for (name, by_lane) in lanes {
		let lanes = by_lane
			.into_values()
			.map(lane)
			.collect::<Result<Vec<_>, _>>()?;
		samples.push(Sample {
			name,
			source: Source::Reads(lanes.into_iter().flatten().collect()),
		});
	}
	Ok(samples)
}

fn lane(mut reads: Vec<Read>) -> Result<Vec<Lane>, Error> {
	reads.sort_by_key(|read| read.mate.map(|mate| mate == Mate::Second));
	match reads.as_slice() {
		[
			Read {
				mate: Some(Mate::First),
				..
			},
			Read {
				mate: Some(Mate::Second),
				..
			},
		] => {
			let mut pair = reads.into_iter().map(|read| read.input);
			Ok(pair
				.next()
				.zip(pair.next())
				.map(|(first, second)| Lane::Paired(first, second))
				.into_iter()
				.collect())
		}
		_ => reads
			.into_iter()
			.map(|read| match read.mate {
				None => Ok(Lane::Single(read.input)),
				Some(_) => Err(Error::Unpaired(read.input.path)),
			})
			.collect(),
	}
}

impl Read {
	fn new(input: Input) -> Self {
		let stem = stem(&input.path);
		let tokens: Vec<(usize, &str)> = stem
			.match_indices(SEPARATORS)
			.map(|(at, _)| at + 1)
			.map(|start| {
				(
					start,
					stem[start..].split(SEPARATORS).next().unwrap_or_default(),
				)
			})
			.collect();
		let found = tokens.iter().rev().find_map(|&(start, token)| match token {
			"R1" | "1" => Some((start, token, Mate::First)),
			"R2" | "2" => Some((start, token, Mate::Second)),
			_ => None,
		});
		let sample_end = tokens
			.iter()
			.find(|(start, token)| {
				found.is_some_and(|(at, ..)| at == *start) || is_lane(token) || is_index(token)
			})
			.map_or(stem.len(), |&(start, _)| start - 1);
		let sample = stem[..sample_end].to_owned();
		match found {
			Some((start, token, mate)) => {
				let lane = format!("{}{}", &stem[..start], &stem[start + token.len()..]);
				Self {
					sample,
					lane,
					mate: Some(mate),
					input,
				}
			}
			None => Self {
				sample,
				lane: stem.clone(),
				mate: None,
				input,
			},
		}
	}
}

fn is_lane(token: &str) -> bool {
	token
		.strip_prefix('L')
		.is_some_and(|digits| digits.len() == 3 && digits.bytes().all(|byte| byte.is_ascii_digit()))
}

fn is_index(token: &str) -> bool {
	token.strip_prefix('S').is_some_and(|digits| {
		!digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit())
	})
}

fn stem(path: &Path) -> String {
	let mut name = path
		.file_name()
		.map(|name| name.to_string_lossy().into_owned())
		.unwrap_or_default();
	while let Some(stripped) = EXTENSIONS
		.iter()
		.find_map(|extension| name.strip_suffix(extension))
	{
		name = stripped.to_owned();
	}
	name
}
