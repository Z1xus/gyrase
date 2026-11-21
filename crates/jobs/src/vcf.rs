use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use gyrase_call::{Call, Genotype};
use gyrase_index::Format as IndexFormat;
use gyrase_io::{BgzfWriter, Level, Workers};
use gyrase_variants::{
	Encoded, Format, Formatter, Header, IndexCollector, IndexKind, Record, Writer,
};

use crate::{Alignments, Error};

type Sink = BgzfWriter<BufWriter<File>>;

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Debug)]
pub struct Vcf {
	path: PathBuf,
	header: Header,
	keys: Keys,
	writer: Writer<Sink>,
	formatter: Formatter,
	record: Record,
	encoded: Encoded,
	collector: IndexCollector,
	text: Vec<u8>,
}

#[derive(Debug, Clone, Copy)]
struct Keys {
	pass: u32,
	depth: u32,
	genotype: u32,
	genotype_quality: u32,
	allele_depths: u32,
	likelihoods: u32,
}

impl Vcf {
	pub fn create(path: &Path, alignments: &Alignments, workers: &Workers) -> Result<Self, Error> {
		let header = header(alignments)?;
		let keys = Keys {
			pass: key(&header, "PASS")?,
			depth: key(&header, "DP")?,
			genotype: key(&header, "GT")?,
			genotype_quality: key(&header, "GQ")?,
			allele_depths: key(&header, "AD")?,
			likelihoods: key(&header, "PL")?,
		};
		let file = File::create(path).map_err(|source| Error::Open {
			path: path.to_owned(),
			source,
		})?;
		let sink = BgzfWriter::new(BufWriter::new(file), Level::DEFAULT, workers);
		Ok(Self {
			path: path.to_owned(),
			writer: Writer::new(sink, &header, Format::Vcf)?,
			header,
			keys,
			formatter: Formatter::new(Format::Vcf),
			record: Record::new(),
			encoded: Encoded::default(),
			collector: IndexCollector::new(),
			text: Vec::new(),
		})
	}

	pub fn write(&mut self, contig: u32, calls: &[Call]) -> Result<(), Error> {
		self.encoded.clear();
		for call in calls {
			self.fill(contig, call)?;
			self.formatter
				.push_record(&self.header, &self.record, &mut self.encoded)?;
		}
		let collector = &mut self.collector;
		self.writer
			.write(&self.encoded, |span, placed| collector.push(span, placed))?;
		Ok(())
	}

	pub fn finish(self, workers: &Workers) -> Result<PathBuf, Error> {
		let finished = self.writer.finish()?;
		let index = self
			.collector
			.build(&finished, &self.header, IndexKind::Vcf)?;
		let mut index_path = self.path.clone().into_os_string();
		index_path.push(".tbi");
		let index_path = PathBuf::from(index_path);
		let file = File::create(&index_path).map_err(|source| Error::Open {
			path: index_path,
			source,
		})?;
		index.write(IndexFormat::Tbi, BufWriter::new(file), workers)?;
		Ok(self.path)
	}

	fn fill(&mut self, contig: u32, call: &Call) -> Result<(), Error> {
		let keys = self.keys;
		let record = &mut self.record;
		record.clear();
		record.set_site(
			contig,
			u32::try_from(call.position).map_err(|_| gyrase_variants::Error::Limits("position"))?,
		);
		record.set_alleles(
			std::iter::once(call.reference.as_slice())
				.chain(call.alternates.iter().map(Vec::as_slice)),
		);
		record.set_qual(Some(single(call.quality)));
		record.push_filter(keys.pass);
		record.push_info_ints(keys.depth, &[saturate(call.depth)]);
		self.text.clear();
		genotype_text(call.genotype, &mut self.text);
		record.push_format_text(keys.genotype, self.text.len(), &self.text);
		record.push_format_ints(
			keys.genotype_quality,
			1,
			&[i32::from(call.genotype_quality)],
		);
		let depths: Vec<i32> = call.allele_depths.iter().copied().map(saturate).collect();
		record.push_format_ints(keys.allele_depths, depths.len(), &depths);
		let likelihoods: Vec<i32> = call.likelihoods.iter().copied().map(saturate).collect();
		record.push_format_ints(keys.likelihoods, likelihoods.len(), &likelihoods);
		Ok(())
	}
}

fn header(alignments: &Alignments) -> Result<Header, Error> {
	let source = format!("##source=gyrase {VERSION}");
	let contigs: Vec<String> = alignments
		.header()
		.references()
		.iter()
		.map(|reference| {
			format!(
				"##contig=<ID={},length={}>",
				reference.name(),
				reference.length()
			)
		})
		.collect();
	let lines: Vec<&str> = ["##fileformat=VCFv4.2", source.as_str()]
		.into_iter()
		.chain(contigs.iter().map(String::as_str))
		.chain([
			"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
			"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
			"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">",
			"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">",
		])
		.collect();
	Ok(Header::from_lines(&lines, vec![alignments.sample()])?)
}

fn key(header: &Header, name: &str) -> Result<u32, Error> {
	header
		.key(name)
		.ok_or(Error::Variants(gyrase_variants::Error::Header(
			"missing a built-in tag",
		)))
}

fn genotype_text(genotype: Genotype, out: &mut Vec<u8>) {
	match genotype {
		Genotype::Diploid(first, second) => {
			out.extend_from_slice(format!("{first}/{second}").as_bytes());
		}
		Genotype::Haploid(allele) => out.extend_from_slice(allele.to_string().as_bytes()),
		Genotype::Missing(gyrase_call::Ploidy::Diploid) => out.extend_from_slice(b"./."),
		Genotype::Missing(gyrase_call::Ploidy::Haploid) => out.push(b'.'),
	}
}

fn saturate(value: u32) -> i32 {
	i32::try_from(value).unwrap_or(i32::MAX)
}

#[allow(clippy::cast_possible_truncation)]
fn single(value: f64) -> f32 {
	value as f32
}
