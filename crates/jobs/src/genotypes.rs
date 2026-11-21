use std::fs::File;
use std::ops::Range;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use gyrase_extract::{
	Allele, CallBuild, CallThresholds, Genotype, Orientation, Panel, Ploidy, Sex, SexCall,
	SexEvidence, Site as PanelSite, SiteCall, expected_ploidy, in_x_par, infer_sex,
};
use gyrase_genome::{Base, Build, Chromosome};
use gyrase_io::Workers;
use gyrase_variants::{Alleles, END_INT, Header, LineBatch, Loaded, Reader, RecordBatch, Site};

use crate::counted::Counted;
use crate::{Error, Input};

const DEFAULT_BUILD: Build = Build::Grch38;
const BATCH_BYTES: usize = 1 << 22;
const MIN_Y_CALLS: usize = 10;
const PLOIDY_LIMIT: usize = 2;

#[derive(Debug, Clone, PartialEq)]
pub struct Genotypes {
	pub build: Build,
	pub calls: Vec<SiteCall>,
	pub sex: SexCall,
	pub records: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Source {
	None,
	Fill,
	Block,
	Record,
}

#[derive(Debug, Clone, Copy)]
enum Verdict {
	Skip,
	Missing,
	Call(SiteCall, Option<Allele>),
}

impl Verdict {
	fn of(call: Option<SiteCall>, reference: Option<Allele>) -> Self {
		call.map_or(Self::Missing, |call| Self::Call(call, reference))
	}
}

fn combine(kept: SiteCall, added: SiteCall, reference: Option<Allele>) -> SiteCall {
	let pick = |current: Allele, other: Allele| {
		if Some(current) == reference {
			other
		} else {
			current
		}
	};
	let genotype = match (kept.genotype, added.genotype) {
		(Genotype::Pair(first, second), Genotype::Pair(other, then)) => {
			Genotype::Pair(pick(first, other), pick(second, then))
		}
		(Genotype::Hemizygous(only), Genotype::Hemizygous(other)) => {
			Genotype::Hemizygous(pick(only, other))
		}
		(genotype, _) => genotype,
	};
	SiteCall { genotype, ..kept }
}

#[derive(Debug, Clone, Copy)]
struct Genotyped {
	alleles: [Option<usize>; PLOIDY_LIMIT],
	ploidy: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct Support {
	depth: Option<u32>,
	quality: Option<u32>,
}

#[derive(Debug, Clone, Copy)]
struct Keys {
	genotype: Option<u32>,
	format_depth: Option<u32>,
	info_depth: Option<u32>,
	minimum_depth: Option<u32>,
	allele_depths: Option<u32>,
	quality: Option<u32>,
}

#[derive(Debug, Default)]
struct Scratch {
	alleles: Alleles,
	loaded: Loaded,
	filters: Vec<u32>,
}

#[derive(Debug)]
struct Collector<'a> {
	panel: &'a Panel,
	build: CallBuild,
	keys: Keys,
	table: Vec<(u64, usize)>,
	calls: Vec<SiteCall>,
	sources: Vec<Source>,
	seen: [bool; Chromosome::COUNT],
	gvcf: bool,
	records: usize,
}

fn open(input: &Input, workers: &Workers, bytes: &Arc<AtomicU64>) -> Result<Reader, Error> {
	let file = File::open(&input.path).map_err(|source| Error::Open {
		path: input.path.clone(),
		source,
	})?;
	Ok(Reader::open(Box::new(Counted::new(file, bytes)), workers)?)
}

fn detect(header: &Header, path: &Path) -> Result<Build, Error> {
	if header.contig_count() == 0 {
		return Ok(DEFAULT_BUILD);
	}
	Build::detect(&header.contigs()?).ok_or_else(|| Error::UnknownBuild(path.to_owned()))
}

pub fn build(input: &Input, workers: &Workers) -> Result<Build, Error> {
	let reader = open(input, workers, &Arc::default())?;
	detect(reader.header(), &input.path)
}

pub fn read(
	input: &Input,
	panel: &Panel,
	workers: &Workers,
	mut progress: impl FnMut(u64) -> Result<(), Error>,
) -> Result<Genotypes, Error> {
	let bytes = Arc::new(AtomicU64::new(0));
	let reader = open(input, workers, &bytes)?;
	let header = reader.header().clone();
	if header.sample_count() == 0 {
		return Err(Error::NoSamples(input.path.clone()));
	}
	let build = detect(&header, &input.path)?;
	let mut collector = Collector::new(panel, CallBuild::try_from(build)?, &header);
	let mut scratch = Scratch::default();
	match reader {
		Reader::Vcf(mut vcf) => {
			let mut batch = LineBatch::default();
			let mut names = Names::new(&header);
			while vcf.read_batch(&mut batch, BATCH_BYTES)? {
				for line in batch.lines() {
					let line = line?;
					let chromosome = names.resolve(line.chrom());
					collector.record(&header, &line, chromosome, &mut scratch)?;
				}
				progress(bytes.load(Ordering::Relaxed))?;
			}
		}
		Reader::Bcf(mut bcf) => {
			let mut batch = RecordBatch::default();
			let chromosomes: Vec<Option<Chromosome>> =
				header.contig_names().map(Chromosome::parse).collect();
			while bcf.read_batch(&mut batch, BATCH_BYTES)? {
				for record in batch.records() {
					let record = record?;
					let chromosome = usize::try_from(record.contig()?)
						.ok()
						.and_then(|index| chromosomes.get(index).copied().flatten());
					collector.record(&header, &record, chromosome, &mut scratch)?;
				}
				progress(bytes.load(Ordering::Relaxed))?;
			}
		}
	}
	progress(bytes.load(Ordering::Relaxed))?;
	Ok(collector.finish(build))
}

#[derive(Debug)]
struct Names {
	known: Vec<(Vec<u8>, Option<Chromosome>)>,
	last: usize,
}

impl Names {
	fn new(header: &Header) -> Self {
		Self {
			known: header
				.contig_names()
				.map(|name| (name.as_bytes().to_vec(), Chromosome::parse(name)))
				.collect(),
			last: 0,
		}
	}

	fn resolve(&mut self, name: &[u8]) -> Option<Chromosome> {
		if let Some((known, chromosome)) = self.known.get(self.last)
			&& known == name
		{
			return *chromosome;
		}
		if let Some(index) = self.known.iter().position(|(known, _)| known == name) {
			self.last = index;
			return self
				.known
				.get(index)
				.and_then(|(_, chromosome)| *chromosome);
		}
		let chromosome = std::str::from_utf8(name).ok().and_then(Chromosome::parse);
		self.known.push((name.to_vec(), chromosome));
		self.last = self.known.len() - 1;
		chromosome
	}
}

fn key(chromosome: Chromosome, position: u32) -> u64 {
	(u64::try_from(chromosome.index()).unwrap_or_default() << 32) | u64::from(position)
}

fn count(value: usize) -> f64 {
	f64::from(u32::try_from(value).unwrap_or(u32::MAX))
}

fn unsigned(value: i32) -> Option<u32> {
	u32::try_from(value).ok()
}

fn allele_of(byte: u8) -> Option<Allele> {
	Allele::from_base(Base::from_ascii(byte))
}

fn symbolic(allele: &[u8]) -> bool {
	matches!(allele.first(), Some(b'<' | b'*')) || allele.is_empty()
}

fn reference_only(alleles: &Alleles) -> bool {
	alleles
		.iter()
		.skip(1)
		.all(|alt| matches!(alt, b"<NON_REF>" | b"<*>" | b"."))
}

fn genotyped(values: &[i32]) -> Option<Genotyped> {
	let mut alleles = [None; PLOIDY_LIMIT];
	let mut ploidy = 0;
	for &value in values {
		if value == END_INT {
			break;
		}
		let code = value >> 1;
		let allele = if code > 0 {
			usize::try_from(code - 1).ok()
		} else {
			None
		};
		*alleles.get_mut(ploidy)? = allele;
		ploidy += 1;
	}
	(ploidy > 0).then_some(Genotyped { alleles, ploidy })
}

fn make_call(
	genotype: Genotyped,
	allele: impl Fn(usize) -> Option<Allele>,
	support: Support,
) -> Option<SiteCall> {
	let [first, second] = genotype.alleles;
	let genotype = match genotype.ploidy {
		1 => Genotype::Hemizygous(allele(first?)?),
		_ => Genotype::Pair(allele(first?)?, allele(second?)?),
	};
	Some(SiteCall {
		genotype,
		depth: support.depth.unwrap_or(CallThresholds::DEFAULT.min_depth),
		quality: support
			.quality
			.unwrap_or(CallThresholds::DEFAULT.min_quality),
	})
}

impl Keys {
	fn new(header: &Header) -> Self {
		let key = |name: &str| header.key(name);
		let depth = key("DP");
		Self {
			genotype: key("GT").filter(|&key| header.format_def(key).is_some()),
			format_depth: depth.filter(|&key| header.format_def(key).is_some()),
			info_depth: depth.filter(|&key| header.info_def(key).is_some()),
			minimum_depth: key("MIN_DP").filter(|&key| header.format_def(key).is_some()),
			allele_depths: key("AD").filter(|&key| header.format_def(key).is_some()),
			quality: key("GQ").filter(|&key| header.format_def(key).is_some()),
		}
	}
}

fn first_sample<'a, S: Site>(
	site: &S,
	header: &Header,
	key: Option<u32>,
	loaded: &'a mut Loaded,
) -> Result<Option<&'a [i32]>, Error> {
	let Some(key) = key else {
		return Ok(None);
	};
	Ok(site
		.format(header, key, loaded)?
		.then(|| loaded.sample_ints(0)))
}

impl<'a> Collector<'a> {
	fn new(panel: &'a Panel, build: CallBuild, header: &Header) -> Self {
		let mut table: Vec<(u64, usize)> = panel
			.sites()
			.iter()
			.enumerate()
			.filter_map(|(index, site)| {
				let position = match build {
					CallBuild::Grch37 => Some(site.grch37()),
					CallBuild::Grch38 => site.grch38(),
				}?;
				Some((key(site.chromosome(), position), index))
			})
			.collect();
		table.sort_unstable();
		let gvcf = header.meta_lines().any(|line| {
			line.contains("NON_REF") || line.starts_with("##GVCFBlock") || line.contains("ID=*,")
		});
		Self {
			panel,
			build,
			keys: Keys::new(header),
			table,
			calls: vec![SiteCall::MISSING; panel.len()],
			sources: vec![Source::None; panel.len()],
			seen: [false; Chromosome::COUNT],
			gvcf,
			records: 0,
		}
	}

	fn sites(&self, chromosome: Chromosome, start: u32, end: u32) -> Range<usize> {
		let first = key(chromosome, start.saturating_add(1));
		let last = key(chromosome, end.max(start.saturating_add(1)));
		let begin = self.table.partition_point(|&(at, _)| at < first);
		let finish = self.table.partition_point(|&(at, _)| at <= last);
		begin..finish.max(begin)
	}

	fn expected(&self, site: &PanelSite) -> (Option<Allele>, Option<Allele>) {
		let flip = self.build == CallBuild::Grch38 && site.orientation() == Orientation::Opposite;
		let orient = |base: Option<Base>| {
			base.and_then(Allele::from_base)
				.map(|allele| if flip { allele.complement() } else { allele })
		};
		(orient(site.reference()), orient(site.alternate()))
	}

	fn record<S: Site>(
		&mut self,
		header: &Header,
		site: &S,
		chromosome: Option<Chromosome>,
		scratch: &mut Scratch,
	) -> Result<(), Error> {
		self.records += 1;
		let Some(chromosome) = chromosome else {
			return Ok(());
		};
		if let Some(seen) = self.seen.get_mut(chromosome.index()) {
			*seen = true;
		}
		let start = site.start()?;
		let end = site.end(header)?;
		let sites = self.sites(chromosome, start, end);
		if sites.is_empty() && self.gvcf {
			return Ok(());
		}
		site.alleles(&mut scratch.alleles)?;
		let reference_length = scratch.alleles.get(0).map_or(0, <[u8]>::len);
		let block = reference_only(&scratch.alleles);
		let non_ref = scratch
			.alleles
			.iter()
			.skip(1)
			.any(|alt| matches!(alt, b"<NON_REF>" | b"<*>"));
		let span = usize::try_from(end.saturating_sub(start)).unwrap_or_default();
		self.gvcf |= non_ref || (block && span > reference_length);
		if sites.is_empty() {
			return Ok(());
		}
		let passed = passed(site, header, scratch)?;
		let genotype = first_sample(site, header, self.keys.genotype, &mut scratch.loaded)?
			.and_then(genotyped);
		if genotype.is_none() && !block {
			return Ok(());
		}
		let support = self.support(site, header, scratch, block)?;
		for slot in sites {
			let Some(&(_, index)) = self.table.get(slot) else {
				continue;
			};
			let Some(panel_site) = self.panel.sites().get(index) else {
				continue;
			};
			let offset = usize::try_from(
				self.position(panel_site)
					.saturating_sub(1)
					.saturating_sub(start),
			)
			.unwrap_or_default();
			let call = if block {
				self.block_call(
					panel_site,
					offset,
					passed,
					genotype,
					support,
					&scratch.alleles,
				)
			} else {
				self.variant_call(
					panel_site,
					offset,
					passed,
					genotype,
					support,
					&scratch.alleles,
				)
			};
			self.store(index, block, call);
		}
		Ok(())
	}

	fn position(&self, site: &PanelSite) -> u32 {
		match self.build {
			CallBuild::Grch37 => site.grch37(),
			CallBuild::Grch38 => site.grch38().unwrap_or_default(),
		}
	}

	fn store(&mut self, index: usize, block: bool, verdict: Verdict) {
		let (Some(source), Some(slot)) = (self.sources.get_mut(index), self.calls.get_mut(index))
		else {
			return;
		};
		let (call, reference) = match verdict {
			Verdict::Skip => return,
			Verdict::Missing => (SiteCall::MISSING, None),
			Verdict::Call(call, reference) => (call, reference),
		};
		if block {
			if *source == Source::None {
				*source = Source::Block;
				*slot = call;
			}
		} else if *source != Source::Record || *slot == SiteCall::MISSING {
			*source = Source::Record;
			*slot = call;
		} else if call != SiteCall::MISSING {
			*slot = combine(*slot, call, reference);
		}
	}

	fn support<S: Site>(
		&self,
		site: &S,
		header: &Header,
		scratch: &mut Scratch,
		block: bool,
	) -> Result<Support, Error> {
		let quality = first_sample(site, header, self.keys.quality, &mut scratch.loaded)?
			.and_then(|values| values.first().copied())
			.and_then(unsigned);
		let mut depth = None;
		if block {
			depth = first_sample(site, header, self.keys.minimum_depth, &mut scratch.loaded)?
				.and_then(|values| values.first().copied())
				.and_then(unsigned);
		}
		if depth.is_none() {
			depth = first_sample(site, header, self.keys.format_depth, &mut scratch.loaded)?
				.and_then(|values| values.first().copied())
				.and_then(unsigned);
		}
		if depth.is_none()
			&& let Some(key) = self.keys.info_depth
			&& site.info(header, key, &mut scratch.loaded)?
		{
			depth = scratch.loaded.ints().first().copied().and_then(unsigned);
		}
		if depth.is_none() {
			depth = first_sample(site, header, self.keys.allele_depths, &mut scratch.loaded)?
				.map(|values| values.iter().copied().filter_map(unsigned).sum());
		}
		Ok(Support { depth, quality })
	}

	fn block_call(
		&self,
		site: &PanelSite,
		offset: usize,
		passed: bool,
		genotype: Option<Genotyped>,
		support: Support,
		alleles: &Alleles,
	) -> Verdict {
		let (reference, alternate) = self.expected(site);
		let Some(reference) = reference else {
			return Verdict::Skip;
		};
		if offset == 0 {
			let base = alleles.get(0).and_then(|bytes| bytes.first()).copied();
			let known = base.and_then(allele_of);
			if known.is_some() && known != Some(reference) && known != alternate {
				return Verdict::Missing;
			}
		}
		if !passed {
			return Verdict::Missing;
		}
		let genotype = genotype.unwrap_or(Genotyped {
			alleles: [Some(0); PLOIDY_LIMIT],
			ploidy: PLOIDY_LIMIT,
		});
		Verdict::of(
			make_call(genotype, |index| (index == 0).then_some(reference), support),
			Some(reference),
		)
	}

	fn variant_call(
		&self,
		site: &PanelSite,
		offset: usize,
		passed: bool,
		genotype: Option<Genotyped>,
		support: Support,
		alleles: &Alleles,
	) -> Verdict {
		let Some(reference_bytes) = alleles.get(0) else {
			return Verdict::Skip;
		};
		let (reference, alternate) = self.expected(site);
		let Some(reference) = reference else {
			return indel_call(offset, passed, genotype, support, alleles);
		};
		if !passed {
			return Verdict::Missing;
		}
		let Some(genotype) = genotype else {
			return Verdict::Skip;
		};
		let base = reference_bytes.get(offset).copied().and_then(allele_of);
		if base.is_none() || (base != Some(reference) && base != alternate) {
			return Verdict::Missing;
		}
		let allele = |index: usize| {
			let bytes = alleles.get(index)?;
			if index == 0 {
				return base;
			}
			if symbolic(bytes) {
				return None;
			}
			if offset == 0 {
				bytes.first().copied().and_then(allele_of)
			} else if bytes.len() == reference_bytes.len() {
				bytes.get(offset).copied().and_then(allele_of)
			} else {
				None
			}
		};
		Verdict::of(make_call(genotype, allele, support), base)
	}

	fn finish(mut self, build: Build) -> Genotypes {
		let sex = self.sex();
		let assumed = match sex {
			SexCall::Determined(sex) => sex,
			SexCall::Undetermined => Sex::Female,
		};
		for position in 0..self.table.len() {
			let Some(&(_, index)) = self.table.get(position) else {
				continue;
			};
			let Some(site) = self.panel.sites().get(index) else {
				continue;
			};
			let ploidy = expected_ploidy(assumed, site.chromosome(), site.grch37());
			let (reference, _) = self.expected(site);
			let seen = self
				.seen
				.get(site.chromosome().index())
				.copied()
				.unwrap_or(false);
			let (Some(source), Some(call)) =
				(self.sources.get_mut(index), self.calls.get_mut(index))
			else {
				continue;
			};
			match (*source, reference) {
				(Source::None, Some(reference))
					if seen && !self.gvcf && ploidy != Ploidy::Absent =>
				{
					*source = Source::Fill;
					*call = SiteCall {
						genotype: Genotype::Pair(reference, reference),
						depth: CallThresholds::DEFAULT.min_depth,
						quality: CallThresholds::DEFAULT.min_quality,
					};
				}
				_ => {}
			}
			if matches!(*source, Source::Fill | Source::Block)
				&& ploidy == Ploidy::Haploid
				&& let Genotype::Pair(first, second) = call.genotype
				&& first == second
			{
				call.genotype = Genotype::Hemizygous(first);
			}
		}
		Genotypes {
			build,
			calls: self.calls,
			sex,
			records: self.records,
		}
	}

	fn sex(&self) -> SexCall {
		let (mut diploid, mut heterozygous, mut haploid, mut unseen, mut called_y) =
			(0usize, 0usize, 0usize, 0usize, 0usize);
		let x_seen = self
			.seen
			.get(Chromosome::X.index())
			.copied()
			.unwrap_or(false);
		for &(_, index) in &self.table {
			let (Some(site), Some(call), Some(source)) = (
				self.panel.sites().get(index),
				self.calls.get(index),
				self.sources.get(index),
			) else {
				continue;
			};
			let observed = matches!(source, Source::Block | Source::Record);
			match site.chromosome() {
				Chromosome::X if !in_x_par(site.grch37()) => match call.genotype {
					Genotype::Pair(first, second) if observed => {
						diploid += 1;
						heterozygous += usize::from(first != second);
					}
					Genotype::Hemizygous(_) if observed => haploid += 1,
					Genotype::Missing if *source == Source::None && x_seen && !self.gvcf => {
						unseen += 1;
					}
					_ => {}
				},
				Chromosome::Y => {
					called_y += usize::from(observed && call.genotype != Genotype::Missing);
				}
				_ => {}
			}
		}
		if haploid > diploid {
			return SexCall::Determined(Sex::Male);
		}
		let total = diploid + unseen;
		let fraction = if total == 0 {
			f64::NAN
		} else {
			count(heterozygous) / count(total)
		};
		let inferred = infer_sex(SexEvidence {
			x_depth_ratio: f64::NAN,
			y_depth_ratio: f64::NAN,
			x_panel_heterozygous_fraction: fraction,
		});
		match inferred {
			SexCall::Undetermined if called_y >= MIN_Y_CALLS => SexCall::Determined(Sex::Male),
			other => other,
		}
	}
}

fn passed<S: Site>(site: &S, header: &Header, scratch: &mut Scratch) -> Result<bool, Error> {
	site.filters(header, &mut scratch.filters)?;
	Ok(scratch
		.filters
		.iter()
		.all(|&key| header.name(key) == Some("PASS")))
}

fn indel_call(
	offset: usize,
	passed: bool,
	genotype: Option<Genotyped>,
	support: Support,
	alleles: &Alleles,
) -> Verdict {
	let Some(reference_length) = alleles.get(0).map(<[u8]>::len) else {
		return Verdict::Skip;
	};
	let kinds: Vec<Option<Allele>> = alleles
		.iter()
		.skip(1)
		.map(|alt| {
			if symbolic(alt) {
				None
			} else if alt.len() > reference_length {
				Some(Allele::Insertion)
			} else if alt.len() < reference_length {
				Some(Allele::Deletion)
			} else {
				None
			}
		})
		.collect();
	if offset != 0 || kinds.iter().all(Option::is_none) {
		return Verdict::Skip;
	}
	if !passed {
		return Verdict::Missing;
	}
	let Some(genotype) = genotype else {
		return Verdict::Skip;
	};
	let present = |kind: Allele| kinds.contains(&Some(kind));
	let reference = match (present(Allele::Insertion), present(Allele::Deletion)) {
		(true, false) => Some(Allele::Deletion),
		(false, true) => Some(Allele::Insertion),
		_ => None,
	};
	let allele = |index: usize| match index.checked_sub(1) {
		None => reference,
		Some(alt) => kinds.get(alt).copied().flatten(),
	};
	Verdict::of(make_call(genotype, allele, support), reference)
}
