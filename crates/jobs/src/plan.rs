use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use gyrase_extract::{
	CallBuild, CallThresholds, Format, Panel, Settings, Sex, SexCall, SiteCall, Summary, Timestamp,
	write_raw_data,
};
use gyrase_genome::{Build, Reference};
use gyrase_io::Workers;
use gyrase_reads::Header;

use crate::align::{self, Alignment};
use crate::alignments::{header, indexed};
use crate::bam::{Output, Sorting, coordinate_sorted, hold_budget, index, sort_budget};
use crate::genotypes;
use crate::partial::{Partial, suffixed};
use crate::{
	Alignments, Calling, Cancel, Error, Estimate, Input, Kind, Lane, Progress, Resources, Sample,
	Source, Stage, Store, Targets, Vcf, shards,
};

const IO_THREADS: NonZeroUsize = NonZeroUsize::MIN.saturating_add(1);
const THREAD_MEMORY: u64 = 512 << 20;
const ALIGNED_BUILD: Build = Build::Grch38;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Request {
	pub output: PathBuf,
	pub formats: Vec<Format>,
	pub bam: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Plan {
	pub sample: Sample,
	pub build: Build,
	pub stages: Vec<Stage>,
	pub outputs: Outputs,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Outputs {
	pub alignments: Option<PathBuf>,
	pub variants: Option<PathBuf>,
	pub raw: Vec<(Format, PathBuf)>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Report {
	pub sex: SexCall,
	pub variants: usize,
	pub raw: Vec<(Format, Summary)>,
}

struct Pools {
	compute: Workers,
	io: Workers,
}

impl Plan {
	pub fn new(sample: Sample, store: &Store, request: &Request) -> Result<Self, Error> {
		let file = |suffix: &str| request.output.join(format!("{}.{suffix}", sample.name));
		let raw = |exported: bool| -> Vec<(Format, PathBuf)> {
			request
				.formats
				.iter()
				.filter(|_| exported)
				.map(|&format| {
					(
						format,
						file(&format!("{}.{}", format.name(), format.extension())),
					)
				})
				.collect()
		};
		if let Source::Variants(input) = &sample.source {
			if request.formats.is_empty() {
				return Err(Error::NoFormats(input.path.clone()));
			}
			let build = genotypes::build(input, &Workers::new(IO_THREADS)?)?;
			CallBuild::try_from(build)?;
			return Ok(Self {
				build,
				stages: vec![Stage::Export],
				outputs: Outputs {
					alignments: None,
					variants: None,
					raw: raw(true),
				},
				sample,
			});
		}
		let (build, steps, alignments) = match &sample.source {
			Source::Alignments(input) => {
				let header = header(input, None, &Workers::new(IO_THREADS)?)?;
				let build = Build::detect(&header.contigs()?)
					.ok_or_else(|| Error::UnknownBuild(input.path.clone()))?;
				match (input.kind, indexed(input), coordinate_sorted(&header)) {
					(Kind::Bam | Kind::Cram, true, _) => (build, Vec::new(), None),
					(Kind::Bam | Kind::Cram, false, true) => (build, vec![Stage::Index], None),
					_ => (build, vec![Stage::Sort], Some(file("sorted.bam"))),
				}
			}
			Source::Reads(_) => (
				ALIGNED_BUILD,
				vec![Stage::Align, Stage::Sort],
				request.bam.then(|| file("bam")),
			),
			Source::Variants(_) | Source::Reference(_) => {
				return Err(Error::NotAlignments(sample.path().to_owned()));
			}
		};
		let exported = CallBuild::try_from(build).is_ok();
		let aligned = matches!(sample.source, Source::Reads(_));
		let mut stages = Vec::new();
		if store.reference(build).is_none()
			|| (aligned && !store.aligner_index_path(build).is_file())
		{
			stages.push(Stage::Reference);
		}
		stages.extend(steps);
		stages.push(Stage::Call);
		if exported && !request.formats.is_empty() {
			stages.push(Stage::Export);
		}
		let outputs = Outputs {
			alignments,
			variants: Some(file("vcf.gz")),
			raw: raw(exported),
		};
		Ok(Self {
			sample,
			build,
			stages,
			outputs,
		})
	}

	pub fn run(
		&self,
		store: &Store,
		resources: &Resources,
		cancel: &Cancel,
		mut report: impl FnMut(Progress),
	) -> Result<Report, Error> {
		self.create_folder()?;
		let pools = Pools::new(resources)?;
		if let Source::Variants(input) = &self.sample.source {
			return self.convert(input, &pools, cancel, &mut report);
		}
		let reference = Arc::new(Reference::open(&self.reference(
			store,
			cancel,
			&mut report,
		)?)?);
		let prepared;
		let mut _scratch = None;
		let alignments = match (&self.sample.source, &self.outputs.alignments) {
			(Source::Alignments(input), None) => {
				if self.stages.contains(&Stage::Index) {
					index(input, &pools.compute, cancel, &mut report)?;
				}
				Alignments::open(input, reference, &pools.io)?
			}
			(Source::Alignments(input), Some(path)) => {
				Sorting {
					reference: &reference,
					budget: sort_budget(resources.memory),
					workers: &pools.compute,
					cancel,
				}
				.run(input, path, &mut report)?;
				prepared = Input::open(path)?;
				Alignments::open(&prepared, reference, &pools.io)?
			}
			(Source::Reads(lanes), path) => {
				let (found, guard) = self.align_reads(
					lanes,
					path.as_deref(),
					reference,
					store,
					resources,
					&pools,
					cancel,
					&mut report,
				)?;
				_scratch = guard;
				found
			}
			_ => return Err(Error::NotAlignments(self.sample.path().to_owned())),
		};
		let call_build = CallBuild::try_from(self.build).ok();
		let panel = Panel::embedded(&pools.io)?;
		let targets = call_build.map_or_else(Targets::default, |build| {
			Targets::new(&panel, &alignments, build)
		});
		let shards = shards(&alignments, self.build);
		let calling = Calling {
			alignments: &alignments,
			workers: &pools.compute,
			io: &pools.io,
			cancel,
		};
		let total = shards.iter().map(|shard| shard.cost).sum();
		let mut estimate = Estimate::new(Stage::Call, total, Instant::now());
		report(estimate.progress());
		let survey = calling.survey(&shards)?;
		let variants = self
			.outputs
			.variants
			.as_ref()
			.ok_or_else(|| Error::NotAlignments(self.sample.path().to_owned()))?;
		let mut vcf = Vcf::create(variants, &alignments, &pools.io)?;
		let mut done = 0;
		let mut variants = 0;
		let genotyped = calling.call(&shards, &survey, &targets.targets, |shard, calls| {
			done += shard.cost;
			variants += calls.len();
			report(estimate.advance(done, Instant::now()));
			vcf.write(shard.interval.contig.get(), calls)
		})?;
		vcf.finish(&pools.io)?;
		let mut raw = Vec::new();
		if call_build.is_some() {
			let calls = targets.site_calls(&panel, &genotyped);
			let mut estimate =
				Estimate::new(Stage::Export, count(self.outputs.raw.len()), Instant::now());
			raw = self.write_raw(&panel, &calls, survey.sex, cancel, |done| {
				report(estimate.advance(done, Instant::now()));
			})?;
		}
		Ok(Report {
			sex: survey.sex,
			variants,
			raw,
		})
	}

	fn create_folder(&self) -> Result<(), Error> {
		let first = self
			.outputs
			.variants
			.iter()
			.chain(self.outputs.raw.iter().map(|(_, path)| path))
			.next();
		if let Some(folder) = first.and_then(|path| path.parent()) {
			std::fs::create_dir_all(folder).map_err(|source| Error::Write {
				path: folder.to_owned(),
				source,
			})?;
		}
		Ok(())
	}

	fn convert(
		&self,
		input: &Input,
		pools: &Pools,
		cancel: &Cancel,
		report: &mut impl FnMut(Progress),
	) -> Result<Report, Error> {
		let panel = Panel::embedded(&pools.io)?;
		let formats = count(self.outputs.raw.len());
		let mut estimate = Estimate::new(Stage::Export, input.size + formats, Instant::now());
		report(estimate.progress());
		let genotypes = genotypes::read(input, &panel, &pools.compute, |bytes| {
			cancel.check()?;
			report(estimate.advance(bytes, Instant::now()));
			Ok(())
		})?;
		let raw = self.write_raw(&panel, &genotypes.calls, genotypes.sex, cancel, |done| {
			report(estimate.advance(input.size + done, Instant::now()));
		})?;
		Ok(Report {
			sex: genotypes.sex,
			variants: genotypes.records,
			raw,
		})
	}

	fn write_raw(
		&self,
		panel: &Panel,
		calls: &[SiteCall],
		sex: SexCall,
		cancel: &Cancel,
		mut advance: impl FnMut(u64),
	) -> Result<Vec<(Format, Summary)>, Error> {
		let build = CallBuild::try_from(self.build)?;
		let mut raw = Vec::new();
		for (index, (format, path)) in self.outputs.raw.iter().enumerate() {
			cancel.check()?;
			let settings = Settings {
				format: *format,
				sex: match sex {
					SexCall::Determined(sex) => sex,
					SexCall::Undetermined => Sex::Female,
				},
				thresholds: CallThresholds::DEFAULT,
				build,
				generated: now()?,
			};
			let (partial, file) = Partial::create(path)?;
			raw.push((*format, write_raw_data(panel, calls, &settings, file)?));
			partial.commit()?;
			advance(count(index + 1));
		}
		Ok(raw)
	}

	#[allow(clippy::too_many_arguments)]
	fn align_reads(
		&self,
		lanes: &[Lane],
		path: Option<&Path>,
		reference: Arc<Reference>,
		store: &Store,
		resources: &Resources,
		pools: &Pools,
		cancel: &Cancel,
		report: &mut impl FnMut(Progress),
	) -> Result<(Alignments, Option<Scratch>), Error> {
		let output = path.map_or_else(|| self.scratch(), Path::to_owned);
		let guard = path.is_none().then(|| Scratch(output.clone()));
		let (header, held) = self.align(
			lanes,
			&output,
			path.is_none(),
			&reference,
			store,
			resources,
			pools,
			cancel,
			report,
		)?;
		match held {
			Output::Memory(records) => Ok((Alignments::held(header, reference, records), None)),
			Output::File => {
				let alignments = Alignments::open(&Input::open(&output)?, reference, &pools.io)?;
				Ok((alignments, guard))
			}
		}
	}

	#[allow(clippy::too_many_arguments)]
	fn align(
		&self,
		lanes: &[Lane],
		path: &Path,
		hold: bool,
		reference: &Reference,
		store: &Store,
		resources: &Resources,
		pools: &Pools,
		cancel: &Cancel,
		report: &mut impl FnMut(Progress),
	) -> Result<(Header, Output), Error> {
		let preparing = self.stages.contains(&Stage::Reference);
		let mut estimate = Estimate::new(Stage::Reference, 1, Instant::now());
		let index = align::index(
			reference,
			&store.aligner_index_path(self.build),
			&pools.compute,
			cancel,
		)?;
		if preparing {
			report(estimate.advance(1, Instant::now()));
		}
		let resident = u64::try_from(index.heap_bytes()).unwrap_or(u64::MAX);
		let memory = resources.memory.saturating_sub(resident);
		Alignment {
			index: &index,
			compute: &pools.compute,
			gpu: resources.gpu,
			io: &pools.io,
			cancel,
		}
		.run(
			lanes,
			&self.sample.name,
			path,
			hold,
			if hold {
				hold_budget(memory)
			} else {
				sort_budget(memory)
			},
			report,
		)
	}

	fn scratch(&self) -> PathBuf {
		let folder = self
			.outputs
			.variants
			.as_deref()
			.and_then(Path::parent)
			.unwrap_or(Path::new("."));
		folder.join(format!(
			"{}.{}.tmp.bam",
			self.sample.name,
			std::process::id()
		))
	}

	fn reference(
		&self,
		store: &Store,
		cancel: &Cancel,
		report: &mut impl FnMut(Progress),
	) -> Result<PathBuf, Error> {
		if let Some(path) = store.reference(self.build) {
			return Ok(path);
		}
		let mut estimate: Option<Estimate> = None;
		store.fetch(self.build, cancel, |download| {
			let now = Instant::now();
			let estimate = estimate.get_or_insert_with(|| {
				Estimate::new(
					Stage::Reference,
					download.total.unwrap_or(download.received),
					now,
				)
			});
			report(estimate.advance(download.received, now));
		})
	}
}

struct Scratch(PathBuf);

impl Drop for Scratch {
	fn drop(&mut self) {
		for path in [
			self.0.clone(),
			suffixed(&self.0, ".bai"),
			suffixed(&self.0, ".csi"),
		] {
			let _ = std::fs::remove_file(path);
		}
	}
}

impl Pools {
	fn new(resources: &Resources) -> Result<Self, Error> {
		let fitting = usize::try_from(resources.memory / THREAD_MEMORY).unwrap_or(usize::MAX);
		let threads =
			NonZeroUsize::new(resources.threads.get().min(fitting)).unwrap_or(NonZeroUsize::MIN);
		Ok(Self {
			compute: Workers::new(threads)?,
			io: Workers::new(IO_THREADS)?,
		})
	}
}

fn now() -> Result<Timestamp, Error> {
	let seconds = SystemTime::now()
		.duration_since(UNIX_EPOCH)
		.map_or(0, |elapsed| elapsed.as_secs());
	Ok(Timestamp::from_unix(seconds)?)
}

fn count(value: usize) -> u64 {
	u64::try_from(value).unwrap_or(u64::MAX)
}
