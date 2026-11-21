use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use gyrase_extract::{
	CallBuild, CallThresholds, Format, Panel, Settings, Sex, SexCall, Summary, Timestamp,
	write_raw_data,
};
use gyrase_genome::{Build, Reference};
use gyrase_io::Workers;

use crate::align::{self, Alignment};
use crate::alignments::{header, indexed};
use crate::bam::{coordinate_sorted, index, sort, sort_budget};
use crate::partial::Partial;
use crate::{
	Alignments, Calling, Cancel, Error, Estimate, Input, Kind, Lane, Progress, Resources, Sample,
	Source, Stage, Store, Targets, Vcf, shards,
};

const IO_THREADS: NonZeroUsize = NonZeroUsize::MIN.saturating_add(1);
const THREAD_MEMORY: u64 = 512 << 20;
const DEFAULT_BUILD: Build = Build::Grch38;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Request {
	pub output: PathBuf,
	pub formats: Vec<Format>,
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
	pub variants: PathBuf,
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
		let (build, steps, alignments) = match &sample.source {
			Source::Alignments(input) => {
				let header = header(input, None, &Workers::new(IO_THREADS)?)?;
				let build = Build::detect(&header.contigs()?).unwrap_or(DEFAULT_BUILD);
				match (indexed(input), input.kind, coordinate_sorted(&header)) {
					(true, _, _) => (build, Vec::new(), None),
					(false, Kind::Bam, true) => (build, vec![Stage::Index], None),
					(false, Kind::Bam, false) => {
						(build, vec![Stage::Sort], Some(file("sorted.bam")))
					}
					(false, _, _) => return Err(Error::MissingIndex(input.path.clone())),
				}
			}
			Source::Reads(_) => (
				DEFAULT_BUILD,
				vec![Stage::Align, Stage::Sort],
				Some(file("bam")),
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
			variants: file("vcf.gz"),
			raw: request
				.formats
				.iter()
				.filter(|_| exported)
				.map(|&format| {
					(
						format,
						file(&format!("{}.{}", format.name(), format.extension())),
					)
				})
				.collect(),
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
		let pools = Pools::new(resources)?;
		let reference = Arc::new(Reference::open(&self.reference(
			store,
			cancel,
			&mut report,
		)?)?);
		let prepared;
		let input = match (&self.sample.source, &self.outputs.alignments) {
			(Source::Alignments(input), None) => {
				if self.stages.contains(&Stage::Index) {
					index(input, &pools.compute, cancel, &mut report)?;
				}
				input
			}
			(Source::Alignments(input), Some(path)) => {
				let budget = sort_budget(resources.memory);
				sort(input, path, budget, &pools.compute, cancel, &mut report)?;
				prepared = Input::open(path)?;
				&prepared
			}
			(Source::Reads(lanes), Some(path)) => {
				self.align(
					lanes,
					path,
					&reference,
					store,
					resources,
					&pools,
					cancel,
					&mut report,
				)?;
				prepared = Input::open(path)?;
				&prepared
			}
			_ => return Err(Error::NotAlignments(self.sample.path().to_owned())),
		};
		let alignments = Alignments::open(input, reference, &pools.io)?;
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
		let mut vcf = Vcf::create(&self.outputs.variants, &alignments, &pools.io)?;
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
		if let Some(build) = call_build {
			let calls = targets.site_calls(&panel, &genotyped);
			let mut estimate =
				Estimate::new(Stage::Export, count(self.outputs.raw.len()), Instant::now());
			for (index, (format, path)) in self.outputs.raw.iter().enumerate() {
				cancel.check()?;
				let settings = Settings {
					format: *format,
					sex: match survey.sex {
						SexCall::Determined(sex) => sex,
						SexCall::Undetermined => Sex::Female,
					},
					thresholds: CallThresholds::DEFAULT,
					build,
					generated: now()?,
				};
				let (partial, file) = Partial::create(path)?;
				raw.push((*format, write_raw_data(&panel, &calls, &settings, file)?));
				partial.commit()?;
				report(estimate.advance(count(index + 1), Instant::now()));
			}
		}
		Ok(Report {
			sex: survey.sex,
			variants,
			raw,
		})
	}

	#[allow(clippy::too_many_arguments)]
	fn align(
		&self,
		lanes: &[Lane],
		path: &Path,
		reference: &Reference,
		store: &Store,
		resources: &Resources,
		pools: &Pools,
		cancel: &Cancel,
		report: &mut impl FnMut(Progress),
	) -> Result<(), Error> {
		let preparing = self.stages.contains(&Stage::Reference);
		let mut estimate = Estimate::new(Stage::Reference, 1, Instant::now());
		if preparing {
			report(estimate.progress());
		}
		let index = align::index(
			reference,
			&store.aligner_index_path(self.build),
			&pools.compute,
		)?;
		if preparing {
			report(estimate.advance(1, Instant::now()));
		}
		let resident = u64::try_from(index.heap_bytes()).unwrap_or(u64::MAX);
		Alignment {
			index: &index,
			compute: &pools.compute,
			io: &pools.io,
			cancel,
		}
		.run(
			lanes,
			&self.sample.name,
			path,
			sort_budget(resources.memory.saturating_sub(resident)),
			report,
		)?;
		Ok(())
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
