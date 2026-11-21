use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use gyrase_extract::{
	CallBuild, CallThresholds, Format, Panel, Settings, Sex, SexCall, Summary, Timestamp,
	write_raw_data,
};
use gyrase_genome::{Build, Reference};
use gyrase_io::Workers;

use crate::alignments::{header, indexed};
use crate::bam::{coordinate_sorted, index, sort, sort_budget};
use crate::partial::Partial;
use crate::{
	Alignments, Calling, Cancel, Error, Estimate, Input, Kind, Progress, Resources, Sample, Source,
	Stage, Store, Targets, Vcf, shards,
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
		let io = Workers::new(IO_THREADS)?;
		let Source::Alignments(input) = &sample.source else {
			return Err(Error::NotAlignments(sample.path().to_owned()));
		};
		let header = header(input, None, &io)?;
		let build = Build::detect(&header.contigs()?).unwrap_or(DEFAULT_BUILD);
		let exported = CallBuild::try_from(build).is_ok();
		let mut stages = Vec::new();
		if store.reference(build).is_none() {
			stages.push(Stage::Reference);
		}
		let file = |suffix: &str| request.output.join(format!("{}.{suffix}", sample.name));
		let mut alignments = None;
		if !indexed(input) {
			if input.kind != Kind::Bam {
				return Err(Error::MissingIndex(input.path.clone()));
			}
			if coordinate_sorted(&header) {
				stages.push(Stage::Index);
			} else {
				stages.push(Stage::Sort);
				alignments = Some(file("sorted.bam"));
			}
		}
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
		let Source::Alignments(input) = &self.sample.source else {
			return Err(Error::NotAlignments(self.sample.path().to_owned()));
		};
		let pools = Pools::new(resources)?;
		let reference = Arc::new(Reference::open(&self.reference(
			store,
			cancel,
			&mut report,
		)?)?);
		let sorted;
		let input = match &self.outputs.alignments {
			Some(path) => {
				sort(
					input,
					path,
					sort_budget(resources.memory),
					&pools.compute,
					cancel,
					&mut report,
				)?;
				sorted = Input::open(path)?;
				&sorted
			}
			None => input,
		};
		if self.stages.contains(&Stage::Index) {
			index(input, &pools.compute, cancel, &mut report)?;
		}
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
