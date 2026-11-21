use std::collections::HashSet;
use std::ffi::OsString;
use std::io::{IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Duration;

use gyrase_extract::{Format, Sex, SexCall};
use gyrase_jobs::{
	Cancel, Input, Plan, Progress, Report, Request, Resources, Source, Stage, Store, samples,
};
use lexopt::prelude::*;

use crate::failure::Failure;

const USAGE: &str = "\
Usage: gyr <file>... [options]
       gyr <samtools|bcftools> <command> [options]

Options:
  -o, --output DIR     write files to DIR [next to the input]
  -f, --formats LIST   23andme,ancestry,ftdna,myheritage,livingdna or none [all]
  -b, --bam            also write a BAM file from FASTQ input
  -h, --help           show this help
";

pub(crate) fn main(args: Vec<OsString>) -> ExitCode {
	match run(args) {
		Ok(()) => ExitCode::SUCCESS,
		Err(Failure::Help) => {
			let _ = write!(std::io::stdout(), "{USAGE}");
			ExitCode::SUCCESS
		}
		Err(Failure::Failed(code)) => ExitCode::from(code),
		Err(Failure::Usage(message)) => {
			let _ = write!(std::io::stderr(), "gyr: {message}\n{USAGE}");
			ExitCode::FAILURE
		}
		Err(Failure::Error(message)) => {
			let _ = writeln!(std::io::stderr(), "gyr: {message}");
			ExitCode::FAILURE
		}
	}
}

fn run(args: Vec<OsString>) -> Result<(), Failure> {
	let mut parser = lexopt::Parser::from_args(args);
	let mut output = None;
	let mut formats = Format::ALL.to_vec();
	let mut bam = false;
	let mut gpu = true;
	let mut paths = Vec::new();
	while let Some(arg) = parser.next()? {
		match arg {
			Short('o') | Long("output") => output = Some(PathBuf::from(parser.value()?)),
			Short('f') | Long("formats") => formats = parse_formats(&parser.value()?.string()?)?,
			Short('b') | Long("bam") => bam = true,
			Long("no-gpu") => gpu = false,
			Short('h') | Long("help") => return Err(Failure::Help),
			Value(value) => paths.push(PathBuf::from(value)),
			_ => return Err(arg.unexpected().into()),
		}
	}
	if paths.is_empty() {
		return Err(Failure::usage("no input files"));
	}
	let inputs = paths
		.iter()
		.map(|path| Input::open(path))
		.collect::<Result<Vec<_>, _>>()?;
	let store = Store::locate()?;
	let resources = Resources {
		gpu,
		..Resources::detect()?
	};
	let cancel = Cancel::default();
	let handler = cancel.clone();
	ctrlc::set_handler(move || handler.cancel())
		.map_err(|error| Failure::new(error.to_string()))?;
	let plans = samples(inputs)?
		.into_iter()
		.map(|sample| {
			let folder = output.clone().unwrap_or_else(|| {
				sample
					.path()
					.parent()
					.filter(|parent| !parent.as_os_str().is_empty())
					.map_or_else(|| PathBuf::from("."), Path::to_owned)
			});
			Plan::new(
				sample,
				&store,
				&Request {
					output: folder,
					formats: formats.clone(),
					bam,
				},
			)
		})
		.collect::<Result<Vec<_>, _>>()?;
	let mut used: HashSet<PathBuf> = paths.iter().map(|path| identity(path)).collect();
	if let Some(path) = plans
		.iter()
		.flat_map(written)
		.find(|path| !used.insert(identity(path)))
	{
		return Err(Failure::new(format!(
			"{}: would overwrite an input or another output",
			path.display()
		)));
	}
	for plan in &plans {
		let mut line = Line::new();
		let report = plan.run(&store, &resources, &cancel, |progress| line.show(progress));
		line.clear();
		summarize(plan, &report?)?;
	}
	Ok(())
}

fn written(plan: &Plan) -> impl Iterator<Item = PathBuf> {
	let outputs = &plan.outputs;
	let alignments = outputs
		.alignments
		.iter()
		.flat_map(|path| [path.clone(), sidecar(path, ".bai"), sidecar(path, ".csi")]);
	let variants = outputs
		.variants
		.iter()
		.flat_map(|path| [path.clone(), sidecar(path, ".tbi")]);
	alignments
		.chain(variants)
		.chain(outputs.raw.iter().map(|(_, path)| path.clone()))
}

fn sidecar(path: &Path, suffix: &str) -> PathBuf {
	let mut name = path.as_os_str().to_owned();
	name.push(suffix);
	PathBuf::from(name)
}

fn identity(path: &Path) -> PathBuf {
	path.canonicalize().unwrap_or_else(|_| {
		let parent = path
			.parent()
			.filter(|parent| !parent.as_os_str().is_empty())
			.unwrap_or(Path::new("."));
		match (parent.canonicalize(), path.file_name()) {
			(Ok(folder), Some(name)) => folder.join(name),
			_ => std::path::absolute(path).unwrap_or_else(|_| path.to_owned()),
		}
	})
}

fn parse_formats(list: &str) -> Result<Vec<Format>, Failure> {
	if list == "none" {
		return Ok(Vec::new());
	}
	list.split(',')
		.map(|name| {
			Format::ALL
				.into_iter()
				.find(|format| format.name() == name.trim())
				.ok_or_else(|| Failure::usage(format!("unknown format '{name}'")))
		})
		.collect()
}

struct Line {
	terminal: bool,
	shown: Option<(Stage, u64)>,
}

impl Line {
	fn new() -> Self {
		Self {
			terminal: std::io::stderr().is_terminal(),
			shown: None,
		}
	}

	fn show(&mut self, progress: Progress) {
		let percent = progress.done.saturating_mul(100) / progress.total.max(1);
		let key = (progress.stage, if self.terminal { percent } else { 0 });
		if self.shown == Some(key) {
			return;
		}
		self.shown = Some(key);
		let eta = progress
			.eta
			.map_or_default(|eta| format!("  {} left", duration(eta)));
		let text = format!("{:<20} {percent:>3}%{eta}", stage(progress.stage));
		let mut stderr = std::io::stderr();
		let _ = if self.terminal {
			write!(stderr, "\r\x1b[2K{text}")
		} else {
			writeln!(stderr, "{}", stage(progress.stage))
		};
	}

	fn clear(&self) {
		if self.terminal && self.shown.is_some() {
			let _ = write!(std::io::stderr(), "\r\x1b[2K");
		}
	}
}

fn summarize(plan: &Plan, report: &Report) -> Result<(), Failure> {
	let mut out = std::io::stdout().lock();
	let sex = match report.sex {
		SexCall::Determined(Sex::Male) => "  male",
		SexCall::Determined(Sex::Female) => "  female",
		SexCall::Undetermined => "",
	};
	writeln!(out, "{}  {}{sex}", plan.sample.name, plan.build.name())?;
	let prepared = match plan.sample.source {
		Source::Reads(_) => "aligned",
		_ => "sorted",
	};
	let rows = plan
		.outputs
		.alignments
		.iter()
		.map(|path| (path, String::from(prepared)))
		.chain(
			plan.outputs
				.variants
				.iter()
				.map(|path| (path, format!("{} variants", grouped(report.variants)))),
		)
		.chain(
			plan.outputs
				.raw
				.iter()
				.zip(&report.raw)
				.map(|((_, path), (_, summary))| {
					(
						path,
						format!(
							"{} of {} sites",
							grouped(summary.called),
							grouped(summary.sites)
						),
					)
				}),
		);
	let rows: Vec<(String, String)> = rows
		.map(|(path, detail)| (path.display().to_string(), detail))
		.collect();
	let width = rows.iter().map(|(path, _)| path.len()).max().unwrap_or(0);
	for (path, detail) in rows {
		writeln!(out, "  {path:<width$}  {detail}")?;
	}
	Ok(())
}

fn stage(stage: Stage) -> &'static str {
	match stage {
		Stage::Reference => "Preparing reference",
		Stage::Index => "Indexing",
		Stage::Align => "Aligning",
		Stage::Sort => "Sorting",
		Stage::Call => "Calling variants",
		Stage::Export => "Writing raw data",
	}
}

fn duration(duration: Duration) -> String {
	let seconds = duration.as_secs();
	match seconds {
		0..60 => format!("{seconds}s"),
		60..3_600 => format!("{}m {:02}s", seconds / 60, seconds % 60),
		_ => format!("{}h {:02}m", seconds / 3_600, seconds / 60 % 60),
	}
}

fn grouped(value: usize) -> String {
	let digits = value.to_string();
	let mut text = String::new();
	for (index, digit) in digits.chars().enumerate() {
		if index > 0 && (digits.len() - index).is_multiple_of(3) {
			text.push(',');
		}
		text.push(digit);
	}
	text
}
