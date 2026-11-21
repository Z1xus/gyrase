use std::path::Path;
use std::time::Duration;

use gyrase_extract::Format;
use gyrase_jobs::{Kind, Stage};

use crate::StepKind;

const SIZE_UNITS: [(u64, &str); 3] = [(1_000_000_000, "GB"), (1_000_000, "MB"), (1_000, "KB")];
const MINUTE: u64 = 60;
const HOUR: u64 = 3600;

pub const fn format(format: Format) -> &'static str {
	match format {
		Format::TwentyThreeAndMeV5 => "23andMe",
		Format::AncestryV2 => "Ancestry",
		Format::FamilyFinder => "FTDNA",
		Format::MyHeritage => "MyHeritage",
		Format::LivingDna => "LivingDNA",
	}
}

pub const fn kind(kind: Kind) -> &'static str {
	match kind {
		Kind::Fastq => "FASTQ",
		Kind::Fasta => "FASTA",
		Kind::Sam => "SAM",
		Kind::Bam => "BAM",
		Kind::Cram => "CRAM",
		Kind::Vcf => "VCF",
		Kind::Bcf => "BCF",
	}
}

pub const fn step_kind(stage: Stage) -> StepKind {
	match stage {
		Stage::Reference => StepKind::Reference,
		Stage::Index => StepKind::Index,
		Stage::Align => StepKind::Align,
		Stage::Sort => StepKind::Sort,
		Stage::Call => StepKind::Call,
		Stage::Export => StepKind::Export,
	}
}

pub const fn step(stage: Stage) -> &'static str {
	match stage {
		Stage::Reference => "Download reference",
		Stage::Index => "Index file",
		Stage::Align => "Align reads",
		Stage::Sort => "Sort reads",
		Stage::Call => "Find variants",
		Stage::Export => "Write raw data",
	}
}

pub const fn activity(stage: Stage) -> &'static str {
	match stage {
		Stage::Reference => "Downloading reference…",
		Stage::Index => "Indexing…",
		Stage::Align => "Aligning reads…",
		Stage::Sort => "Sorting reads…",
		Stage::Call => "Finding variants…",
		Stage::Export => "Writing raw data…",
	}
}

pub fn size(bytes: u64) -> String {
	SIZE_UNITS
		.iter()
		.find(|(unit, _)| bytes >= *unit)
		.map_or_else(
			|| format!("{bytes} B"),
			|(unit, name)| {
				let tenths = bytes / (unit / 10);
				format!("{}.{} {name}", tenths / 10, tenths % 10)
			},
		)
}

pub fn eta(remaining: Duration) -> String {
	let seconds = remaining.as_secs();
	if seconds < MINUTE {
		format!("{seconds} s left")
	} else if seconds < HOUR {
		format!("{} min left", seconds.div_ceil(MINUTE))
	} else {
		format!("{} h {} min left", seconds / HOUR, seconds % HOUR / MINUTE)
	}
}

pub fn grouped(value: usize) -> String {
	let digits = value.to_string();
	let mut grouped = String::new();
	for (index, digit) in digits.chars().enumerate() {
		if index > 0 && (digits.len() - index).is_multiple_of(3) {
			grouped.push(',');
		}
		grouped.push(digit);
	}
	grouped
}

pub fn file_name(path: &Path) -> String {
	path.file_name()
		.map(|name| name.to_string_lossy().into_owned())
		.unwrap_or_default()
}
