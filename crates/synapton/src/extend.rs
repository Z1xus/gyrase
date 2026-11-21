use std::ops::Range;

use crate::alignment::CigarOp;
use crate::num::{u16_of, u32_of};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExtensionJob {
	pub(crate) query: Range<u32>,
	pub(crate) window: Range<u32>,
	pub(crate) band: u16,
}

impl ExtensionJob {
	pub(crate) fn band(&self) -> usize {
		usize::from(self.band)
	}
}

#[derive(Debug, Clone, Default)]
pub(crate) struct ExtensionBatch {
	queries: Vec<u8>,
	windows: Vec<u8>,
	jobs: Vec<ExtensionJob>,
}

impl ExtensionBatch {
	pub(crate) fn jobs(&self) -> &[ExtensionJob] {
		&self.jobs
	}

	pub(crate) fn query(&self, job: &ExtensionJob) -> &[u8] {
		self.queries
			.get(job.query.start as usize..job.query.end as usize)
			.unwrap_or_default()
	}

	pub(crate) fn window(&self, job: &ExtensionJob) -> &[u8] {
		self.windows
			.get(job.window.start as usize..job.window.end as usize)
			.unwrap_or_default()
	}

	pub(crate) fn clear(&mut self) {
		self.queries.clear();
		self.windows.clear();
		self.jobs.clear();
	}

	pub(crate) fn len(&self) -> usize {
		self.jobs.len()
	}

	pub(crate) fn push(&mut self, query: &[u8], window: &[u8], band: usize) {
		let query_start = u32_of(self.queries.len());
		self.queries.extend_from_slice(query);
		let window_start = u32_of(self.windows.len());
		self.windows.extend_from_slice(window);
		self.jobs.push(ExtensionJob {
			query: query_start..u32_of(self.queries.len()),
			window: window_start..u32_of(self.windows.len()),
			band: u16_of(band),
		});
	}
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Extension {
	pub(crate) score: i32,
	pub(crate) query_start: u32,
	pub(crate) query_end: u32,
	pub(crate) window_start: u32,
	pub(crate) window_end: u32,
	pub(crate) operations: Range<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Outcome {
	Unaligned,
	Aligned(Extension),
}

#[derive(Debug, Clone, Default)]
pub(crate) struct ExtensionResults {
	outcomes: Vec<Outcome>,
	operations: Vec<CigarOp>,
}

impl ExtensionResults {
	pub(crate) fn outcomes(&self) -> &[Outcome] {
		&self.outcomes
	}

	pub(crate) fn operations(&self, extension: &Extension) -> &[CigarOp] {
		self.operations
			.get(extension.operations.start as usize..extension.operations.end as usize)
			.unwrap_or_default()
	}

	pub(crate) fn reset(&mut self, jobs: usize) {
		self.outcomes.clear();
		self.outcomes.resize(jobs, Outcome::Unaligned);
		self.operations.clear();
	}

	pub(crate) fn set(&mut self, job: usize, outcome: Outcome) {
		if let Some(slot) = self.outcomes.get_mut(job) {
			*slot = outcome;
		}
	}

	pub(crate) fn push_operations(&mut self, operations: &[CigarOp]) -> Range<u32> {
		let start = u32_of(self.operations.len());
		self.operations.extend_from_slice(operations);
		start..u32_of(self.operations.len())
	}
}

impl Extension {
	pub(crate) fn new(
		score: i32,
		query: Range<u32>,
		window: Range<u32>,
		operations: Range<u32>,
	) -> Self {
		Self {
			score,
			query_start: query.start,
			query_end: query.end,
			window_start: window.start,
			window_end: window.end,
			operations,
		}
	}
}
