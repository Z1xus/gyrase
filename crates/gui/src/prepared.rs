use std::path::{Path, PathBuf};

use gyrase_extract::Format;
use gyrase_jobs::{Input, Plan, Request, Sample, Store, samples};

use crate::message;

const SINGLE_SAMPLE: &str = "Add one sample at a time.";

pub(crate) struct Prepared {
	pub(crate) plan: Plan,
	pub(crate) formats: Vec<Format>,
	pub(crate) output: PathBuf,
}

impl Prepared {
	pub(crate) fn open(paths: &[PathBuf]) -> Result<Self, String> {
		let inputs = paths
			.iter()
			.map(|path| Input::open(path))
			.collect::<Result<Vec<_>, _>>()
			.map_err(|error| message::text(&error))?;
		let found = samples(inputs).map_err(|error| message::text(&error))?;
		let [sample] = <[Sample; 1]>::try_from(found).map_err(|_| SINGLE_SAMPLE.to_owned())?;
		let output = paths
			.first()
			.and_then(|path| path.parent())
			.filter(|folder| !folder.as_os_str().is_empty())
			.unwrap_or(Path::new("."));
		Self::plan(sample, Format::ALL.to_vec(), output.to_owned())
	}

	pub(crate) fn toggled(&self, format: Format) -> Vec<Format> {
		Format::ALL
			.into_iter()
			.filter(|&candidate| (candidate == format) ^ self.formats.contains(&candidate))
			.collect()
	}

	pub(crate) fn with(&self, formats: Vec<Format>, output: PathBuf) -> Result<Self, String> {
		Self::plan(self.plan.sample.clone(), formats, output)
	}

	pub(crate) fn replanned(self) -> Self {
		self.with(self.formats.clone(), self.output.clone())
			.unwrap_or(self)
	}

	fn plan(sample: Sample, formats: Vec<Format>, output: PathBuf) -> Result<Self, String> {
		let store = Store::locate().map_err(|error| message::text(&error))?;
		let request = Request {
			output: output.clone(),
			formats: formats.clone(),
		};
		let plan = Plan::new(sample, &store, &request).map_err(|error| message::text(&error))?;
		Ok(Self {
			plan,
			formats,
			output,
		})
	}
}
