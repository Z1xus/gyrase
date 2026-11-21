use std::fs::File;
use std::path::Path;

use gyrase_genome::Interval;
use gyrase_index::Index;
use gyrase_io::Workers;
use gyrase_variants::{ChunkReader, Header, Reader, Regions, open_regions};

use crate::failure::{Failure, file};
use crate::shared::open;

pub enum Source {
	Whole(Reader),
	Ranged(Reader<ChunkReader>),
}

macro_rules! each {
	($source:expr, $reader:ident => $body:expr) => {
		match $source {
			$crate::bcftools::input::Source::Whole($reader) => $body,
			$crate::bcftools::input::Source::Ranged($reader) => $body,
		}
	};
}
pub(crate) use each;

impl Source {
	pub fn header(&self) -> &Header {
		match self {
			Self::Whole(reader) => reader.header(),
			Self::Ranged(reader) => reader.header(),
		}
	}
}

#[derive(Default)]
pub struct Selection {
	regions: Vec<String>,
	targets: Vec<String>,
}

impl Selection {
	pub fn regions(&mut self, list: &str) {
		self.regions.extend(list.split(',').map(single_position));
	}

	pub fn regions_file(&mut self, path: &str) -> Result<(), Failure> {
		self.regions.extend(read_regions(path)?);
		Ok(())
	}

	pub fn targets(&mut self, list: &str) {
		self.targets.extend(list.split(',').map(single_position));
	}

	pub fn targets_file(&mut self, path: &str) -> Result<(), Failure> {
		self.targets.extend(read_regions(path)?);
		Ok(())
	}

	pub fn open(
		&self,
		path: &str,
		workers: &Workers,
	) -> Result<(Source, Option<Regions>), Failure> {
		if !self.regions.is_empty() && !self.targets.is_empty() {
			return Err(Failure::new("-r or -R cannot be combined with -t or -T"));
		}
		if self.regions.is_empty() {
			let reader = Reader::open(open(path)?, workers)?;
			let regions = if self.targets.is_empty() {
				None
			} else {
				let contigs = reader.header().contigs()?;
				let intervals = self
					.targets
					.iter()
					.map(|text| Interval::parse(text, &contigs))
					.collect::<Result<Vec<_>, _>>()?;
				Some(Regions::positions(&intervals))
			};
			return Ok((Source::Whole(reader), regions));
		}
		let index_path = ["tbi", "csi"]
			.iter()
			.map(|extension| format!("{path}.{extension}"))
			.find(|candidate| Path::new(candidate).exists())
			.ok_or_else(|| Failure::new(format!("failed to load the index of \"{path}\"")))?;
		let index = Index::read(File::open(&index_path)?, workers)?;
		let texts: Vec<&str> = self.regions.iter().map(String::as_str).collect();
		let (reader, regions) = open_regions(Path::new(path), &index, &texts, workers)?;
		Ok((Source::Ranged(reader), Some(regions)))
	}
}

fn single_position(region: &str) -> String {
	match region.rsplit_once(':') {
		Some((_, tail)) if !tail.is_empty() && tail.bytes().all(|byte| byte.is_ascii_digit()) => {
			format!("{region}-{tail}")
		}
		_ => region.to_owned(),
	}
}

fn read_regions(path: &str) -> Result<Vec<String>, Failure> {
	let text = std::fs::read_to_string(path).map_err(|error| file(path, "reading", &error))?;
	let bed = Path::new(path)
		.extension()
		.is_some_and(|extension| extension.eq_ignore_ascii_case("bed"));
	let mut regions = Vec::new();
	for line in text.lines() {
		if line.is_empty() || line.starts_with(['#', '@']) || line.starts_with("track") {
			continue;
		}
		let mut fields = line.split('\t');
		let (Some(name), Some(first)) = (fields.next(), fields.next()) else {
			return Err(Failure::new(format!(
				"could not parse the line \"{line}\" in \"{path}\""
			)));
		};
		let last = fields.next().filter(|last| !last.is_empty());
		let parse = |text: &str| {
			text.parse::<u64>().map_err(|_| {
				Failure::new(format!("could not parse the line \"{line}\" in \"{path}\""))
			})
		};
		let (first, last) = if bed {
			let end = parse(last.unwrap_or(first))?;
			(parse(first)? + 1, end)
		} else {
			let first = parse(first)?;
			(first, last.map_or(Ok(first), parse)?)
		};
		if last >= first {
			regions.push(format!("{{{name}}}:{first}-{last}"));
		}
	}
	Ok(regions)
}
