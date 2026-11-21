mod ancestry;
mod ftdna;
mod layout;
mod livingdna;
mod myheritage;
mod resolve;
mod twentythree;

use std::io::{BufWriter, Write};

use gyrase_genome::Build;

use self::layout::{Layout, Resolved};
use crate::{CallThresholds, Chip, Error, Panel, Sex, SiteCall, Timestamp};

const WRITE_BUFFER_BYTES: usize = 1 << 20;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
	TwentyThreeAndMeV5,
	AncestryV2,
	FamilyFinder,
	MyHeritage,
	LivingDna,
}

impl Format {
	pub const fn chip(self) -> Chip {
		match self {
			Self::TwentyThreeAndMeV5 => Chip::TwentyThreeAndMeV5,
			Self::AncestryV2 => Chip::AncestryV2,
			Self::FamilyFinder => Chip::FtdnaV4,
			Self::MyHeritage => Chip::MyHeritageGsa,
			Self::LivingDna => Chip::LivingDnaV2,
		}
	}

	pub const fn extension(self) -> &'static str {
		match self {
			Self::TwentyThreeAndMeV5 | Self::AncestryV2 | Self::LivingDna => "txt",
			Self::FamilyFinder | Self::MyHeritage => "csv",
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallBuild {
	Grch37,
	Grch38,
}

impl TryFrom<Build> for CallBuild {
	type Error = Error;

	fn try_from(build: Build) -> Result<Self, Error> {
		match build {
			Build::Grch37 => Ok(Self::Grch37),
			Build::Grch38 => Ok(Self::Grch38),
			Build::Chm13 => Err(Error::UnsupportedBuild),
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Settings {
	pub format: Format,
	pub sex: Sex,
	pub thresholds: CallThresholds,
	pub build: CallBuild,
	pub generated: Timestamp,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Summary {
	pub sites: usize,
	pub called: usize,
}

pub fn write_raw_data(
	panel: &Panel,
	calls: &[SiteCall],
	settings: &Settings,
	sink: impl Write,
) -> Result<Summary, Error> {
	if calls.len() != panel.len() {
		return Err(Error::CallCount {
			expected: panel.len(),
			actual: calls.len(),
		});
	}
	match settings.format {
		Format::TwentyThreeAndMeV5 => {
			export(&twentythree::TwentyThree, panel, calls, settings, sink)
		}
		Format::AncestryV2 => export(&ancestry::Ancestry, panel, calls, settings, sink),
		Format::FamilyFinder => export(&ftdna::FamilyFinder, panel, calls, settings, sink),
		Format::MyHeritage => export(&myheritage::MyHeritage, panel, calls, settings, sink),
		Format::LivingDna => export(&livingdna::LivingDna, panel, calls, settings, sink),
	}
}

fn export<L: Layout>(
	layout: &L,
	panel: &Panel,
	calls: &[SiteCall],
	settings: &Settings,
	sink: impl Write,
) -> Result<Summary, Error> {
	let mut order: Vec<usize> = panel
		.chip_sites(settings.format.chip())
		.map(|(index, _)| index)
		.collect();
	order.sort_by_key(|&index| {
		panel
			.sites()
			.get(index)
			.map_or((u8::MAX, 0, index), |site| {
				let (rank, position) = layout.sort_key(site);
				(rank, position, index)
			})
	});
	let mut sink = BufWriter::with_capacity(WRITE_BUFFER_BYTES, sink);
	layout.header(&mut sink, settings.generated)?;
	let mut called = 0;
	for &index in &order {
		let (Some(site), Some(call)) = (panel.sites().get(index), calls.get(index)) else {
			continue;
		};
		let genotype = resolve::resolve(call, site, settings);
		called += usize::from(!matches!(genotype, Resolved::NoCall));
		layout.row(&mut sink, panel.name(site), site, genotype)?;
	}
	sink.flush()?;
	Ok(Summary {
		sites: order.len(),
		called,
	})
}
