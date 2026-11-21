mod bytes;
mod read;
mod write;

use std::collections::{HashMap, HashSet};

use gyrase_genome::Base;

use crate::site::{NameCode, rs_number};
use crate::{Chip, Error, Site, SiteEntry, SiteName};

const MAX_SITES: usize = 8_000_000;
const MAX_NAME_LEN: usize = 255;

#[derive(Debug, Clone)]
pub struct Panel {
	sites: Vec<Site>,
	names: Vec<Box<str>>,
}

impl Panel {
	pub fn build(entries: impl IntoIterator<Item = SiteEntry>) -> Result<Self, Error> {
		let mut entries: Vec<SiteEntry> = entries.into_iter().collect();
		if entries.len() > MAX_SITES {
			return Err(Error::TooManySites(MAX_SITES));
		}
		entries.sort_by(|left, right| {
			let left_key = (left.chrom, left.grch38.unwrap_or(0), left.grch37);
			let right_key = (right.chrom, right.grch38.unwrap_or(0), right.grch37);
			left_key
				.cmp(&right_key)
				.then_with(|| left.name.cmp(&right.name))
		});
		let mut seen = HashSet::with_capacity(entries.len());
		let mut interned: HashMap<String, u32> = HashMap::new();
		let mut names: Vec<Box<str>> = Vec::new();
		let mut sites = Vec::with_capacity(entries.len());
		for entry in &entries {
			if entry.name.is_empty() || entry.name.len() > MAX_NAME_LEN {
				return Err(Error::InvalidName(entry.name.clone()));
			}
			if !seen.insert((entry.chrom, entry.grch37, entry.name.as_str())) {
				return Err(Error::DuplicateSite {
					name: entry.name.clone(),
					chrom: entry.chrom.label(),
					position: entry.grch37,
				});
			}
			let name = match rs_number(&entry.name) {
				Some(number) => NameCode::Rs(number),
				None => NameCode::Named(intern(&mut interned, &mut names, &entry.name)?),
			};
			sites.push(Site {
				chrom: entry.chrom,
				grch37: entry.grch37,
				grch38: entry.grch38.unwrap_or(0),
				name,
				reference: entry.reference.filter(|base| *base != Base::N),
				alternate: entry.alternate.filter(|base| *base != Base::N),
				orientation: entry.orientation,
				chips: entry.chips,
			});
		}
		Ok(Self { sites, names })
	}

	pub fn sites(&self) -> &[Site] {
		&self.sites
	}

	pub fn len(&self) -> usize {
		self.sites.len()
	}

	pub fn is_empty(&self) -> bool {
		self.sites.is_empty()
	}

	pub fn name(&self, site: &Site) -> SiteName<'_> {
		match site.name {
			NameCode::Rs(number) => SiteName::Rs(number),
			NameCode::Named(index) => SiteName::Named(
				usize::try_from(index)
					.ok()
					.and_then(|index| self.names.get(index))
					.map_or("", |name| &**name),
			),
		}
	}

	pub fn chip_sites(&self, chip: Chip) -> impl Iterator<Item = (usize, &Site)> {
		self.sites
			.iter()
			.enumerate()
			.filter(move |(_, site)| site.chips.contains(chip))
	}
}

fn intern(
	interned: &mut HashMap<String, u32>,
	names: &mut Vec<Box<str>>,
	name: &str,
) -> Result<u32, Error> {
	if let Some(&index) = interned.get(name) {
		return Ok(index);
	}
	let index = u32::try_from(names.len()).map_err(|_| Error::TooManySites(MAX_SITES))?;
	interned.insert(name.to_owned(), index);
	names.push(name.into());
	Ok(index)
}
