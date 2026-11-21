use std::io::{self, Write};

use gyrase_genome::Chromosome;

use super::layout::{Layout, Resolved, doubled_or_ordered};
use crate::{Site, SiteName, Timestamp, in_x_par};

pub(super) struct FamilyFinder;

impl Layout for FamilyFinder {
	fn header(&self, sink: &mut impl Write, _generated: Timestamp) -> io::Result<()> {
		sink.write_all(b"RSID,CHROMOSOME,POSITION,RESULT\n")
	}

	fn sort_key(&self, site: &Site) -> (u8, u32) {
		(
			u8::try_from(site.chromosome().index()).unwrap_or(u8::MAX),
			site.grch37(),
		)
	}

	fn row(
		&self,
		sink: &mut impl Write,
		name: SiteName<'_>,
		site: &Site,
		genotype: Resolved,
	) -> io::Result<()> {
		let label = if site.chromosome() == Chromosome::X && in_x_par(site.grch37()) {
			"XY"
		} else {
			site.chromosome().label()
		};
		write!(sink, "{name},{label},{},", site.grch37())?;
		sink.write_all(&doubled_or_ordered(genotype).unwrap_or(*b"--"))?;
		sink.write_all(b"\n")
	}
}
