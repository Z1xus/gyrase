use std::io::{self, Write};

use super::layout::{Layout, Resolved, doubled_or_ordered};
use crate::{Site, SiteName, Timestamp};

pub(super) struct MyHeritage;

const HEADER_AFTER_TIMESTAMP: &str = " UTC\r\n##reference=build37\r\n#\r\n# MyHeritage DNA raw data.\r\n# For each SNP,  we provide the identifier, chromosome number, base pair position and genotype.\r\n# The genotype is reported on the forward (+) strand with respect to the human reference build 37.\r\n# THIS INFORMATION IS FOR YOUR PERSONAL USE AND IS INTENDED FOR GENEALOGICAL RESEARCH\r\n# ONLY. IT IS NOT INTENDED FOR MEDICAL OR HEALTH PURPOSES. PLEASE BE AWARE THAT THE\r\n# DOWNLOADED DATA WILL NO LONGER BE PROTECTED BY OUR SECURITY MEASURES.\r\nRSID,CHROMOSOME,POSITION,RESULT\r\n";

impl Layout for MyHeritage {
	fn header(&self, sink: &mut impl Write, generated: Timestamp) -> io::Result<()> {
		let stamp = generated.iso_utc();
		let stamp = stamp.strip_suffix(" UTC").unwrap_or(&stamp);
		write!(
			sink,
			"##fileformat=MyHeritage\r\n##format=MHv1.0\r\n##chip=GSA\r\n##timestamp={stamp}{HEADER_AFTER_TIMESTAMP}"
		)
	}

	fn sort_key(&self, site: &Site) -> (u8, u32) {
		(
			u8::try_from(site.chrom().index()).unwrap_or(u8::MAX),
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
		write!(
			sink,
			"\"{name}\",\"{}\",\"{}\",\"",
			site.chrom().label(),
			site.grch37()
		)?;
		sink.write_all(&doubled_or_ordered(genotype).unwrap_or(*b"--"))?;
		sink.write_all(b"\"\r\n")
	}
}
