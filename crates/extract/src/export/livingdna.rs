use std::io::{self, Write};

use super::layout::{Layout, Resolved, doubled_or_ordered};
use crate::{Site, SiteName, Timestamp};

pub(super) struct LivingDna;

const HEADER_AFTER_DATE: &str = "\t\t\t\r\n# Genotype chip: Sirius\t\t\t\r\n# The content of this file is subject to updates and changes depending on the time of download.\t\t\t\r\n# This genotype data should be treated as personal information.\t\t\t\r\n# This genotype data is not suitable for clinical/medical research or diagnosis.\t\t\t\r\n# The user assumes all responsibility for the security of this file.\t\t\t\r\n# Please refer to the Living DNA Terms and Conditions on our website (www.livingdna.com) for more information.\t\t\t\r\n# Human Genome Reference Build 37 (GRCh37.p13).\t\t\t\r\n# Genotypes are presented on the forward strand.\t\t\t\r\n#\t\t\t\r\n# rsid\tchromosome\tposition\tgenotype\r\n";

impl Layout for LivingDna {
	fn header(&self, sink: &mut impl Write, generated: Timestamp) -> io::Result<()> {
		write!(
			sink,
			"# Living DNA customer genotype data download file version: 1.0.2\t\t\t\r\n# File creation date: {}{HEADER_AFTER_DATE}",
			generated.month_day_year_plain()
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
			"{name}\t{}\t{}\t",
			site.chrom().label(),
			site.grch37()
		)?;
		sink.write_all(&doubled_or_ordered(genotype).unwrap_or(*b"--"))?;
		sink.write_all(b"\r\n")
	}
}
