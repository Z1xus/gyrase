use std::io::{self, Write};

use crate::{Allele, Site, SiteName, Timestamp};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Resolved {
	NoCall,
	One(Allele),
	Two(Allele, Allele),
}

pub(super) trait Layout {
	fn header(&self, sink: &mut impl Write, generated: Timestamp) -> io::Result<()>;
	fn sort_key(&self, site: &Site) -> (u8, u32);
	fn row(
		&self,
		sink: &mut impl Write,
		name: SiteName<'_>,
		site: &Site,
		genotype: Resolved,
	) -> io::Result<()>;
}

pub(super) fn ascii_ordered(first: Allele, second: Allele) -> [u8; 2] {
	let (first, second) = (first.ascii(), second.ascii());
	if first <= second {
		[first, second]
	} else {
		[second, first]
	}
}

pub(super) fn doubled_or_ordered(genotype: Resolved) -> Option<[u8; 2]> {
	match genotype {
		Resolved::NoCall => None,
		Resolved::One(allele) => Some([allele.ascii(); 2]),
		Resolved::Two(first, second) => Some(ascii_ordered(first, second)),
	}
}
