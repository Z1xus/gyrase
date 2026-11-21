use crate::{ContigId, Contigs, Error};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Interval {
	pub contig: ContigId,
	pub start: u64,
	pub end: u64,
}

impl Interval {
	pub const fn len(&self) -> u64 {
		self.end.saturating_sub(self.start)
	}

	pub const fn is_empty(&self) -> bool {
		self.end <= self.start
	}

	pub const fn contains(&self, contig: ContigId, position: u64) -> bool {
		self.contig.get() == contig.get() && self.start <= position && position < self.end
	}

	pub fn whole(contig: ContigId, contigs: &Contigs) -> Option<Self> {
		Some(Self {
			contig,
			start: 0,
			end: contigs.get(contig)?.length,
		})
	}

	pub fn parse(text: &str, contigs: &Contigs) -> Result<Self, Error> {
		let invalid = || Error::InvalidRegion(text.to_owned());
		let trimmed = text.trim();
		if let Some(contig) = contigs.resolve(trimmed) {
			return Self::whole(contig, contigs).ok_or_else(invalid);
		}
		let (name, range) = match trimmed
			.strip_prefix('{')
			.and_then(|rest| rest.split_once('}'))
		{
			Some((name, "")) => {
				let contig = contigs
					.resolve(name)
					.ok_or_else(|| Error::UnknownContig(name.to_owned()))?;
				return Self::whole(contig, contigs).ok_or_else(invalid);
			}
			Some((name, rest)) => (name, rest.strip_prefix(':').ok_or_else(invalid)?),
			None => trimmed
				.rsplit_once(':')
				.ok_or_else(|| Error::UnknownContig(trimmed.to_owned()))?,
		};
		let contig = contigs
			.resolve(name)
			.ok_or_else(|| Error::UnknownContig(name.to_owned()))?;
		let length = contigs.get(contig).ok_or_else(invalid)?.length;
		let (first, last) = match range.split_once('-') {
			Some((first, last)) => (first, last),
			None => (range, ""),
		};
		let start = match first.trim() {
			"" if range.contains('-') => 0,
			first => parse_coordinate(first)
				.ok_or_else(invalid)?
				.saturating_sub(1),
		};
		let end = match last.trim() {
			"" => length,
			last => parse_coordinate(last).ok_or_else(invalid)?.min(length),
		};
		if start >= end {
			return Err(invalid());
		}
		Ok(Self { contig, start, end })
	}
}

fn parse_coordinate(text: &str) -> Option<u64> {
	let digits: String = text
		.trim()
		.chars()
		.filter(|character| *character != ',')
		.collect();
	digits.parse().ok().filter(|&value| value > 0)
}
