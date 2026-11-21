use std::io::{Cursor, Read};

use gyrase_genome::{Base, Chromosome};
use gyrase_io::{Workers, decode};

use super::bytes::{Bytes, unzigzag};
use super::write::{ALLELE_SHIFT, MAGIC, ORIENTATION_BIT, RESERVED_BITS, VERSION};
use super::{MAX_SITES, Panel};

use crate::site::NameCode;
use crate::{Chip, ChipSet, Error, Orientation, Site};

const EMBEDDED: &[u8] = include_bytes!("../../data/panel.bgz");
const MAX_PANEL_BYTES: u64 = 256 << 20;
const MAX_NAMES: usize = 1 << 24;
const ALLELE_MASK: u8 = 0b111;

const fn base_from_code(code: u8) -> Result<Option<Base>, Error> {
	match code {
		0 => Ok(None),
		1 => Ok(Some(Base::A)),
		2 => Ok(Some(Base::C)),
		3 => Ok(Some(Base::G)),
		4 => Ok(Some(Base::T)),
		_ => Err(Error::Corrupt("allele code")),
	}
}

impl Panel {
	pub fn embedded(workers: &Workers) -> Result<Self, Error> {
		Self::load(EMBEDDED, workers)
	}

	pub fn load(compressed: &'static [u8], workers: &Workers) -> Result<Self, Error> {
		Self::read(decode(Box::new(Cursor::new(compressed)), workers)?)
	}

	pub fn read(source: impl Read) -> Result<Self, Error> {
		let mut data = Vec::new();
		source.take(MAX_PANEL_BYTES + 1).read_to_end(&mut data)?;
		if data.len() as u64 > MAX_PANEL_BYTES {
			return Err(Error::Corrupt("panel too large"));
		}
		Self::parse(Bytes::new(&data))
	}

	fn parse(mut bytes: Bytes<'_>) -> Result<Self, Error> {
		if bytes.array::<8>().map_err(|_| Error::NotPanel)? != MAGIC {
			return Err(Error::NotPanel);
		}
		let version = bytes.u16()?;
		if version != VERSION {
			return Err(Error::UnsupportedVersion(version));
		}
		if usize::from(bytes.u16()?) != Chip::ALL.len() {
			return Err(Error::Corrupt("chip count"));
		}
		let site_count = bytes.length()?;
		if site_count > MAX_SITES {
			return Err(Error::TooManySites(MAX_SITES));
		}
		let name_count = bytes.length()?;
		if name_count > MAX_NAMES {
			return Err(Error::Corrupt("name count"));
		}
		let mut names = Vec::new();
		for _ in 0..name_count {
			let len = usize::from(bytes.u8()?);
			let name = std::str::from_utf8(bytes.take(len)?)
				.map_err(|_| Error::Corrupt("name encoding"))?;
			names.push(Box::from(name));
		}
		let mut counts = [0usize; Chromosome::COUNT];
		let mut total = 0usize;
		for count in &mut counts {
			*count = bytes.length()?;
			total = total
				.checked_add(*count)
				.ok_or(Error::Corrupt("chromosome counts"))?;
		}
		if total != site_count {
			return Err(Error::Corrupt("chromosome counts"));
		}
		let mut positions = bytes.section()?;
		let mut offsets = bytes.section()?;
		let mut name_codes = bytes.section()?;
		let mut alleles = bytes.section()?;
		let mut chips = bytes.section()?;
		if !bytes.is_empty() {
			return Err(Error::Corrupt("trailing bytes"));
		}
		if alleles.len() != site_count || chips.len() != site_count * 2 {
			return Err(Error::Corrupt("column length"));
		}
		let mut sites = Vec::with_capacity(site_count);
		let mut offset = 0i64;
		for (chromosome_index, count) in counts.into_iter().enumerate() {
			let chromosome =
				Chromosome::from_index(chromosome_index).ok_or(Error::Corrupt("chromosome"))?;
			let mut grch38 = 0u32;
			for _ in 0..count {
				let delta = u32::try_from(positions.varint()?)
					.map_err(|_| Error::Corrupt("position delta"))?;
				grch38 = grch38
					.checked_add(delta)
					.ok_or(Error::Corrupt("position"))?;
				offset = offset
					.checked_add(unzigzag(offsets.varint()?))
					.ok_or(Error::Corrupt("offset"))?;
				let grch37 = u32::try_from(i64::from(grch38) + offset)
					.map_err(|_| Error::Corrupt("grch37 position"))?;
				let name = decode_name(name_codes.varint()?, names.len())?;
				let packed = alleles.u8()?;
				if packed & RESERVED_BITS != 0 {
					return Err(Error::Corrupt("allele flags"));
				}
				let orientation = if packed & ORIENTATION_BIT == 0 {
					Orientation::Same
				} else {
					Orientation::Opposite
				};
				sites.push(Site {
					chromosome,
					grch37,
					grch38,
					name,
					reference: base_from_code(packed & ALLELE_MASK)?,
					alternate: base_from_code((packed >> ALLELE_SHIFT) & ALLELE_MASK)?,
					orientation,
					chips: ChipSet::from_bits(chips.u16()?).ok_or(Error::Corrupt("chip bits"))?,
				});
			}
		}
		if !(positions.is_empty() && offsets.is_empty() && name_codes.is_empty()) {
			return Err(Error::Corrupt("column length"));
		}
		Ok(Self { sites, names })
	}
}

fn decode_name(code: u64, name_count: usize) -> Result<NameCode, Error> {
	let value = u32::try_from(code >> 1).map_err(|_| Error::Corrupt("name code"))?;
	if code & 1 == 0 {
		return Ok(NameCode::Rs(value));
	}
	if usize::try_from(value).is_ok_and(|index| index < name_count) {
		Ok(NameCode::Named(value))
	} else {
		Err(Error::Corrupt("name index"))
	}
}
