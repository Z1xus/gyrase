use std::io::Write;

use gyrase_genome::{Base, Chromosome};

use super::bytes::{put_varint, zigzag};
use super::{MAX_NAME_LEN, Panel};

use crate::site::NameCode;
use crate::{Chip, Error, Orientation};

pub(super) const MAGIC: [u8; 8] = *b"GYRPANEL";
pub(super) const VERSION: u16 = 1;
pub(super) const ORIENTATION_BIT: u8 = 1 << 6;
pub(super) const ALLELE_SHIFT: u32 = 3;
pub(super) const RESERVED_BITS: u8 = 1 << 7;

pub(super) const fn base_code(base: Option<Base>) -> u8 {
	match base {
		Some(Base::A) => 1,
		Some(Base::C) => 2,
		Some(Base::G) => 3,
		Some(Base::T) => 4,
		Some(Base::N) | None => 0,
	}
}

impl Panel {
	pub fn write(&self, mut sink: impl Write) -> Result<(), Error> {
		let mut head = Vec::new();
		head.extend_from_slice(&MAGIC);
		head.extend_from_slice(&VERSION.to_le_bytes());
		let chip_count = u16::try_from(Chip::ALL.len()).unwrap_or(0);
		head.extend_from_slice(&chip_count.to_le_bytes());
		let site_count = u32::try_from(self.sites.len()).map_err(|_| Error::Corrupt("sites"))?;
		head.extend_from_slice(&site_count.to_le_bytes());
		let name_count = u32::try_from(self.names.len()).map_err(|_| Error::Corrupt("names"))?;
		head.extend_from_slice(&name_count.to_le_bytes());
		for name in &self.names {
			let len = u8::try_from(name.len()).map_err(|_| Error::InvalidName(name.to_string()))?;
			if usize::from(len) > MAX_NAME_LEN {
				return Err(Error::InvalidName(name.to_string()));
			}
			head.push(len);
			head.extend_from_slice(name.as_bytes());
		}
		let mut counts = [0u32; Chromosome::COUNT];
		for site in &self.sites {
			if let Some(count) = counts.get_mut(site.chromosome.index()) {
				*count += 1;
			}
		}
		for count in counts {
			head.extend_from_slice(&count.to_le_bytes());
		}
		sink.write_all(&head)?;
		for section in self.columns()? {
			let len = u32::try_from(section.len()).map_err(|_| Error::Corrupt("section"))?;
			sink.write_all(&len.to_le_bytes())?;
			sink.write_all(&section)?;
		}
		Ok(())
	}

	fn columns(&self) -> Result<[Vec<u8>; 5], Error> {
		let mut positions = Vec::new();
		let mut offsets = Vec::new();
		let mut names = Vec::new();
		let mut alleles = Vec::with_capacity(self.sites.len());
		let mut chips = Vec::with_capacity(self.sites.len() * 2);
		let mut previous_chromosome = None;
		let mut previous_position = 0u32;
		let mut previous_offset = 0i64;
		for site in &self.sites {
			if previous_chromosome != Some(site.chromosome) {
				previous_chromosome = Some(site.chromosome);
				previous_position = 0;
			}
			let delta = site
				.grch38
				.checked_sub(previous_position)
				.ok_or(Error::Corrupt("unsorted sites"))?;
			put_varint(&mut positions, u64::from(delta));
			previous_position = site.grch38;
			let offset = i64::from(site.grch37) - i64::from(site.grch38);
			put_varint(&mut offsets, zigzag(offset - previous_offset));
			previous_offset = offset;
			let code = match site.name {
				NameCode::Rs(number) => u64::from(number) << 1,
				NameCode::Named(index) => (u64::from(index) << 1) | 1,
			};
			put_varint(&mut names, code);
			let orientation = match site.orientation {
				Orientation::Same => 0,
				Orientation::Opposite => ORIENTATION_BIT,
			};
			alleles.push(
				base_code(site.reference)
					| (base_code(site.alternate) << ALLELE_SHIFT)
					| orientation,
			);
			chips.extend_from_slice(&site.chips.bits().to_le_bytes());
		}
		Ok([positions, offsets, names, alleles, chips])
	}
}
