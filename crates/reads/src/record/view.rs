use gyrase_genome::ContigId;

use super::cigar::{Cigar, CigarKind};
use super::fixed::{self, CORE_LEN, Core, Layout};
use super::flags::Flags;
use super::sequence::{Qualities, Sequence};
use super::tags::{Array, Tags, Value};
use crate::error::RecordError;

const LONG_CIGAR_TAG: [u8; 2] = *b"CG";

#[derive(Debug, Clone, Copy)]
pub struct RecordRef<'a> {
	core: &'a Core,
	rest: &'a [u8],
	layout: Layout,
}

impl<'a> RecordRef<'a> {
	pub fn new(bytes: &'a [u8]) -> Result<Self, RecordError> {
		let (core, rest) = bytes
			.split_first_chunk::<CORE_LEN>()
			.ok_or(RecordError::TooShort { len: bytes.len() })?;
		let layout = Layout::parse(core, rest.len())?;
		let name = rest.get(..layout.name_end).ok_or(RecordError::Truncated)?;
		if name.last() != Some(&0) {
			return Err(RecordError::InvalidName);
		}
		let cigar = rest
			.get(layout.name_end..layout.cigar_end)
			.ok_or(RecordError::Truncated)?;
		Cigar::new(cigar)?;
		Ok(Self { core, rest, layout })
	}

	pub(crate) const fn from_parts(core: &'a Core, rest: &'a [u8], layout: Layout) -> Self {
		Self { core, rest, layout }
	}

	pub(crate) const fn parts(&self) -> (&'a Core, &'a [u8]) {
		(self.core, self.rest)
	}

	pub fn name(&self) -> &'a [u8] {
		let field = self.rest.get(..self.layout.name_end).unwrap_or_default();
		let end = memchr::memchr(0, field).unwrap_or(field.len());
		field.get(..end).unwrap_or_default()
	}

	pub fn flags(&self) -> Flags {
		Flags::new(fixed::flags(self.core))
	}

	pub fn contig(&self) -> Option<ContigId> {
		contig_id(fixed::contig(self.core))
	}

	pub fn position(&self) -> Option<u32> {
		u32::try_from(fixed::position(self.core)).ok()
	}

	pub fn mapq(&self) -> u8 {
		fixed::mapq(self.core)
	}

	pub fn bin(&self) -> u16 {
		fixed::bin(self.core)
	}

	pub fn mate_contig(&self) -> Option<ContigId> {
		contig_id(fixed::mate_contig(self.core))
	}

	pub fn mate_position(&self) -> Option<u32> {
		u32::try_from(fixed::mate_position(self.core)).ok()
	}

	pub fn template_length(&self) -> i32 {
		fixed::template_len(self.core)
	}

	pub fn stored_cigar(&self) -> Cigar<'a> {
		let bytes = self
			.rest
			.get(self.layout.name_end..self.layout.cigar_end)
			.unwrap_or_default();
		Cigar::trusted(bytes)
	}

	pub fn cigar(&self) -> Cigar<'a> {
		self.long_cigar().unwrap_or_else(|| self.stored_cigar())
	}

	pub fn has_long_cigar(&self) -> bool {
		self.long_cigar().is_some()
	}

	pub fn sequence(&self) -> Sequence<'a> {
		let packed = self
			.rest
			.get(self.layout.cigar_end..self.layout.seq_end)
			.unwrap_or_default();
		Sequence::new(packed, self.layout.seq_len)
	}

	pub fn qualities(&self) -> Qualities<'a> {
		let phred = self
			.rest
			.get(self.layout.seq_end..self.layout.qual_end)
			.unwrap_or_default();
		Qualities::new(phred)
	}

	pub fn tags(&self) -> Tags<'a> {
		Tags::new(self.rest.get(self.layout.qual_end..).unwrap_or_default())
	}

	pub fn reference_end(&self) -> Option<u64> {
		let start = u64::from(self.position()?);
		let length = if self.flags().is_unmapped() {
			0
		} else {
			self.cigar().reference_len()
		};
		Some(start + length.max(1))
	}

	fn long_cigar(&self) -> Option<Cigar<'a>> {
		let stored = self.stored_cigar();
		if stored.len() != 2 {
			return None;
		}
		let mut ops = stored.ops();
		let (first, second) = (ops.next()?, ops.next()?);
		let placeholder = first.kind() == CigarKind::SoftClip
			&& usize::try_from(first.len()).ok() == Some(self.layout.seq_len)
			&& second.kind() == CigarKind::Skip;
		if !placeholder {
			return None;
		}
		match self.tags().get(LONG_CIGAR_TAG) {
			Ok(Some(Value::Array(Array::UInt32(items)))) => Cigar::new(items.as_bytes()).ok(),
			_ => None,
		}
	}
}

fn contig_id(raw: i32) -> Option<ContigId> {
	u32::try_from(raw).ok().map(ContigId::new)
}
