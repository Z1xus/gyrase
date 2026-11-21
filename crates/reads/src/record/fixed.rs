use crate::bytes::{i32_le, u16_le, u32_le};
use crate::error::RecordError;

pub(crate) const CORE_LEN: usize = 32;
pub(crate) type Core = [u8; CORE_LEN];

const CONTIG: usize = 0;
const POSITION: usize = 4;
const NAME_LEN: usize = 8;
const MAPQ: usize = 9;
const BIN: usize = 10;
const CIGAR_OPS: usize = 12;
const FLAGS: usize = 14;
const SEQ_LEN: usize = 16;
const MATE_CONTIG: usize = 20;
const MATE_POSITION: usize = 24;
const TEMPLATE_LEN: usize = 28;

pub(crate) fn contig(core: &Core) -> i32 {
	i32_le(core, CONTIG).unwrap_or(-1)
}

pub(crate) fn position(core: &Core) -> i32 {
	i32_le(core, POSITION).unwrap_or(-1)
}

pub(crate) fn name_len(core: &Core) -> usize {
	usize::from(core.get(NAME_LEN).copied().unwrap_or(0))
}

pub(crate) fn mapq(core: &Core) -> u8 {
	core.get(MAPQ).copied().unwrap_or(0)
}

pub(crate) fn bin(core: &Core) -> u16 {
	u16_le(core, BIN).unwrap_or(0)
}

pub(crate) fn cigar_ops(core: &Core) -> usize {
	usize::from(u16_le(core, CIGAR_OPS).unwrap_or(0))
}

pub(crate) fn flags(core: &Core) -> u16 {
	u16_le(core, FLAGS).unwrap_or(0)
}

pub(crate) fn seq_len(core: &Core) -> usize {
	u32_le(core, SEQ_LEN).map_or(0, |len| len as usize)
}

pub(crate) fn mate_contig(core: &Core) -> i32 {
	i32_le(core, MATE_CONTIG).unwrap_or(-1)
}

pub(crate) fn mate_position(core: &Core) -> i32 {
	i32_le(core, MATE_POSITION).unwrap_or(-1)
}

pub(crate) fn template_len(core: &Core) -> i32 {
	i32_le(core, TEMPLATE_LEN).unwrap_or(0)
}

pub(crate) fn put(core: &mut Core, at: usize, value: &[u8]) {
	if let Some(slot) = core.get_mut(at..at + value.len()) {
		slot.copy_from_slice(value);
	}
}

pub(crate) fn set_contig(core: &mut Core, value: i32) {
	put(core, CONTIG, &value.to_le_bytes());
}

pub(crate) fn set_position(core: &mut Core, value: i32) {
	put(core, POSITION, &value.to_le_bytes());
}

pub(crate) fn set_name_len(core: &mut Core, value: u8) {
	put(core, NAME_LEN, &[value]);
}

pub(crate) fn set_mapq(core: &mut Core, value: u8) {
	put(core, MAPQ, &[value]);
}

pub(crate) fn set_bin(core: &mut Core, value: u16) {
	put(core, BIN, &value.to_le_bytes());
}

pub(crate) fn set_cigar_ops(core: &mut Core, value: u16) {
	put(core, CIGAR_OPS, &value.to_le_bytes());
}

pub(crate) fn set_flags(core: &mut Core, value: u16) {
	put(core, FLAGS, &value.to_le_bytes());
}

pub(crate) fn set_seq_len(core: &mut Core, value: u32) {
	put(core, SEQ_LEN, &value.to_le_bytes());
}

pub(crate) fn set_mate_contig(core: &mut Core, value: i32) {
	put(core, MATE_CONTIG, &value.to_le_bytes());
}

pub(crate) fn set_mate_position(core: &mut Core, value: i32) {
	put(core, MATE_POSITION, &value.to_le_bytes());
}

pub(crate) fn set_template_len(core: &mut Core, value: i32) {
	put(core, TEMPLATE_LEN, &value.to_le_bytes());
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Layout {
	pub name_end: usize,
	pub cigar_end: usize,
	pub seq_end: usize,
	pub qual_end: usize,
	pub seq_len: usize,
}

impl Layout {
	pub(crate) fn from_counts(name_len: usize, cigar_ops: usize, seq_len: usize) -> Option<Self> {
		let cigar_end = name_len.checked_add(cigar_ops.checked_mul(4)?)?;
		let seq_end = cigar_end.checked_add(seq_len.checked_add(1)? / 2)?;
		let qual_end = seq_end.checked_add(seq_len)?;
		Some(Self {
			name_end: name_len,
			cigar_end,
			seq_end,
			qual_end,
			seq_len,
		})
	}

	pub(crate) fn parse(core: &Core, rest_len: usize) -> Result<Self, RecordError> {
		if name_len(core) == 0 {
			return Err(RecordError::InvalidName);
		}
		let layout = Self::from_counts(name_len(core), cigar_ops(core), seq_len(core))
			.filter(|layout| layout.qual_end <= rest_len)
			.ok_or(RecordError::Truncated)?;
		if i32::try_from(layout.seq_len).is_err()
			|| contig(core) < -1
			|| position(core) < -1
			|| mate_contig(core) < -1
			|| mate_position(core) < -1
		{
			return Err(RecordError::InvalidCoordinate);
		}
		Ok(layout)
	}
}
