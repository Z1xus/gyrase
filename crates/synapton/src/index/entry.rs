use crate::num::{u32_from, usize_of};
use crate::seeding::AUX_BITS;

pub(crate) const TAG_BITS: u32 = 23;
pub(crate) const MAIN_TAG_BITS: u32 = TAG_BITS - AUX_BITS;
pub(crate) const MAX_BUCKET_BITS: u32 = 30;
const ORIENTATION_BITS: u32 = 2;
const OFFSET_BITS: u32 = 7;
const POSITION_BITS: u32 = 32;
const PREFIX_SHIFT: u32 = OFFSET_BITS + POSITION_BITS;
const FIRST_ORIENTATION_BIT: u64 = 0b10;
const OFFSET_MASK: u64 = (1 << OFFSET_BITS) - 1;
const AUX_MASK: u64 = (1 << AUX_BITS) - 1;
const MAIN_TAG_MASK: u64 = (1 << MAIN_TAG_BITS) - 1;
const KEY_MAIN_BITS: u32 = 64 - AUX_BITS;

pub(crate) fn split_key(key: u64, bucket_bits: u32) -> (usize, u64) {
	let bucket = usize_of(key >> (64 - bucket_bits));
	let main_tag =
		(key >> (KEY_MAIN_BITS - bucket_bits - MAIN_TAG_BITS + AUX_BITS)) & MAIN_TAG_MASK;
	(bucket, (main_tag << AUX_BITS) | (key & AUX_MASK))
}

pub(crate) fn prefix(tag: u64, orientation: u8) -> u64 {
	(tag << ORIENTATION_BITS) | u64::from(orientation)
}

pub(crate) fn pack(prefix: u64, offset: u32, position: u32) -> u64 {
	(prefix << PREFIX_SHIFT) | (u64::from(offset) << POSITION_BITS) | u64::from(position)
}

pub(crate) const fn prefix_of(entry: u64) -> u64 {
	entry >> PREFIX_SHIFT
}

pub(crate) fn position_of(entry: u64) -> u32 {
	u32_from(entry & u64::from(u32::MAX))
}

pub(crate) fn offset_of(entry: u64) -> u32 {
	u32_from((entry >> POSITION_BITS) & OFFSET_MASK)
}

pub(crate) const fn main_tag_of(prefix: u64) -> u64 {
	prefix >> (ORIENTATION_BITS + AUX_BITS)
}

pub(crate) const fn first_orientation_of(prefix: u64) -> u64 {
	prefix & FIRST_ORIENTATION_BIT
}
