use std::num::NonZeroU32;

use crate::Strand;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Change {
	Insertion { start: u32, len: u32 },
	Deletion { len: u32 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Indel {
	pub anchor: u64,
	pub change: Change,
	pub mapping_quality: u8,
	pub strand: Strand,
	pub flank: u8,
	pub fragment: Option<NonZeroU32>,
}

pub(crate) fn left_align_deletion(reference: &[u8], mut anchor: usize, len: usize) -> usize {
	while anchor > 0
		&& reference.get(anchor) == reference.get(anchor + len)
		&& reference.get(anchor).is_some()
	{
		anchor -= 1;
	}
	anchor
}

pub(crate) fn left_align_insertion(
	reference: &[u8],
	mut anchor: usize,
	inserted: &mut [u8],
) -> usize {
	while anchor > 0
		&& inserted
			.last()
			.is_some_and(|last| reference.get(anchor) == Some(last))
	{
		inserted.rotate_right(1);
		anchor -= 1;
	}
	anchor
}
