use gyrase_io::{Block, VirtualOffset};

use crate::error::RecordError;
use crate::record::RecordRef;

#[derive(Debug, Clone, Copy)]
pub(super) enum Location {
	Block { index: u32, offset: u32, len: u32 },
	Stitched { offset: u32, len: u32 },
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Entry {
	pub location: Location,
	pub start: VirtualOffset,
	pub end: VirtualOffset,
}

#[derive(Debug, Clone, Copy)]
pub struct Located<'a> {
	pub record: RecordRef<'a>,
	pub start: VirtualOffset,
	pub end: VirtualOffset,
}

#[derive(Debug, Default)]
pub struct RecordBatch {
	blocks: Vec<Block>,
	stitched: Vec<u8>,
	entries: Vec<Entry>,
	reference_count: u32,
}

impl RecordBatch {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn len(&self) -> usize {
		self.entries.len()
	}

	pub fn is_empty(&self) -> bool {
		self.entries.is_empty()
	}

	pub fn at(&self, index: usize) -> Result<Located<'_>, RecordError> {
		let entry = self.entries.get(index).ok_or(RecordError::Truncated)?;
		let bytes =
			slice(&self.blocks, &self.stitched, entry.location).ok_or(RecordError::Truncated)?;
		let record = RecordRef::new(bytes)?;
		for contig in [record.contig(), record.mate_contig()]
			.into_iter()
			.flatten()
		{
			if contig.get() >= self.reference_count {
				return Err(RecordError::ContigOutOfRange {
					contig: contig.get(),
				});
			}
		}
		Ok(Located {
			record,
			start: entry.start,
			end: entry.end,
		})
	}

	pub fn iter(&self) -> impl Iterator<Item = Result<Located<'_>, RecordError>> {
		(0..self.len()).map(|index| self.at(index))
	}

	pub(super) fn reset(&mut self, blocks: Vec<Block>, reference_count: u32) {
		self.blocks = blocks;
		self.stitched.clear();
		self.entries.clear();
		self.reference_count = reference_count;
	}

	pub(super) fn parts_mut(&mut self) -> (&[Block], &mut Vec<u8>, &mut Vec<Entry>) {
		(&self.blocks, &mut self.stitched, &mut self.entries)
	}

	pub(super) fn retain(&mut self, mut keep: impl FnMut(RecordRef<'_>) -> bool) {
		let (blocks, stitched) = (&self.blocks, &self.stitched);
		self.entries.retain(|entry| {
			slice(blocks, stitched, entry.location)
				.and_then(|bytes| RecordRef::new(bytes).ok())
				.is_none_or(&mut keep)
		});
	}
}

fn slice<'a>(blocks: &'a [Block], stitched: &'a [u8], location: Location) -> Option<&'a [u8]> {
	match location {
		Location::Block { index, offset, len } => {
			let start = offset as usize;
			blocks
				.get(index as usize)?
				.data
				.get(start..start.checked_add(len as usize)?)
		}
		Location::Stitched { offset, len } => {
			let start = offset as usize;
			stitched.get(start..start.checked_add(len as usize)?)
		}
	}
}
