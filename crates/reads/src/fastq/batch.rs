const MAX_BATCH_BYTES: usize = 1 << 30;

#[derive(Debug, Clone, Copy)]
pub(super) struct Span {
	pub name_start: u32,
	pub name_end: u32,
	pub header_end: u32,
	pub bases_start: u32,
	pub bases_end: u32,
	pub qualities_start: u32,
	pub qualities_end: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadRef<'a> {
	pub name: &'a [u8],
	pub comment: &'a [u8],
	pub bases: &'a [u8],
	pub qualities: &'a [u8],
}

#[derive(Debug, Default)]
pub struct FastqBatch {
	pub(super) data: Vec<u8>,
	pub(super) spans: Vec<Span>,
	pub(super) consumed: usize,
}

impl FastqBatch {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn len(&self) -> usize {
		self.spans.len()
	}

	pub fn is_empty(&self) -> bool {
		self.spans.is_empty()
	}

	pub fn get(&self, index: usize) -> Option<ReadRef<'_>> {
		let span = self.spans.get(index)?;
		let range = |start: u32, end: u32| self.data.get(start as usize..end as usize);
		let comment_start = self
			.data
			.get(span.name_end as usize..span.header_end as usize)
			.map_or(span.header_end, |tail| {
				let blanks = tail
					.iter()
					.take_while(|byte| matches!(byte, b' ' | b'\t'))
					.count();
				span.name_end + u32::try_from(blanks).unwrap_or(0)
			});
		Some(ReadRef {
			name: range(span.name_start, span.name_end)?,
			comment: range(comment_start, span.header_end)?,
			bases: range(span.bases_start, span.bases_end)?,
			qualities: range(span.qualities_start, span.qualities_end)?,
		})
	}

	pub fn iter(&self) -> impl Iterator<Item = ReadRef<'_>> {
		(0..self.len()).filter_map(|index| self.get(index))
	}

	pub(super) fn is_full(&self) -> bool {
		!self.spans.is_empty() && self.data.len() >= MAX_BATCH_BYTES
	}

	pub(super) fn restart(&mut self) {
		self.data.copy_within(self.consumed.., 0);
		self.data.truncate(self.data.len() - self.consumed);
		self.consumed = 0;
		self.spans.clear();
	}
}
