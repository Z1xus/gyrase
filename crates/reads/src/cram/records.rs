#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Span {
	pub(crate) start: usize,
	pub(crate) end: usize,
}

impl Span {
	pub(crate) fn since(start: usize, end: usize) -> Self {
		Self { start, end }
	}

	pub(crate) fn of<T>(self, arena: &[T]) -> &[T] {
		arena.get(self.start..self.end).unwrap_or_default()
	}

	pub(crate) fn is_empty(self) -> bool {
		self.end <= self.start
	}
}

#[derive(Debug, Default)]
pub(crate) struct Record {
	pub(crate) flags: u16,
	pub(crate) cram_flags: u32,
	pub(crate) reference: i32,
	pub(crate) position: i32,
	pub(crate) end: i32,
	pub(crate) length: usize,
	pub(crate) mapq: u8,
	pub(crate) read_group: i32,
	pub(crate) mate_line: Option<usize>,
	pub(crate) mate_flags: u8,
	pub(crate) mate_reference: i32,
	pub(crate) mate_position: i32,
	pub(crate) template_length: Option<i32>,
	pub(crate) explicit_template_length: Option<i32>,
	pub(crate) name: Span,
	pub(crate) cigar: Span,
	pub(crate) bases: Span,
	pub(crate) qualities: Span,
	pub(crate) tags: Span,
}

#[derive(Debug, Default)]
pub(crate) struct RecordTable {
	pub(crate) records: Vec<Record>,
	pub(crate) names: Vec<u8>,
	pub(crate) cigars: Vec<u32>,
	pub(crate) bases: Vec<u8>,
	pub(crate) qualities: Vec<u8>,
	pub(crate) tags: Vec<u8>,
}

pub(crate) const FLAG_PAIRED: u16 = 0x1;
pub(crate) const FLAG_UNMAPPED: u16 = 0x4;
pub(crate) const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub(crate) const FLAG_REVERSE: u16 = 0x10;
pub(crate) const FLAG_MATE_REVERSE: u16 = 0x20;
pub(crate) const FLAG_FIRST: u16 = 0x40;

pub(crate) const CF_QUALITIES: u32 = 0x1;
pub(crate) const CF_DETACHED: u32 = 0x2;
pub(crate) const CF_MATE_DOWNSTREAM: u32 = 0x4;
pub(crate) const CF_NO_SEQUENCE: u32 = 0x8;
pub(crate) const CF_EXPLICIT_LENGTH: u32 = 0x10;

pub(crate) const MATE_REVERSE: u8 = 0x1;
pub(crate) const MATE_UNMAPPED: u8 = 0x2;

pub(crate) const CIGAR_MATCH: u32 = 0;
pub(crate) const CIGAR_INSERTION: u32 = 1;
pub(crate) const CIGAR_DELETION: u32 = 2;
pub(crate) const CIGAR_SKIP: u32 = 3;
pub(crate) const CIGAR_SOFT_CLIP: u32 = 4;
pub(crate) const CIGAR_HARD_CLIP: u32 = 5;
pub(crate) const CIGAR_PADDING: u32 = 6;
