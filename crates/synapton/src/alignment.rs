use std::ops::Range;

use gyrase_genome::ContigId;

use crate::num::u32_of;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
	Forward,
	Reverse,
}

impl Strand {
	pub(crate) const fn index(self) -> usize {
		match self {
			Self::Forward => 0,
			Self::Reverse => 1,
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarKind {
	Match,
	Insertion,
	Deletion,
	SoftClip,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp(u32);

impl CigarOp {
	const KIND_BITS: u32 = 4;

	pub const fn new(kind: CigarKind, length: u32) -> Self {
		let code = match kind {
			CigarKind::Match => 0,
			CigarKind::Insertion => 1,
			CigarKind::Deletion => 2,
			CigarKind::SoftClip => 4,
		};
		Self((length << Self::KIND_BITS) | code)
	}

	pub const fn kind(self) -> CigarKind {
		match self.0 & ((1 << Self::KIND_BITS) - 1) {
			1 => CigarKind::Insertion,
			2 => CigarKind::Deletion,
			4 => CigarKind::SoftClip,
			_ => CigarKind::Match,
		}
	}

	pub const fn length(self) -> u32 {
		self.0 >> Self::KIND_BITS
	}

	pub const fn packed(self) -> u32 {
		self.0
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Flag {
	Paired,
	ProperPair,
	Unmapped,
	MateUnmapped,
	Reverse,
	MateReverse,
	First,
	Second,
	Secondary,
	Supplementary,
}

impl Flag {
	const fn bit(self) -> u16 {
		match self {
			Self::Paired => 0x1,
			Self::ProperPair => 0x2,
			Self::Unmapped => 0x4,
			Self::MateUnmapped => 0x8,
			Self::Reverse => 0x10,
			Self::MateReverse => 0x20,
			Self::First => 0x40,
			Self::Second => 0x80,
			Self::Secondary => 0x100,
			Self::Supplementary => 0x800,
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct FlagSet(u16);

impl FlagSet {
	pub const fn with(self, flag: Flag) -> Self {
		Self(self.0 | flag.bit())
	}

	pub(crate) const fn without(self, flag: Flag) -> Self {
		Self(self.0 & !flag.bit())
	}

	pub const fn contains(self, flag: Flag) -> bool {
		self.0 & flag.bit() != 0
	}

	pub const fn bits(self) -> u16 {
		self.0
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum End {
	Single,
	First,
	Second,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadId {
	pub index: u32,
	pub end: End,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MateLocation {
	pub contig: ContigId,
	pub position: u32,
	pub strand: Strand,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mate {
	Unpaired,
	Unmapped,
	Mapped(MateLocation),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
	pub contig: ContigId,
	pub position: u32,
	pub strand: Strand,
	pub flags: FlagSet,
	pub mapq: u8,
	pub score: i32,
	pub suboptimal_score: Option<i32>,
	pub edit_distance: u32,
	pub mate: Mate,
	pub template_length: i32,
	pub(crate) cigar: Range<u32>,
	pub(crate) md: Range<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Placement {
	Mapped(Alignment),
	Unmapped { flags: FlagSet, mate: Mate },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
	pub read: ReadId,
	pub placement: Placement,
}

#[derive(Debug, Clone, Default)]
pub struct AlignedBatch {
	records: Vec<Record>,
	cigar: Vec<CigarOp>,
	md: Vec<u8>,
}

impl AlignedBatch {
	pub fn records(&self) -> &[Record] {
		&self.records
	}

	pub fn cigar(&self, alignment: &Alignment) -> &[CigarOp] {
		self.cigar
			.get(alignment.cigar.start as usize..alignment.cigar.end as usize)
			.unwrap_or_default()
	}

	pub fn md(&self, alignment: &Alignment) -> &[u8] {
		self.md
			.get(alignment.md.start as usize..alignment.md.end as usize)
			.unwrap_or_default()
	}

	pub(crate) fn push_cigar(&mut self, operations: &[CigarOp]) -> Range<u32> {
		let start = u32_of(self.cigar.len());
		self.cigar.extend_from_slice(operations);
		start..u32_of(self.cigar.len())
	}

	pub(crate) fn push_md(&mut self, md: &[u8]) -> Range<u32> {
		let start = u32_of(self.md.len());
		self.md.extend_from_slice(md);
		start..u32_of(self.md.len())
	}

	pub(crate) fn push_record(&mut self, record: Record) {
		self.records.push(record);
	}

	pub(crate) fn absorb(&mut self, other: Self) {
		let cigar_base = u32_of(self.cigar.len());
		let md_base = u32_of(self.md.len());
		self.cigar.extend(other.cigar);
		self.md.extend(other.md);
		self.records
			.extend(other.records.into_iter().map(|mut record| {
				if let Placement::Mapped(alignment) = &mut record.placement {
					alignment.cigar =
						alignment.cigar.start + cigar_base..alignment.cigar.end + cigar_base;
					alignment.md = alignment.md.start + md_base..alignment.md.end + md_base;
				}
				record
			}));
	}
}
