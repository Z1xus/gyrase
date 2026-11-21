use crate::bytes::u32_le;
use crate::error::RecordError;

const KIND_MASK: u32 = 0xf;
const LENGTH_SHIFT: u32 = 4;
const SYMBOLS: &[u8; 9] = b"MIDNSHP=X";
pub(crate) const MAX_OP_LEN: u32 = (1 << 28) - 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CigarKind {
	Match,
	Insertion,
	Deletion,
	Skip,
	SoftClip,
	HardClip,
	Padding,
	SeqMatch,
	SeqMismatch,
}

impl CigarKind {
	pub const fn from_code(code: u8) -> Option<Self> {
		Some(match code {
			0 => Self::Match,
			1 => Self::Insertion,
			2 => Self::Deletion,
			3 => Self::Skip,
			4 => Self::SoftClip,
			5 => Self::HardClip,
			6 => Self::Padding,
			7 => Self::SeqMatch,
			8 => Self::SeqMismatch,
			_ => return None,
		})
	}

	pub fn from_symbol(symbol: u8) -> Option<Self> {
		let code = SYMBOLS.iter().position(|&candidate| candidate == symbol)?;
		Self::from_code(u8::try_from(code).ok()?)
	}

	pub const fn code(self) -> u8 {
		match self {
			Self::Match => 0,
			Self::Insertion => 1,
			Self::Deletion => 2,
			Self::Skip => 3,
			Self::SoftClip => 4,
			Self::HardClip => 5,
			Self::Padding => 6,
			Self::SeqMatch => 7,
			Self::SeqMismatch => 8,
		}
	}

	pub const fn symbol(self) -> u8 {
		match self {
			Self::Match => b'M',
			Self::Insertion => b'I',
			Self::Deletion => b'D',
			Self::Skip => b'N',
			Self::SoftClip => b'S',
			Self::HardClip => b'H',
			Self::Padding => b'P',
			Self::SeqMatch => b'=',
			Self::SeqMismatch => b'X',
		}
	}

	pub const fn consumes_reference(self) -> bool {
		matches!(
			self,
			Self::Match | Self::Deletion | Self::Skip | Self::SeqMatch | Self::SeqMismatch
		)
	}

	pub const fn consumes_query(self) -> bool {
		matches!(
			self,
			Self::Match | Self::Insertion | Self::SoftClip | Self::SeqMatch | Self::SeqMismatch
		)
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CigarOp {
	len: u32,
	kind: CigarKind,
}

impl CigarOp {
	pub const fn new(len: u32, kind: CigarKind) -> Option<Self> {
		if len > MAX_OP_LEN {
			return None;
		}
		Some(Self { len, kind })
	}

	pub const fn len(self) -> u32 {
		self.len
	}

	pub const fn is_empty(self) -> bool {
		self.len == 0
	}

	pub const fn kind(self) -> CigarKind {
		self.kind
	}

	pub const fn raw(self) -> u32 {
		self.len << LENGTH_SHIFT | self.kind.code() as u32
	}

	fn from_raw(raw: u32) -> Option<Self> {
		let code = u8::try_from(raw & KIND_MASK).ok()?;
		Some(Self {
			len: raw >> LENGTH_SHIFT,
			kind: CigarKind::from_code(code)?,
		})
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cigar<'a> {
	bytes: &'a [u8],
}

impl<'a> Cigar<'a> {
	pub(crate) fn new(bytes: &'a [u8]) -> Result<Self, RecordError> {
		if !bytes.len().is_multiple_of(4) {
			return Err(RecordError::Truncated);
		}
		for raw in raw_ops(bytes) {
			let code = u8::try_from(raw & KIND_MASK).map_err(|_| RecordError::Truncated)?;
			CigarKind::from_code(code).ok_or(RecordError::InvalidCigarOp { code })?;
		}
		Ok(Self { bytes })
	}

	pub(crate) const fn trusted(bytes: &'a [u8]) -> Self {
		Self { bytes }
	}

	pub fn len(&self) -> usize {
		self.bytes.len() / 4
	}

	pub fn is_empty(&self) -> bool {
		self.bytes.is_empty()
	}

	pub fn ops(&self) -> impl Iterator<Item = CigarOp> + 'a {
		raw_ops(self.bytes).filter_map(CigarOp::from_raw)
	}

	pub fn reference_len(&self) -> u64 {
		self.ops()
			.filter(|op| op.kind.consumes_reference())
			.map(|op| u64::from(op.len))
			.sum()
	}

	pub fn query_len(&self) -> u64 {
		self.ops()
			.filter(|op| op.kind.consumes_query())
			.map(|op| u64::from(op.len))
			.sum()
	}
}

fn raw_ops(bytes: &[u8]) -> impl Iterator<Item = u32> + '_ {
	(0..bytes.len() / 4).filter_map(|index| u32_le(bytes, index * 4))
}
