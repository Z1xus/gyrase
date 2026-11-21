use crate::cram::Error;
use crate::cram::cursor::Cursor;

const MAX_ARRAY_LEN: usize = 1 << 28;
const FIRST_RESERVE: usize = 1 << 16;
const MAX_NESTING: u32 = 1;
const MAX_CODE_BITS: u32 = 31;
const TABLE_LEN: usize = 32;
const MAX_SYMBOLS: usize = 1 << 20;

fn malformed() -> Error {
	Error::Malformed("encoding")
}

#[derive(Debug)]
pub(crate) struct BitReader<'a> {
	data: &'a [u8],
	position: usize,
}

impl<'a> BitReader<'a> {
	pub(crate) fn new(data: &'a [u8]) -> Self {
		Self { data, position: 0 }
	}

	fn peek(&self, count: u32) -> u32 {
		let rest = self.data.get(self.position / 8..).unwrap_or_default();
		let mut padded = [0u8; 8];
		if let Some(chunk) = rest.first_chunk::<8>() {
			padded = *chunk;
		} else {
			for (slot, byte) in padded.iter_mut().zip(rest) {
				*slot = *byte;
			}
		}
		let window = u64::from_be_bytes(padded);
		let used = u32::try_from(self.position % 8).unwrap_or_default();
		let shift = 64 - used - count;
		u32::try_from((window >> shift) & ((1u64 << count) - 1)).unwrap_or_default()
	}

	fn skip(&mut self, count: u32) -> Result<(), Error> {
		let end = self.position + count as usize;
		if end > self.data.len() * 8 {
			return Err(Error::Truncated);
		}
		self.position = end;
		Ok(())
	}

	fn bits(&mut self, count: u32) -> Result<u32, Error> {
		if count == 0 {
			return Ok(0);
		}
		let value = self.peek(count);
		self.skip(count)?;
		Ok(value)
	}

	fn bit(&mut self) -> Result<u32, Error> {
		self.bits(1)
	}
}

#[derive(Debug, Default)]
pub(crate) struct External<'a> {
	data: &'a [u8],
	position: usize,
	stop: Option<u8>,
}

impl<'a> External<'a> {
	pub(crate) fn new(data: &'a [u8], stop: Option<u8>) -> Self {
		Self {
			data,
			position: 0,
			stop,
		}
	}

	fn rest(&self) -> &'a [u8] {
		self.data.get(self.position..).unwrap_or_default()
	}

	fn byte(&mut self) -> Result<u8, Error> {
		let Some(&byte) = self.data.get(self.position) else {
			return Err(Error::Truncated);
		};
		self.position += 1;
		Ok(byte)
	}

	fn itf8(&mut self) -> Result<i32, Error> {
		let mut cursor = Cursor::new(self.rest());
		let value = cursor.itf8()?;
		self.position += cursor.position();
		Ok(value)
	}

	fn take(&mut self, count: usize) -> Result<&'a [u8], Error> {
		let Some(end) = self.position.checked_add(count) else {
			return Err(Error::Truncated);
		};
		let Some(bytes) = self.data.get(self.position..end) else {
			return Err(Error::Truncated);
		};
		self.position = end;
		Ok(bytes)
	}

	fn until(&mut self, stop: u8) -> Result<&'a [u8], Error> {
		let rest = self.rest();
		let stop = self.stop.unwrap_or(stop);
		let length = memchr::memchr(stop, rest).ok_or(Error::Truncated)?;
		self.position += length + 1;
		rest.get(..length).ok_or(Error::Truncated)
	}
}

#[derive(Debug)]
pub(crate) struct Data<'a> {
	core: BitReader<'a>,
	externals: Vec<External<'a>>,
}

impl<'a> Data<'a> {
	pub(crate) fn new(core: BitReader<'a>, externals: Vec<External<'a>>) -> Self {
		Self { core, externals }
	}

	fn external(&mut self, slot: usize) -> Result<&mut External<'a>, Error> {
		match self.externals.get_mut(slot) {
			Some(external) => Ok(external),
			None => Err(Error::Truncated),
		}
	}
}

#[derive(Debug, Default)]
pub(crate) struct Slots(Vec<i32>);

impl Slots {
	fn intern(&mut self, id: i32) -> usize {
		if let Some(slot) = self.0.iter().position(|known| *known == id) {
			return slot;
		}
		self.0.push(id);
		self.0.len() - 1
	}

	pub(crate) fn ids(&self) -> &[i32] {
		&self.0
	}
}

#[derive(Debug)]
pub(crate) enum Huffman {
	Constant(i32),
	Canonical(Box<Canonical>),
}

#[derive(Debug)]
pub(crate) struct Canonical {
	symbols: Box<[i32]>,
	limit: [u32; TABLE_LEN],
	first: [u32; TABLE_LEN],
	offset: [u32; TABLE_LEN],
	min_len: u32,
	max_len: u32,
}

impl Huffman {
	fn read(params: &mut Cursor<'_>) -> Result<Self, Error> {
		let count = usize::try_from(params.itf8()?)
			.ok()
			.filter(|count| *count <= MAX_SYMBOLS)
			.ok_or_else(malformed)?;
		let mut codes = Vec::with_capacity(count.min(1024));
		for _ in 0..count {
			codes.push((0u32, params.itf8()?));
		}
		let length_count = usize::try_from(params.itf8()?).map_err(|_| malformed())?;
		if length_count != count {
			return Err(malformed());
		}
		for code in &mut codes {
			code.0 = u32::try_from(params.itf8()?).map_err(|_| malformed())?;
		}
		match codes.as_slice() {
			[] => return Err(malformed()),
			[(0, symbol)] => return Ok(Self::Constant(*symbol)),
			_ => {}
		}
		codes.sort_unstable();
		Canonical::build(&codes).map(|table| Self::Canonical(Box::new(table)))
	}

	fn decode(&self, core: &mut BitReader<'_>) -> Result<i32, Error> {
		match self {
			Self::Constant(symbol) => Ok(*symbol),
			Self::Canonical(table) => table.decode(core),
		}
	}
}

impl Canonical {
	fn build(codes: &[(u32, i32)]) -> Result<Self, Error> {
		let mut counts = [0u32; TABLE_LEN];
		for (length, _) in codes {
			let slot = counts
				.get_mut(*length as usize)
				.filter(|_| (1..=MAX_CODE_BITS).contains(length))
				.ok_or_else(malformed)?;
			*slot += 1;
		}
		let mut table = Self {
			symbols: codes.iter().map(|(_, symbol)| *symbol).collect(),
			limit: [0; TABLE_LEN],
			first: [0; TABLE_LEN],
			offset: [0; TABLE_LEN],
			min_len: 0,
			max_len: 0,
		};
		let mut code = 0u32;
		let mut index = 0u32;
		for length in 1..=MAX_CODE_BITS {
			let count = counts.get(length as usize).copied().unwrap_or(0);
			let end = code.checked_add(count).ok_or_else(malformed)?;
			if end > 1 << length {
				return Err(malformed());
			}
			if let (Some(first), Some(limit), Some(offset)) = (
				table.first.get_mut(length as usize),
				table.limit.get_mut(length as usize),
				table.offset.get_mut(length as usize),
			) {
				*first = code;
				*limit = end;
				*offset = index;
			}
			if count > 0 {
				if table.min_len == 0 {
					table.min_len = length;
				}
				table.max_len = length;
			}
			index += count;
			code = end << 1;
		}
		Ok(table)
	}

	fn decode(&self, core: &mut BitReader<'_>) -> Result<i32, Error> {
		let window = core.peek(self.max_len);
		for length in self.min_len..=self.max_len {
			let prefix = window >> (self.max_len - length);
			let slot = length as usize;
			if prefix < self.limit.get(slot).copied().unwrap_or(0) {
				let first = self.first.get(slot).copied().unwrap_or(0);
				let base = self.offset.get(slot).copied().unwrap_or(0);
				let index = (base + prefix.checked_sub(first).ok_or_else(malformed)?) as usize;
				let symbol = *self.symbols.get(index).ok_or_else(malformed)?;
				core.skip(length)?;
				return Ok(symbol);
			}
		}
		Err(malformed())
	}
}

#[derive(Debug, Default)]
pub(crate) enum Encoding {
	#[default]
	Null,
	External {
		slot: usize,
	},
	Huffman(Huffman),
	ByteArrayLen {
		length: Box<Encoding>,
		values: Box<Encoding>,
	},
	ByteArrayStop {
		stop: u8,
		slot: usize,
	},
	Beta {
		offset: i32,
		bits: u32,
	},
	Subexp {
		offset: i32,
		k: u32,
	},
	Gamma {
		offset: i32,
	},
}

impl Encoding {
	pub(crate) fn read(cursor: &mut Cursor<'_>, slots: &mut Slots) -> Result<Self, Error> {
		Self::read_nested(cursor, slots, 0)
	}

	fn read_nested(cursor: &mut Cursor<'_>, slots: &mut Slots, depth: u32) -> Result<Self, Error> {
		let codec = cursor.itf8()?;
		let size = usize::try_from(cursor.itf8()?).map_err(|_| malformed())?;
		let mut params = Cursor::new(cursor.take(size)?);
		let encoding = match codec {
			0 => Self::Null,
			1 => Self::External {
				slot: slots.intern(params.itf8()?),
			},
			3 => Self::Huffman(Huffman::read(&mut params)?),
			4 if depth < MAX_NESTING => Self::ByteArrayLen {
				length: Box::new(Self::read_nested(&mut params, slots, depth + 1)?),
				values: Box::new(Self::read_nested(&mut params, slots, depth + 1)?),
			},
			5 => Self::ByteArrayStop {
				stop: params.u8()?,
				slot: slots.intern(params.itf8()?),
			},
			6 => Self::Beta {
				offset: params.itf8()?,
				bits: u32::try_from(params.itf8()?)
					.ok()
					.filter(|bits| *bits <= 32)
					.ok_or_else(malformed)?,
			},
			7 => Self::Subexp {
				offset: params.itf8()?,
				k: u32::try_from(params.itf8()?)
					.ok()
					.filter(|k| *k <= 30)
					.ok_or_else(malformed)?,
			},
			9 => Self::Gamma {
				offset: params.itf8()?,
			},
			_ => return Err(malformed()),
		};
		Ok(encoding)
	}

	pub(crate) fn is_null(&self) -> bool {
		matches!(self, Self::Null)
	}

	pub(crate) fn int(&self, data: &mut Data<'_>) -> Result<i32, Error> {
		match self {
			Self::External { slot } => data.external(*slot)?.itf8(),
			Self::Huffman(huffman) => huffman.decode(&mut data.core),
			Self::Beta { offset, bits } => {
				let raw = i64::from(data.core.bits(*bits)?);
				i32::try_from(raw - i64::from(*offset)).map_err(|_| Error::IntegerOverflow)
			}
			Self::Subexp { offset, k } => {
				let mut ones = 0u32;
				while data.core.bit()? == 1 {
					ones += 1;
					if ones + k > 31 {
						return Err(malformed());
					}
				}
				let value = if ones == 0 {
					data.core.bits(*k)?
				} else {
					let tail = ones + k - 1;
					(1u32 << tail) + data.core.bits(tail)?
				};
				i32::try_from(i64::from(value) - i64::from(*offset))
					.map_err(|_| Error::IntegerOverflow)
			}
			Self::Gamma { offset } => {
				let mut zeros = 0u32;
				while data.core.bit()? == 0 {
					zeros += 1;
					if zeros > 31 {
						return Err(malformed());
					}
				}
				let value = (1u32 << zeros) + data.core.bits(zeros)?;
				i32::try_from(i64::from(value) - i64::from(*offset))
					.map_err(|_| Error::IntegerOverflow)
			}
			Self::Null => Err(Error::MissingSeries),
			Self::ByteArrayLen { .. } | Self::ByteArrayStop { .. } => Err(malformed()),
		}
	}

	pub(crate) fn byte(&self, data: &mut Data<'_>) -> Result<u8, Error> {
		match self {
			Self::External { slot } => data.external(*slot)?.byte(),
			_ => u8::try_from(self.int(data)?).map_err(|_| Error::IntegerOverflow),
		}
	}

	pub(crate) fn bytes(
		&self,
		count: usize,
		data: &mut Data<'_>,
		out: &mut Vec<u8>,
	) -> Result<(), Error> {
		if count > MAX_ARRAY_LEN {
			return Err(Error::TooLarge { len: count as u64 });
		}
		if let Self::External { slot } = self {
			out.extend_from_slice(data.external(*slot)?.take(count)?);
		} else {
			out.reserve(count.min(FIRST_RESERVE));
			for _ in 0..count {
				out.push(self.byte(data)?);
			}
		}
		Ok(())
	}

	pub(crate) fn array(&self, data: &mut Data<'_>, out: &mut Vec<u8>) -> Result<(), Error> {
		match self {
			Self::ByteArrayStop { stop, slot } => {
				let bytes = data.external(*slot)?.until(*stop)?;
				out.extend_from_slice(bytes);
				Ok(())
			}
			Self::ByteArrayLen { length, values } => {
				let count = usize::try_from(length.int(data)?).map_err(|_| malformed())?;
				values.bytes(count, data, out)
			}
			Self::Null => Err(Error::MissingSeries),
			_ => Err(malformed()),
		}
	}

	pub(crate) fn tag_value(
		&self,
		kind: u8,
		data: &mut Data<'_>,
		out: &mut Vec<u8>,
	) -> Result<(), Error> {
		match self {
			Self::External { .. } => self.bytes(element_size(kind), data, out),
			_ => self.array(data, out),
		}
	}
}

fn element_size(kind: u8) -> usize {
	match kind {
		b's' | b'S' => 2,
		b'i' | b'I' | b'f' => 4,
		b'd' => 8,
		_ => 1,
	}
}
