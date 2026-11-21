use crate::Error;
use crate::value::{BcfType, Cells, END_INT, MISSING_INT, Store};

const LONG_LENGTH: u8 = 15;
const INT8_MISSING: i8 = i8::MIN;
const INT8_END: i8 = i8::MIN + 1;
const INT16_MISSING: i16 = i16::MIN;
const INT16_END: i16 = i16::MIN + 1;

#[derive(Debug, Clone, Copy)]
pub struct Typed<'a> {
	pub bcf_type: BcfType,
	pub len: usize,
	pub bytes: &'a [u8],
}

pub(crate) struct Cursor<'a> {
	data: &'a [u8],
	at: usize,
}

impl<'a> Cursor<'a> {
	pub fn new(data: &'a [u8]) -> Self {
		Self { data, at: 0 }
	}

	pub fn at(data: &'a [u8], at: usize) -> Self {
		Self { data, at }
	}

	pub fn position(&self) -> usize {
		self.at
	}

	pub fn take(&mut self, len: usize) -> Result<&'a [u8], Error> {
		let end = self
			.at
			.checked_add(len)
			.ok_or(Error::Bcf("length overflow"))?;
		let bytes = self
			.data
			.get(self.at..end)
			.ok_or(Error::Bcf("record is truncated"))?;
		self.at = end;
		Ok(bytes)
	}

	fn descriptor(&mut self) -> Result<(BcfType, usize), Error> {
		let byte = *self
			.take(1)?
			.first()
			.ok_or(Error::Bcf("record is truncated"))?;
		let bcf_type = BcfType::from_code(byte & 0x0f).ok_or(Error::Bcf("unknown value type"))?;
		let len = byte >> 4;
		if len != LONG_LENGTH {
			return Ok((bcf_type, usize::from(len)));
		}
		let length = self.typed()?;
		let value = length
			.scalar_int()
			.and_then(|value| usize::try_from(value).ok())
			.ok_or(Error::Bcf("invalid vector length"))?;
		Ok((bcf_type, value))
	}

	pub fn typed(&mut self) -> Result<Typed<'a>, Error> {
		let (bcf_type, len) = self.descriptor()?;
		let width = len
			.checked_mul(bcf_type.width())
			.ok_or(Error::Bcf("length overflow"))?;
		let bytes = self.take(width)?;
		Ok(Typed {
			bcf_type,
			len,
			bytes,
		})
	}

	pub fn typed_repeated(&mut self, repeat: usize) -> Result<Typed<'a>, Error> {
		let (bcf_type, len) = self.descriptor()?;
		let width = len
			.checked_mul(repeat)
			.and_then(|total| total.checked_mul(bcf_type.width()))
			.ok_or(Error::Bcf("length overflow"))?;
		let bytes = self.take(width)?;
		Ok(Typed {
			bcf_type,
			len,
			bytes,
		})
	}
}

impl Typed<'_> {
	pub fn scalar_int(&self) -> Option<i32> {
		if self.len != 1 {
			return None;
		}
		self.ints().next()
	}

	pub fn ints(&self) -> impl Iterator<Item = i32> + '_ {
		let bcf_type = self.bcf_type;
		let step = bcf_type.width();
		self.bytes
			.chunks_exact(step)
			.filter_map(move |chunk| int_from(bcf_type, chunk))
	}

	pub fn floats(&self) -> impl Iterator<Item = f32> + '_ {
		float_values(self.bytes)
	}

	pub fn text(&self) -> &[u8] {
		self.bytes
	}

	pub fn push_to(&self, store: &mut Store) -> Cells {
		match self.bcf_type {
			BcfType::Null => Cells::Flag,
			BcfType::Int8 | BcfType::Int16 | BcfType::Int32 => store.push_ints(self.ints()),
			BcfType::Float => store.push_floats(float_values(self.bytes)),
			BcfType::Char => store.push_text(self.bytes),
		}
	}
}

pub(crate) fn float_values(bytes: &[u8]) -> impl Iterator<Item = f32> + '_ {
	bytes
		.as_chunks::<4>()
		.0
		.iter()
		.map(|chunk| f32::from_bits(u32::from_le_bytes(*chunk)))
}

fn int_from(bcf_type: BcfType, chunk: &[u8]) -> Option<i32> {
	match bcf_type {
		BcfType::Int8 => {
			let value = i8::from_le_bytes([*chunk.first()?]);
			Some(match value {
				INT8_MISSING => MISSING_INT,
				INT8_END => END_INT,
				other => i32::from(other),
			})
		}
		BcfType::Int16 => {
			let value = i16::from_le_bytes(chunk.try_into().ok()?);
			Some(match value {
				INT16_MISSING => MISSING_INT,
				INT16_END => END_INT,
				other => i32::from(other),
			})
		}
		BcfType::Int32 => Some(i32::from_le_bytes(chunk.try_into().ok()?)),
		_ => None,
	}
}
