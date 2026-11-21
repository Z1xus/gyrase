use crate::cram::Error;

#[derive(Debug, Clone)]
pub(crate) struct Cursor<'a> {
	data: &'a [u8],
	position: usize,
}

impl<'a> Cursor<'a> {
	pub(crate) fn new(data: &'a [u8]) -> Self {
		Self { data, position: 0 }
	}

	pub(crate) fn rest(&self) -> &'a [u8] {
		self.data.get(self.position..).unwrap_or_default()
	}

	pub(crate) fn remaining(&self) -> usize {
		self.rest().len()
	}

	pub(crate) fn preceding(&self, start: usize, len: usize) -> Option<&'a [u8]> {
		self.data.get(start..start.checked_add(len)?)
	}

	pub(crate) fn position(&self) -> usize {
		self.position
	}

	pub(crate) fn is_empty(&self) -> bool {
		self.remaining() == 0
	}

	pub(crate) fn peek(&self) -> Result<u8, Error> {
		match self.data.get(self.position) {
			Some(&byte) => Ok(byte),
			None => Err(Error::Truncated),
		}
	}

	pub(crate) fn u8(&mut self) -> Result<u8, Error> {
		let byte = self.peek()?;
		self.position += 1;
		Ok(byte)
	}

	pub(crate) fn take(&mut self, len: usize) -> Result<&'a [u8], Error> {
		let Some(end) = self.position.checked_add(len) else {
			return Err(Error::Truncated);
		};
		let Some(bytes) = self.data.get(self.position..end) else {
			return Err(Error::Truncated);
		};
		self.position = end;
		Ok(bytes)
	}

	pub(crate) fn array<const N: usize>(&mut self) -> Result<[u8; N], Error> {
		self.take(N)?.try_into().map_err(|_| Error::Truncated)
	}

	pub(crate) fn u16_le(&mut self) -> Result<u16, Error> {
		Ok(u16::from_le_bytes(self.array()?))
	}

	pub(crate) fn u32_le(&mut self) -> Result<u32, Error> {
		Ok(u32::from_le_bytes(self.array()?))
	}

	pub(crate) fn i32_le(&mut self) -> Result<i32, Error> {
		Ok(i32::from_le_bytes(self.array()?))
	}

	pub(crate) fn uint7(&mut self) -> Result<u32, Error> {
		let mut value = 0u32;
		for _ in 0..5 {
			let byte = self.u8()?;
			value = value.checked_mul(128).ok_or(Error::IntegerOverflow)? | u32::from(byte & 0x7f);
			if byte < 0x80 {
				return Ok(value);
			}
		}
		Err(Error::IntegerOverflow)
	}

	pub(crate) fn itf8(&mut self) -> Result<i32, Error> {
		let first = self.u8()?;
		let extra = first.leading_ones().min(4) as usize;
		let mask = if extra >= 4 {
			0x0f
		} else {
			0xff >> (extra + 1)
		};
		let mut value = u32::from(first) & mask;
		for index in 0..extra {
			let byte = u32::from(self.u8()?);
			if extra == 4 && index == 3 {
				value = (value << 4) | (byte & 0x0f);
			} else {
				value = (value << 8) | byte;
			}
		}
		Ok(value.cast_signed())
	}

	pub(crate) fn ltf8(&mut self) -> Result<i64, Error> {
		let first = self.u8()?;
		let extra = first.leading_ones() as usize;
		let mut value = if extra >= 7 {
			0
		} else {
			u64::from(first) & (0xff >> (extra + 1))
		};
		for _ in 0..extra {
			value = (value << 8) | u64::from(self.u8()?);
		}
		Ok(value.cast_signed())
	}
}
