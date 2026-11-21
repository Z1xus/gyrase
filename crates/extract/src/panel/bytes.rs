use crate::Error;

const MAX_VARINT_BYTES: usize = 10;

#[derive(Debug, Clone, Copy)]
pub(super) struct Bytes<'a> {
	data: &'a [u8],
}

impl<'a> Bytes<'a> {
	pub(super) const fn new(data: &'a [u8]) -> Self {
		Self { data }
	}

	pub(super) const fn len(&self) -> usize {
		self.data.len()
	}

	pub(super) const fn is_empty(&self) -> bool {
		self.data.is_empty()
	}

	pub(super) fn take(&mut self, len: usize) -> Result<&'a [u8], Error> {
		let (head, tail) = self.data.split_at_checked(len).ok_or(Error::Truncated)?;
		self.data = tail;
		Ok(head)
	}

	pub(super) fn array<const N: usize>(&mut self) -> Result<[u8; N], Error> {
		self.take(N)?.try_into().map_err(|_| Error::Truncated)
	}

	pub(super) fn u8(&mut self) -> Result<u8, Error> {
		let [byte] = self.array()?;
		Ok(byte)
	}

	pub(super) fn u16(&mut self) -> Result<u16, Error> {
		Ok(u16::from_le_bytes(self.array()?))
	}

	pub(super) fn u32(&mut self) -> Result<u32, Error> {
		Ok(u32::from_le_bytes(self.array()?))
	}

	pub(super) fn length(&mut self) -> Result<usize, Error> {
		usize::try_from(self.u32()?).map_err(|_| Error::Corrupt("length"))
	}

	pub(super) fn section(&mut self) -> Result<Self, Error> {
		let len = self.length()?;
		Ok(Self::new(self.take(len)?))
	}

	pub(super) fn varint(&mut self) -> Result<u64, Error> {
		let mut value = 0u64;
		for index in 0..MAX_VARINT_BYTES {
			let byte = self.u8()?;
			value |= u64::from(byte & 0x7f) << (7 * index);
			if byte & 0x80 == 0 {
				return Ok(value);
			}
		}
		Err(Error::Corrupt("varint"))
	}
}

pub(super) fn put_varint(buffer: &mut Vec<u8>, mut value: u64) {
	while value >= 0x80 {
		buffer.push(u8::try_from(value & 0x7f).unwrap_or(0) | 0x80);
		value >>= 7;
	}
	buffer.push(u8::try_from(value).unwrap_or(0));
}

pub(super) const fn zigzag(value: i64) -> u64 {
	((value << 1) ^ (value >> 63)).cast_unsigned()
}

pub(super) const fn unzigzag(value: u64) -> i64 {
	(value >> 1).cast_signed() ^ -((value & 1).cast_signed())
}
