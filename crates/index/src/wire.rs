use std::io::Read;

use crate::Error;

pub(crate) const MAX_INDEX_BYTES: u64 = 1 << 30;
pub(crate) const MAX_REFERENCES: usize = 1 << 20;

pub(crate) fn read_all(reader: impl Read) -> Result<Vec<u8>, Error> {
	let mut bytes = Vec::new();
	reader.take(MAX_INDEX_BYTES + 1).read_to_end(&mut bytes)?;
	if bytes.len() as u64 > MAX_INDEX_BYTES {
		return Err(Error::TooLarge {
			limit: MAX_INDEX_BYTES,
		});
	}
	Ok(bytes)
}

#[derive(Debug)]
pub(crate) struct Cursor<'a> {
	bytes: &'a [u8],
}

impl<'a> Cursor<'a> {
	pub(crate) const fn new(bytes: &'a [u8]) -> Self {
		Self { bytes }
	}

	pub(crate) fn take(&mut self, len: usize) -> Result<&'a [u8], Error> {
		let (head, tail) = self.bytes.split_at_checked(len).ok_or(Error::Truncated)?;
		self.bytes = tail;
		Ok(head)
	}

	pub(crate) fn array<const N: usize>(&mut self) -> Result<[u8; N], Error> {
		self.take(N)?.try_into().map_err(|_| Error::Truncated)
	}

	pub(crate) fn u32(&mut self) -> Result<u32, Error> {
		Ok(u32::from_le_bytes(self.array()?))
	}

	pub(crate) fn u64(&mut self) -> Result<u64, Error> {
		Ok(u64::from_le_bytes(self.array()?))
	}

	pub(crate) fn optional_u64(&mut self) -> u64 {
		self.u64().unwrap_or(0)
	}

	pub(crate) fn count(&mut self, unit: usize) -> Result<usize, Error> {
		let raw = self.u32()?;
		self.validate_count(raw, unit)
	}

	pub(crate) fn validate_count(&self, raw: u32, unit: usize) -> Result<usize, Error> {
		let invalid = || Error::InvalidCount {
			count: u64::from(raw),
		};
		let count = i32::try_from(raw)
			.ok()
			.and_then(|count| usize::try_from(count).ok())
			.ok_or_else(invalid)?;
		let bytes = count.checked_mul(unit).ok_or_else(invalid)?;
		if bytes > self.bytes.len() {
			return Err(invalid());
		}
		Ok(count)
	}
}
