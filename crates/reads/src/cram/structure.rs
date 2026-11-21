use std::borrow::Cow;
use std::io::Read;

use libdeflater::Crc;

use crate::cram::Error;
use crate::cram::codec::{Method, check_len};
use crate::cram::cursor::Cursor;

const MAGIC: [u8; 4] = *b"CRAM";
const FILE_ID_LEN: usize = 20;
const EOF_START: i32 = 4_542_278;
const HEADER_PROBE: usize = 64;
const MAX_LANDMARKS: usize = 1 << 20;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Version {
	V3_0,
	V3_1,
}

impl Version {
	pub(crate) fn read(reader: &mut impl Read) -> Result<Self, Error> {
		let mut bytes = [0; 4 + 2 + FILE_ID_LEN];
		reader.read_exact(&mut bytes)?;
		let mut cursor = Cursor::new(&bytes);
		if cursor.array::<4>()? != MAGIC {
			return Err(Error::BadMagic);
		}
		match (cursor.u8()?, cursor.u8()?) {
			(3, 0) => Ok(Self::V3_0),
			(3, 1) => Ok(Self::V3_1),
			(major, minor) => Err(Error::UnsupportedVersion { major, minor }),
		}
	}
}

#[derive(Debug)]
pub(crate) struct ContainerHeader {
	pub(crate) reference: i32,
	pub(crate) start: i32,
	pub(crate) records: u32,
	pub(crate) landmarks: Vec<u32>,
	pub(crate) length: usize,
}

impl ContainerHeader {
	pub(crate) fn is_eof(&self) -> bool {
		self.reference == -1 && self.start == EOF_START && self.records == 0
	}

	pub(crate) fn read(reader: &mut impl Read) -> Result<Option<(Self, Vec<u8>)>, Error> {
		let mut head = vec![0; HEADER_PROBE];
		let got = read_up_to(reader, &mut head)?;
		if got == 0 {
			return Ok(None);
		}
		head.truncate(got);
		loop {
			match Self::parse(&head) {
				Ok((header, consumed)) => {
					let leftover = head.get(consumed..).ok_or(Error::Truncated)?;
					let mut body = Vec::with_capacity(header.length);
					body.extend_from_slice(
						leftover
							.get(..leftover.len().min(header.length))
							.ok_or(Error::Truncated)?,
					);
					let missing = header.length - body.len();
					reader.take(missing as u64).read_to_end(&mut body)?;
					if body.len() != header.length {
						return Err(Error::Truncated);
					}
					return Ok(Some((header, body)));
				}
				Err(Error::Truncated) => {
					let before = head.len();
					head.resize(before + HEADER_PROBE, 0);
					let got = read_up_to(reader, head.get_mut(before..).ok_or(Error::Truncated)?)?;
					if got == 0 {
						return Err(Error::Truncated);
					}
					head.truncate(before + got);
				}
				Err(error) => return Err(error),
			}
		}
	}

	fn parse(bytes: &[u8]) -> Result<(Self, usize), Error> {
		let mut cursor = Cursor::new(bytes);
		let length = check_len(u64::from(cursor.u32_le()?))?;
		let reference = cursor.itf8()?;
		let start = cursor.itf8()?;
		cursor.itf8()?;
		let records = u32::try_from(cursor.itf8()?).map_err(|_| Error::Malformed("container"))?;
		cursor.ltf8()?;
		let _bases = cursor.ltf8()?;
		let _blocks = cursor.itf8()?;
		let landmark_count = usize::try_from(cursor.itf8()?)
			.ok()
			.filter(|count| *count <= MAX_LANDMARKS)
			.ok_or(Error::Malformed("container"))?;
		let mut landmarks = Vec::with_capacity(landmark_count.min(1024));
		for _ in 0..landmark_count {
			landmarks
				.push(u32::try_from(cursor.itf8()?).map_err(|_| Error::Malformed("container"))?);
		}
		let covered = cursor.position();
		let checksum = cursor.u32_le()?;
		let mut crc = Crc::new();
		crc.update(bytes.get(..covered).ok_or(Error::Truncated)?);
		if crc.sum() != checksum {
			return Err(Error::ChecksumMismatch("container header"));
		}
		Ok((
			Self {
				reference,
				start,
				records,
				landmarks,
				length,
			},
			cursor.position(),
		))
	}
}

fn read_up_to(reader: &mut impl Read, buffer: &mut [u8]) -> Result<usize, Error> {
	let mut filled = 0;
	while filled < buffer.len() {
		let slot = buffer.get_mut(filled..).ok_or(Error::Truncated)?;
		match reader.read(slot) {
			Ok(0) => break,
			Ok(count) => filled += count,
			Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
			Err(error) => return Err(error.into()),
		}
	}
	Ok(filled)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BlockKind {
	FileHeader,
	CompressionHeader,
	SliceHeader,
	ExternalData,
	CoreData,
}

impl BlockKind {
	fn from_code(code: u8) -> Option<Self> {
		Some(match code {
			0 => Self::FileHeader,
			1 => Self::CompressionHeader,
			2 => Self::SliceHeader,
			4 => Self::ExternalData,
			5 => Self::CoreData,
			_ => return None,
		})
	}
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Block<'a> {
	pub(crate) kind: BlockKind,
	pub(crate) method: Method,
	pub(crate) content_id: i32,
	raw_len: usize,
	data: &'a [u8],
}

impl<'a> Block<'a> {
	pub(crate) fn read(cursor: &mut Cursor<'a>) -> Result<Self, Error> {
		let started = cursor.position();
		let method = Method::from_code(cursor.u8()?).ok_or(Error::Malformed("block"))?;
		let kind = BlockKind::from_code(cursor.u8()?).ok_or(Error::Malformed("block"))?;
		let content_id = cursor.itf8()?;
		let compressed_len = len_from(cursor.itf8()?)?;
		let raw_len = len_from(cursor.itf8()?)?;
		let data = cursor.take(compressed_len)?;
		let covered = cursor.position() - started;
		let checksum = cursor.u32_le()?;
		let mut crc = Crc::new();
		crc.update(cursor.preceding(started, covered).ok_or(Error::Truncated)?);
		if crc.sum() != checksum {
			return Err(Error::ChecksumMismatch("block"));
		}
		Ok(Self {
			kind,
			method,
			content_id,
			raw_len,
			data,
		})
	}

	pub(crate) fn decode(&self) -> Result<Cow<'a, [u8]>, Error> {
		if self.method == Method::Raw {
			if self.data.len() != self.raw_len {
				return Err(Error::SizeMismatch {
					expected: self.raw_len,
					actual: self.data.len(),
				});
			}
			return Ok(Cow::Borrowed(self.data));
		}
		self.method.decode(self.data, self.raw_len).map(Cow::Owned)
	}
}

fn len_from(value: i32) -> Result<usize, Error> {
	let value = u64::try_from(value).map_err(|_| Error::Malformed("block"))?;
	check_len(value)
}
