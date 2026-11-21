use libdeflater::{Crc, Decompressor};

use crate::Error;

pub const MAX_BLOCK_LEN: usize = 65536;
pub(crate) const MAX_PAYLOAD_LEN: usize = 65280;
pub(crate) const FIXED_HEADER_LEN: usize = 12;
pub(crate) const FOOTER_LEN: usize = 8;
pub(crate) const EMPTY_BLOCK: [u8; 28] = [
	0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0x1b, 0, 3, 0, 0, 0, 0, 0, 0, 0,
	0, 0,
];
const GZIP_MAGIC: [u8; 3] = [0x1f, 0x8b, 8];
const HEADER_CRC_FLAG: u8 = 2;
const EXTRA_FLAG: u8 = 4;
const NAME_FLAG: u8 = 8;
const COMMENT_FLAG: u8 = 16;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Block {
	pub offset: u64,
	pub compressed_len: u32,
	pub data: Vec<u8>,
}

#[derive(Debug)]
pub(crate) struct Frame {
	pub offset: u64,
	pub bytes: Vec<u8>,
}

pub(crate) fn has_magic(header: &[u8]) -> bool {
	header.len() >= FIXED_HEADER_LEN
		&& header.starts_with(&GZIP_MAGIC)
		&& header[3] & EXTRA_FLAG != 0
}

pub(crate) fn is_bgzf(header: &[u8]) -> bool {
	has_magic(header) && block_size(header).is_some()
}

pub(crate) fn block_len(header: &[u8], offset: u64) -> Result<usize, Error> {
	if !has_magic(header) {
		return Err(Error::NotBgzf { offset });
	}
	let len = block_size(header).ok_or(Error::NotBgzf { offset })?;
	if !(FIXED_HEADER_LEN + FOOTER_LEN..=MAX_BLOCK_LEN).contains(&len) {
		return Err(Error::Corrupt { offset });
	}
	Ok(len)
}

pub(crate) fn extra_len(header: &[u8]) -> usize {
	usize::from(u16::from_le_bytes([header[10], header[11]]))
}

fn payload_start(bytes: &[u8]) -> Option<usize> {
	let flags = bytes[3];
	let mut start = FIXED_HEADER_LEN + extra_len(bytes);
	for flag in [NAME_FLAG, COMMENT_FLAG] {
		if flags & flag != 0 {
			start += bytes.get(start..)?.iter().position(|&byte| byte == 0)? + 1;
		}
	}
	if flags & HEADER_CRC_FLAG != 0 {
		start += 2;
	}
	Some(start)
}

fn block_size(header: &[u8]) -> Option<usize> {
	let extra = header.get(FIXED_HEADER_LEN..FIXED_HEADER_LEN + extra_len(header))?;
	let mut fields = extra;
	while let [first, second, low, high, rest @ ..] = fields {
		let len = usize::from(u16::from_le_bytes([*low, *high]));
		let (value, next) = rest.split_at_checked(len)?;
		if (*first, *second) == (b'B', b'C') {
			let [low, high] = value else { return None };
			return Some(usize::from(u16::from_le_bytes([*low, *high])) + 1);
		}
		fields = next;
	}
	None
}

impl Frame {
	pub(crate) fn inflate(self, decompressor: &mut Decompressor) -> Result<Block, Error> {
		let corrupt = || Error::Corrupt {
			offset: self.offset,
		};
		let bytes = &self.bytes;
		let start = payload_start(bytes).ok_or_else(corrupt)?;
		let footer = bytes
			.len()
			.checked_sub(FOOTER_LEN)
			.filter(|&end| end >= start)
			.ok_or_else(corrupt)?;
		let expected_crc = read_u32(bytes, footer).ok_or_else(corrupt)?;
		let size = read_u32(bytes, footer + 4)
			.and_then(|size| usize::try_from(size).ok())
			.filter(|&size| size <= MAX_BLOCK_LEN)
			.ok_or_else(corrupt)?;
		let mut data = vec![0; size];
		let written = decompressor
			.deflate_decompress(&bytes[start..footer], &mut data)
			.map_err(|_| corrupt())?;
		let mut checksum = Crc::new();
		checksum.update(&data);
		if written != size || checksum.sum() != expected_crc {
			return Err(corrupt());
		}
		let compressed_len = u32::try_from(bytes.len()).map_err(|_| corrupt())?;
		Ok(Block {
			offset: self.offset,
			compressed_len,
			data,
		})
	}
}

fn read_u32(bytes: &[u8], at: usize) -> Option<u32> {
	Some(u32::from_le_bytes(bytes.get(at..at + 4)?.try_into().ok()?))
}
