use libdeflater::Crc;

use crate::cram::Error;
use crate::cram::codec::{self, Method};

pub(super) const EOF_START: i32 = 4_542_278;

#[derive(Debug, Clone, Copy)]
pub(super) enum BlockKind {
	FileHeader,
	CompressionHeader,
	SliceHeader,
	ExternalData,
	CoreData,
}

impl BlockKind {
	fn code(self) -> u8 {
		match self {
			Self::FileHeader => 0,
			Self::CompressionHeader => 1,
			Self::SliceHeader => 2,
			Self::ExternalData => 4,
			Self::CoreData => 5,
		}
	}
}

pub(super) fn int(value: usize) -> Result<i32, Error> {
	i32::try_from(value).map_err(|_| Error::IntegerOverflow)
}

pub(super) fn put_length(out: &mut Vec<u8>, value: usize) -> Result<(), Error> {
	put_itf8(out, int(value)?);
	Ok(())
}

pub(super) fn put_itf8(out: &mut Vec<u8>, value: i32) {
	let bits = value.cast_unsigned();
	let [b0, b1, b2, b3] = bits.to_be_bytes();
	match bits {
		0..0x80 => out.push(b3),
		0x80..0x4000 => out.extend_from_slice(&[0x80 | b2, b3]),
		0x4000..0x20_0000 => out.extend_from_slice(&[0xc0 | b1, b2, b3]),
		0x20_0000..0x1000_0000 => out.extend_from_slice(&[0xe0 | b0, b1, b2, b3]),
		_ => out.extend_from_slice(&[
			0xf0 | (b0 >> 4),
			(b0 << 4) | (b1 >> 4),
			(b1 << 4) | (b2 >> 4),
			(b2 << 4) | (b3 >> 4),
			b3 & 0x0f,
		]),
	}
}

pub(super) fn put_ltf8(out: &mut Vec<u8>, value: u64) {
	let bytes = value.to_be_bytes();
	let extra = match value {
		0..0x80 => 0,
		0x80..0x4000 => 1,
		0x4000..0x20_0000 => 2,
		0x20_0000..0x1000_0000 => 3,
		0x1000_0000..0x8_0000_0000 => 4,
		0x8_0000_0000..0x400_0000_0000 => 5,
		0x400_0000_0000..0x2_0000_0000_0000 => 6,
		0x2_0000_0000_0000..0x100_0000_0000_0000 => 7,
		_ => 8,
	};
	let prefix: u8 = match extra {
		0 => 0,
		8 => 0xff,
		count => 0xffu8 << (8 - count),
	};
	let first = bytes.get(7 - extra).copied().unwrap_or(0);
	if extra == 8 {
		out.push(prefix);
	} else {
		out.push(prefix | first);
	}
	out.extend(bytes.iter().skip(8 - extra));
}

pub(super) fn crc32(bytes: &[u8]) -> u32 {
	let mut crc = Crc::new();
	crc.update(bytes);
	crc.sum()
}

pub(super) fn block(
	out: &mut Vec<u8>,
	kind: BlockKind,
	content_id: i32,
	raw: &[u8],
) -> Result<(), Error> {
	let (method, data) = codec::compress(raw);
	write_block(out, kind, content_id, method, &data, raw.len())
}

pub(super) fn raw_block(
	out: &mut Vec<u8>,
	kind: BlockKind,
	content_id: i32,
	raw: &[u8],
) -> Result<(), Error> {
	write_block(out, kind, content_id, Method::Raw, raw, raw.len())
}

fn write_block(
	out: &mut Vec<u8>,
	kind: BlockKind,
	content_id: i32,
	method: Method,
	data: &[u8],
	raw_len: usize,
) -> Result<(), Error> {
	let started = out.len();
	out.push(method.code());
	out.push(kind.code());
	put_itf8(out, content_id);
	put_length(out, data.len())?;
	put_length(out, raw_len)?;
	out.extend_from_slice(data);
	let sum = out.get(started..).map(crc32).unwrap_or_default();
	out.extend_from_slice(&sum.to_le_bytes());
	Ok(())
}

#[derive(Debug)]
pub(super) struct ContainerFields {
	pub(super) reference: i32,
	pub(super) start: i32,
	pub(super) span: i32,
	pub(super) records: u32,
	pub(super) counter: u64,
	pub(super) bases: u64,
	pub(super) blocks: i32,
	pub(super) landmarks: Vec<u32>,
}

pub(super) fn container_header(
	fields: &ContainerFields,
	body_len: usize,
) -> Result<Vec<u8>, Error> {
	let mut out = Vec::new();
	let body_len = u32::try_from(body_len).map_err(|_| Error::IntegerOverflow)?;
	out.extend_from_slice(&body_len.to_le_bytes());
	put_itf8(&mut out, fields.reference);
	put_itf8(&mut out, fields.start);
	put_itf8(&mut out, fields.span);
	put_itf8(
		&mut out,
		i32::try_from(fields.records).map_err(|_| Error::IntegerOverflow)?,
	);
	put_ltf8(&mut out, fields.counter);
	put_ltf8(&mut out, fields.bases);
	put_itf8(&mut out, fields.blocks);
	put_length(&mut out, fields.landmarks.len())?;
	for landmark in &fields.landmarks {
		put_itf8(
			&mut out,
			i32::try_from(*landmark).map_err(|_| Error::IntegerOverflow)?,
		);
	}
	let sum = crc32(&out);
	out.extend_from_slice(&sum.to_le_bytes());
	Ok(out)
}
