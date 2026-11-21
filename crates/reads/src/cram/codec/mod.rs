mod arith;
mod bzip2;
mod fqzcomp;
mod gzip;
mod interleave;
mod lzma;
mod model;
mod range_coder;
mod rans4x8;
mod rans4x8_encode;
mod rans_nx16;
mod rans_table;
mod stream;
mod tok3;
mod transform;

use crate::cram::Error;

pub(crate) const MAX_RAW_LEN: usize = 1 << 30;
const MAX_EXPANSION: usize = 1 << 19;
const EXPANSION_SLACK: usize = 64;
const FIRST_RESERVE: usize = 1 << 20;

pub(crate) fn check_len(len: u64) -> Result<usize, Error> {
	usize::try_from(len)
		.ok()
		.filter(|len| *len <= MAX_RAW_LEN)
		.ok_or(Error::TooLarge { len })
}

fn check_expansion(len: usize, coded: usize) -> Result<(), Error> {
	if len / MAX_EXPANSION > coded.saturating_add(EXPANSION_SLACK) {
		return Err(Error::TooLarge { len: len as u64 });
	}
	Ok(())
}

fn first_reserve(len: usize) -> Vec<u8> {
	Vec::with_capacity(len.min(FIRST_RESERVE))
}

fn boxed_array<T: Clone, const N: usize>(value: T) -> Result<Box<[T; N]>, Error> {
	Box::<[T; N]>::try_from(vec![value; N].into_boxed_slice())
		.map_err(|_| Error::TooLarge { len: N as u64 })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Method {
	Raw,
	Gzip,
	Bzip2,
	Lzma,
	Rans4x8,
	RansNx16,
	Arith,
	Fqzcomp,
	Tok3,
}

impl Method {
	pub(crate) fn from_code(code: u8) -> Option<Self> {
		Some(match code {
			0 => Self::Raw,
			1 => Self::Gzip,
			2 => Self::Bzip2,
			3 => Self::Lzma,
			4 => Self::Rans4x8,
			5 => Self::RansNx16,
			6 => Self::Arith,
			7 => Self::Fqzcomp,
			8 => Self::Tok3,
			_ => return None,
		})
	}

	pub(crate) fn code(self) -> u8 {
		match self {
			Self::Raw => 0,
			Self::Gzip => 1,
			Self::Bzip2 => 2,
			Self::Lzma => 3,
			Self::Rans4x8 => 4,
			Self::RansNx16 => 5,
			Self::Arith => 6,
			Self::Fqzcomp => 7,
			Self::Tok3 => 8,
		}
	}

	pub(crate) fn decode(self, input: &[u8], raw_len: usize) -> Result<Vec<u8>, Error> {
		let output = match self {
			Self::Raw => input.to_vec(),
			Self::Gzip => return gzip::decode(input, raw_len),
			Self::Bzip2 => return bzip2::decode(input, raw_len),
			Self::Lzma => return lzma::decode(input, raw_len),
			Self::Rans4x8 => rans4x8::decode(input)?,
			Self::RansNx16 => rans_nx16::decode(input)?,
			Self::Arith => arith::decode(input)?,
			Self::Fqzcomp => fqzcomp::decode(input)?,
			Self::Tok3 => tok3::decode(input)?,
		};
		if output.len() == raw_len {
			Ok(output)
		} else {
			Err(Error::SizeMismatch {
				expected: raw_len,
				actual: output.len(),
			})
		}
	}
}

const MIN_COMPRESSED_LEN: usize = 32;

pub(crate) fn compress(raw: &[u8]) -> (Method, Vec<u8>) {
	let mut best = (Method::Raw, raw.to_vec());
	if raw.len() < MIN_COMPRESSED_LEN {
		return best;
	}
	let candidates = [
		(Method::Gzip, gzip::encode(raw)),
		(
			Method::Rans4x8,
			Some(rans4x8_encode::encode(raw, rans4x8_encode::Order::Zero)),
		),
		(
			Method::Rans4x8,
			Some(rans4x8_encode::encode(raw, rans4x8_encode::Order::One)),
		),
	];
	for (method, encoded) in candidates {
		if let Some(encoded) = encoded
			&& encoded.len() < best.1.len()
		{
			best = (method, encoded);
		}
	}
	best
}
