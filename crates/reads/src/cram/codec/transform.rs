use super::{check_len, first_reserve};
use crate::cram::cursor::Cursor;
use crate::cram::{Codec, Error};

const PACK_ALPHABET_LIMIT: usize = 16;

pub(super) type StreamDecoder = fn(&[u8], Option<usize>) -> Result<Vec<u8>, Error>;

#[derive(Debug)]
pub(super) enum Packing {
	Identity,
	Constant(u8),
	Bits { width: u32, symbols: [u8; 16] },
}

impl Packing {
	pub(super) fn read(cursor: &mut Cursor<'_>, codec: Codec) -> Result<Self, Error> {
		let count = match cursor.u8()? {
			0 => 256,
			count => usize::from(count),
		};
		if count > PACK_ALPHABET_LIMIT {
			return Ok(Self::Identity);
		}
		let mut symbols = [0u8; 16];
		for (slot, symbol) in symbols.iter_mut().zip(cursor.take(count)?) {
			*slot = *symbol;
		}
		match count {
			0 => Err(Error::Corrupt(codec)),
			1 => Ok(Self::Constant(symbols[0])),
			2 => Ok(Self::Bits { width: 1, symbols }),
			3 | 4 => Ok(Self::Bits { width: 2, symbols }),
			_ => Ok(Self::Bits { width: 4, symbols }),
		}
	}

	pub(super) fn unpack(&self, packed: &[u8], len: usize, codec: Codec) -> Result<Vec<u8>, Error> {
		match self {
			Self::Identity => Ok(packed.to_vec()),
			Self::Constant(symbol) => Ok(vec![*symbol; len]),
			Self::Bits { width, symbols } => {
				let per_byte = (8 / width) as usize;
				if len.div_ceil(per_byte) > packed.len() {
					return Err(Error::Corrupt(codec));
				}
				let mask = (1u8 << width) - 1;
				let mut output = first_reserve(len);
				'bytes: for byte in packed {
					for shift in (0..8u32).step_by(*width as usize) {
						if output.len() == len {
							break 'bytes;
						}
						let code = (byte >> shift) & mask;
						output.push(symbols[usize::from(code)]);
					}
				}
				Ok(output)
			}
		}
	}
}

pub(super) fn expand_runs(
	literals: &[u8],
	run_symbols: &[u8],
	run_lengths: &[u8],
	len: usize,
	codec: Codec,
) -> Result<Vec<u8>, Error> {
	let mut repeated = [false; 256];
	for symbol in run_symbols {
		repeated[usize::from(*symbol)] = true;
	}
	let mut runs = Cursor::new(run_lengths);
	let mut output = first_reserve(len);
	for literal in literals {
		let copies = if repeated[usize::from(*literal)] {
			runs.uint7()? as usize + 1
		} else {
			1
		};
		if output.len() + copies > len {
			return Err(Error::Corrupt(codec));
		}
		output.resize(output.len() + copies, *literal);
	}
	Ok(output)
}

fn unstripe(streams: &[Vec<u8>], len: usize) -> Vec<u8> {
	let width = streams.len();
	let mut output = vec![0u8; len];
	for (lane, stream) in streams.iter().enumerate() {
		for (slot, byte) in output.iter_mut().skip(lane).step_by(width).zip(stream) {
			*slot = *byte;
		}
	}
	output
}

pub(super) fn decode_stripe(
	cursor: &mut Cursor<'_>,
	codec: Codec,
	decode_stream: StreamDecoder,
) -> Result<Vec<u8>, Error> {
	let len = check_len(cursor.uint7()?.into())?;
	let count = usize::from(cursor.u8()?);
	if count == 0 {
		return Err(Error::Corrupt(codec));
	}
	let mut compressed_lens = Vec::with_capacity(count);
	for _ in 0..count {
		compressed_lens.push(check_len(cursor.uint7()?.into())?);
	}
	let mut streams = Vec::with_capacity(count);
	for (index, compressed_len) in compressed_lens.into_iter().enumerate() {
		let stream_len = len / count + usize::from(len % count > index);
		let stream = decode_stream(cursor.take(compressed_len)?, Some(stream_len))?;
		if stream.len() != stream_len {
			return Err(Error::Corrupt(codec));
		}
		streams.push(stream);
	}
	Ok(unstripe(&streams, len))
}
