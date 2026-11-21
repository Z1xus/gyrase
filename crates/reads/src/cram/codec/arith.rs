use super::bzip2;
use super::model::Model;
use super::range_coder::RangeDecoder;
use super::transform::{Packing, decode_stripe};
use super::{check_expansion, check_len, first_reserve};
use crate::cram::cursor::Cursor;
use crate::cram::{Codec, Error};

const CODEC: Codec = Codec::Arith;
const MAX_RUN: usize = 4;
const RUN_CONTEXTS: usize = 258;
const RUN_ESCAPE_CONTEXT: usize = 256;

#[derive(Debug, Clone, Copy)]
struct Flags(u8);

impl Flags {
	const ORDER_MASK: u8 = 0x03;
	const EXTERNAL: u8 = 0x04;
	const STRIPE: u8 = 0x08;
	const NO_SIZE: u8 = 0x10;
	const CAT: u8 = 0x20;
	const RLE: u8 = 0x40;
	const PACK: u8 = 0x80;

	fn has(self, mask: u8) -> bool {
		self.0 & mask != 0
	}
}

pub(crate) fn decode(input: &[u8]) -> Result<Vec<u8>, Error> {
	decode_with_size(input, None)
}

fn decode_with_size(input: &[u8], known_len: Option<usize>) -> Result<Vec<u8>, Error> {
	let mut cursor = Cursor::new(input);
	let flags = Flags(cursor.u8()?);
	if flags.has(Flags::STRIPE) {
		if known_len.is_some() {
			return Err(Error::Corrupt(CODEC));
		}
		return decode_stripe(&mut cursor, CODEC, decode_with_size);
	}
	let len = if flags.has(Flags::NO_SIZE) {
		known_len.ok_or(Error::Corrupt(CODEC))?
	} else {
		check_len(cursor.uint7()?.into())?
	};
	let packing = if flags.has(Flags::PACK) {
		let packing = Packing::read(&mut cursor, CODEC)?;
		let packed_len = check_len(cursor.uint7()?.into())?;
		if packed_len > len {
			return Err(Error::Corrupt(CODEC));
		}
		Some((packing, packed_len))
	} else {
		None
	};
	let entropy_len = packing.as_ref().map_or(len, |(_, packed_len)| *packed_len);
	let payload = cursor.rest();
	let entropy = if entropy_len == 0 {
		Vec::new()
	} else if flags.has(Flags::CAT) {
		payload
			.get(..entropy_len)
			.ok_or(Error::Corrupt(CODEC))?
			.to_vec()
	} else if flags.has(Flags::EXTERNAL) {
		bzip2::decode(payload, entropy_len)?
	} else {
		decode_modelled(flags, payload, entropy_len)?
	};
	match packing {
		Some((packing, _)) => packing.unpack(&entropy, len, CODEC),
		None if entropy.len() == len => Ok(entropy),
		None => Err(Error::Corrupt(CODEC)),
	}
}

fn decode_modelled(flags: Flags, payload: &[u8], len: usize) -> Result<Vec<u8>, Error> {
	let (symbol_count, coded) = payload.split_first().ok_or(Error::Truncated)?;
	check_expansion(len, coded.len())?;
	let symbol_count = match symbol_count {
		0 => 256,
		count => usize::from(*count),
	};
	let mut decoder = RangeDecoder::new(coded)?;
	let output = match (flags.0 & Flags::ORDER_MASK, flags.has(Flags::RLE)) {
		(0, false) => decode_order0(&mut decoder, symbol_count, len)?,
		(1, false) => decode_order1(&mut decoder, symbol_count, len)?,
		(0, true) => decode_runs(&mut decoder, symbol_count, len, false)?,
		(1, true) => decode_runs(&mut decoder, symbol_count, len, true)?,
		_ => return Err(Error::Corrupt(CODEC)),
	};
	if decoder.is_overrun() {
		return Err(Error::Corrupt(CODEC));
	}
	Ok(output)
}

fn decode_byte(model: &mut Model, decoder: &mut RangeDecoder<'_>) -> Result<u8, Error> {
	u8::try_from(model.decode(decoder, CODEC)?).map_err(|_| Error::Corrupt(CODEC))
}

fn decode_order0(
	decoder: &mut RangeDecoder<'_>,
	symbol_count: usize,
	len: usize,
) -> Result<Vec<u8>, Error> {
	let mut model = Model::new(symbol_count);
	let mut output = first_reserve(len);
	for _ in 0..len {
		output.push(decode_byte(&mut model, decoder)?);
	}
	Ok(output)
}

fn decode_order1(
	decoder: &mut RangeDecoder<'_>,
	symbol_count: usize,
	len: usize,
) -> Result<Vec<u8>, Error> {
	let mut models = vec![Model::new(symbol_count); 256];
	let mut output = first_reserve(len);
	let mut previous = 0u8;
	for _ in 0..len {
		let model = models
			.get_mut(usize::from(previous))
			.ok_or(Error::Corrupt(CODEC))?;
		previous = decode_byte(model, decoder)?;
		output.push(previous);
	}
	Ok(output)
}

fn decode_runs(
	decoder: &mut RangeDecoder<'_>,
	symbol_count: usize,
	len: usize,
	order1: bool,
) -> Result<Vec<u8>, Error> {
	let mut literal_models = vec![Model::new(symbol_count); if order1 { 256 } else { 1 }];
	let mut run_models = vec![Model::new(MAX_RUN); RUN_CONTEXTS];
	let mut output = first_reserve(len);
	let mut previous = 0u8;
	while output.len() < len {
		let literal_context = if order1 { usize::from(previous) } else { 0 };
		let literal_model = literal_models
			.get_mut(literal_context)
			.ok_or(Error::Corrupt(CODEC))?;
		let literal = decode_byte(literal_model, decoder)?;
		previous = literal;
		output.push(literal);
		let mut context = usize::from(literal);
		let mut run = 0usize;
		loop {
			let model = run_models.get_mut(context).ok_or(Error::Corrupt(CODEC))?;
			let part = usize::from(model.decode(decoder, CODEC)?);
			context = if context == usize::from(literal) {
				RUN_ESCAPE_CONTEXT
			} else {
				context + usize::from(context < RUN_CONTEXTS - 1)
			};
			run += part;
			if part != MAX_RUN - 1 || run >= len {
				break;
			}
		}
		let copies = run.min(len - output.len());
		output.resize(output.len() + copies, literal);
	}
	Ok(output)
}
