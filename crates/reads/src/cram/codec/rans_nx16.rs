use super::interleave::deinterleave;
use super::rans_table::{Order0Table, SLOTS, TOTAL_BITS};
use super::stream::Stream;
use super::transform::{Packing, decode_stripe, expand_runs};
use super::{boxed_array, check_expansion, check_len};
use crate::cram::cursor::Cursor;
use crate::cram::{Codec, Error};

const CODEC: Codec = Codec::RansNx16;
const LOWER_BOUND: u32 = 1 << 15;
const FAST_ORDER1_BITS: u32 = 10;
const ALPHABET: usize = 256;

#[derive(Debug, Clone, Copy)]
struct Flags(u8);

impl Flags {
	const ORDER1: u8 = 0x01;
	const WIDE: u8 = 0x04;
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
	let packed_len = packing.as_ref().map_or(len, |(_, packed_len)| *packed_len);
	let runs = if flags.has(Flags::RLE) {
		Some(RunHeader::read(
			&mut cursor,
			packed_len,
			flags.has(Flags::WIDE),
		)?)
	} else {
		None
	};
	let entropy_len = runs.as_ref().map_or(packed_len, |runs| runs.literals_len);
	let payload = cursor.rest();
	let entropy = if entropy_len == 0 {
		Vec::new()
	} else if flags.has(Flags::CAT) {
		payload
			.get(..entropy_len)
			.ok_or(Error::Corrupt(CODEC))?
			.to_vec()
	} else {
		decode_entropy(flags, payload, entropy_len)?
	};
	let unrun = match &runs {
		Some(runs) => expand_runs(&entropy, &runs.symbols, &runs.lengths, packed_len, CODEC)?,
		None => entropy,
	};
	match packing {
		Some((packing, _)) => packing.unpack(&unrun, len, CODEC),
		None if unrun.len() == len => Ok(unrun),
		None => Err(Error::Corrupt(CODEC)),
	}
}

struct RunHeader {
	literals_len: usize,
	symbols: Vec<u8>,
	lengths: Vec<u8>,
}

impl RunHeader {
	fn read(cursor: &mut Cursor<'_>, packed_len: usize, wide: bool) -> Result<Self, Error> {
		let meta_code = cursor.uint7()?;
		let literals_len = check_len(cursor.uint7()?.into())?;
		if literals_len > packed_len {
			return Err(Error::Corrupt(CODEC));
		}
		let meta_len = check_len(u64::from(meta_code / 2))?;
		let meta = if meta_code & 1 == 1 {
			cursor.take(meta_len)?.to_vec()
		} else {
			let compressed_len = check_len(cursor.uint7()?.into())?;
			let compressed = cursor.take(compressed_len)?;
			if wide {
				decode_order0::<32>(compressed, meta_len)?
			} else {
				decode_order0::<4>(compressed, meta_len)?
			}
		};
		let symbol_count = match meta.first().ok_or(Error::Corrupt(CODEC))? {
			0 => ALPHABET,
			count => usize::from(*count),
		};
		let symbols = meta
			.get(1..=symbol_count)
			.ok_or(Error::Corrupt(CODEC))?
			.to_vec();
		let lengths = meta
			.get(symbol_count + 1..)
			.ok_or(Error::Corrupt(CODEC))?
			.to_vec();
		Ok(Self {
			literals_len,
			symbols,
			lengths,
		})
	}
}

fn decode_entropy(flags: Flags, payload: &[u8], len: usize) -> Result<Vec<u8>, Error> {
	match (flags.has(Flags::ORDER1), flags.has(Flags::WIDE)) {
		(false, false) => decode_order0::<4>(payload, len),
		(false, true) => decode_order0::<32>(payload, len),
		(true, false) => decode_order1::<4>(payload, len),
		(true, true) => decode_order1::<32>(payload, len),
	}
}

fn read_alphabet(cursor: &mut Cursor<'_>) -> Result<[bool; ALPHABET], Error> {
	let mut present = [false; ALPHABET];
	let mut symbol = cursor.u8()?;
	let mut remaining = 0u8;
	loop {
		present[usize::from(symbol)] = true;
		let next = symbol.checked_add(1);
		if remaining == 0 && next == Some(cursor.peek()?) {
			symbol = cursor.u8()?;
			remaining = cursor.u8()?;
		} else if remaining > 0 {
			remaining -= 1;
			symbol = next.ok_or(Error::Corrupt(CODEC))?;
		} else {
			symbol = cursor.u8()?;
		}
		if symbol == 0 {
			return Ok(present);
		}
	}
}

fn scale_shift(total: u64, target: u64) -> u32 {
	if total == 0 {
		return 0;
	}
	let mut scaled = total;
	let mut shift = 0;
	while scaled < target {
		scaled *= 2;
		shift += 1;
	}
	shift
}

fn read_order0_table(cursor: &mut Cursor<'_>) -> Result<Order0Table, Error> {
	let present = read_alphabet(cursor)?;
	let mut frequencies = [0u32; ALPHABET];
	let mut total = 0u64;
	for (frequency, _) in frequencies
		.iter_mut()
		.zip(present)
		.filter(|(_, present)| *present)
	{
		*frequency = cursor.uint7()?;
		total += u64::from(*frequency);
	}
	let shift = scale_shift(total, SLOTS as u64);
	let mut table = Order0Table::new()?;
	let mut used = 0u32;
	for (symbol, frequency) in (0u8..=255).zip(frequencies) {
		if frequency != 0 {
			let scaled = frequency.checked_shl(shift).ok_or(Error::Corrupt(CODEC))?;
			used = table
				.fill(symbol, used, scaled)
				.ok_or(Error::Corrupt(CODEC))?;
		}
	}
	if used as usize != SLOTS {
		return Err(Error::Corrupt(CODEC));
	}
	Ok(table)
}

fn read_states<const N: usize>(cursor: &mut Cursor<'_>) -> Result<[u32; N], Error> {
	let mut states = [0u32; N];
	for state in &mut states {
		*state = cursor.u32_le()?;
		if *state < LOWER_BOUND {
			return Err(Error::Corrupt(CODEC));
		}
	}
	Ok(states)
}

fn renormalise(state: u32, stream: &mut Stream<'_>) -> u32 {
	if state < LOWER_BOUND {
		(state << 16) | stream.word()
	} else {
		state
	}
}

pub(super) fn decode_order0<const N: usize>(payload: &[u8], len: usize) -> Result<Vec<u8>, Error> {
	check_expansion(len, payload.len())?;
	let mut cursor = Cursor::new(payload);
	let table = read_order0_table(&mut cursor)?;
	let mut states = read_states::<N>(&mut cursor)?;
	let mut stream = Stream::new(cursor.rest());
	let mut output = vec![0u8; len];
	let (groups, remainder) = output.as_chunks_mut::<N>();
	for group in groups {
		for (state, byte) in states.iter_mut().zip(group) {
			let (symbol, next) = table.advance(*state);
			*byte = symbol;
			*state = renormalise(next, &mut stream);
		}
	}
	for (state, byte) in states.iter().zip(remainder) {
		*byte = table.peek(*state);
	}
	if stream.is_overrun() {
		return Err(Error::Corrupt(CODEC));
	}
	Ok(output)
}

struct Order1Table {
	shift: u32,
	symbols: Box<[[u8; SLOTS]; ALPHABET]>,
	entries: Box<[[u32; ALPHABET]; ALPHABET]>,
}

impl Order1Table {
	fn read(cursor: &mut Cursor<'_>) -> Result<Self, Error> {
		let header = cursor.u8()?;
		let shift = u32::from(header >> 4);
		if shift != TOTAL_BITS && shift != FAST_ORDER1_BITS {
			return Err(Error::Corrupt(CODEC));
		}
		if header & 1 == 0 {
			return Self::read_rows(cursor, shift);
		}
		let raw_len = check_len(cursor.uint7()?.into())?;
		let compressed_len = check_len(cursor.uint7()?.into())?;
		let rows = decode_order0::<4>(cursor.take(compressed_len)?, raw_len)?;
		Self::read_rows(&mut Cursor::new(&rows), shift)
	}

	fn read_rows(cursor: &mut Cursor<'_>, shift: u32) -> Result<Self, Error> {
		let mut table = Self {
			shift,
			symbols: boxed_array([0; SLOTS])?,
			entries: boxed_array([0; ALPHABET])?,
		};
		let present = read_alphabet(cursor)?;
		for (context, _) in (0u8..=255).zip(present).filter(|(_, present)| *present) {
			table.read_context(cursor, context, &present)?;
		}
		Ok(table)
	}

	fn read_context(
		&mut self,
		cursor: &mut Cursor<'_>,
		context: u8,
		present: &[bool; ALPHABET],
	) -> Result<(), Error> {
		let mut frequencies = [0u32; ALPHABET];
		let mut zero_run = 0u32;
		let mut total = 0u64;
		for (frequency, _) in frequencies.iter_mut().zip(present).filter(|(_, p)| **p) {
			if zero_run > 0 {
				zero_run -= 1;
				continue;
			}
			*frequency = cursor.uint7()?;
			if *frequency == 0 {
				zero_run = u32::from(cursor.u8()?);
			}
			total += u64::from(*frequency);
		}
		if total == 0 {
			return Ok(());
		}
		let target = 1u64 << self.shift;
		let shift = scale_shift(total, target);
		let row = &mut self.symbols[usize::from(context)];
		let entries = &mut self.entries[usize::from(context)];
		let mut used = 0u32;
		for (symbol, frequency) in (0u8..=255).zip(frequencies) {
			if frequency == 0 {
				continue;
			}
			let scaled = frequency.checked_shl(shift).ok_or(Error::Corrupt(CODEC))?;
			let end = used
				.checked_add(scaled)
				.filter(|end| u64::from(*end) <= target)
				.ok_or(Error::Corrupt(CODEC))?;
			row.get_mut(used as usize..end as usize)
				.ok_or(Error::Corrupt(CODEC))?
				.fill(symbol);
			entries[usize::from(symbol)] = scaled | (used << 16);
			used = end;
		}
		if u64::from(used) != target {
			return Err(Error::Corrupt(CODEC));
		}
		Ok(())
	}

	fn advance(&self, state: u32, context: u8) -> (u8, u32) {
		let slot = state & ((1 << self.shift) - 1);
		let symbol = self.symbols[usize::from(context)][slot as usize & (SLOTS - 1)];
		let entry = self.entries[usize::from(context)][usize::from(symbol)];
		let next = (entry & 0xffff)
			.wrapping_mul(state >> self.shift)
			.wrapping_add(slot)
			.wrapping_sub(entry >> 16);
		(symbol, next)
	}
}

pub(super) fn decode_order1<const N: usize>(payload: &[u8], len: usize) -> Result<Vec<u8>, Error> {
	check_expansion(len, payload.len())?;
	let mut cursor = Cursor::new(payload);
	let table = Order1Table::read(&mut cursor)?;
	let mut states = read_states::<N>(&mut cursor)?;
	let mut stream = Stream::new(cursor.rest());
	let lane_len = len / N;
	let mut lanes = vec![0u8; lane_len * N];
	let mut contexts = [0u8; N];
	for group in lanes.as_chunks_mut::<N>().0 {
		for ((state, context), byte) in states.iter_mut().zip(&mut contexts).zip(group) {
			let (symbol, next) = table.advance(*state, *context);
			*byte = symbol;
			*context = symbol;
			*state = renormalise(next, &mut stream);
		}
	}
	let mut output = vec![0u8; len];
	if let (Some(last_state), Some(last_context)) = (states.last_mut(), contexts.last_mut()) {
		for byte in output.get_mut(lane_len * N..).unwrap_or_default() {
			let (symbol, next) = table.advance(*last_state, *last_context);
			*byte = symbol;
			*last_context = symbol;
			*last_state = renormalise(next, &mut stream);
		}
	}
	deinterleave::<N>(&lanes, &mut output);
	if stream.is_overrun() {
		return Err(Error::Corrupt(CODEC));
	}
	Ok(output)
}
