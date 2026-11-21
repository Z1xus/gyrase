use super::interleave::deinterleave;
use super::rans_table::{Order0Table, SLOT_MASK, TOTAL_BITS};
use super::stream::Stream;
use super::{boxed_array, check_expansion, check_len};
use crate::cram::cursor::Cursor;
use crate::cram::{Codec, Error};

const CODEC: Codec = Codec::Rans4x8;
const TOTAL: u32 = 1 << TOTAL_BITS;
const SLOTS: usize = 1 << TOTAL_BITS;
const LOWER_BOUND: u32 = 1 << 23;
const WAYS: usize = 4;
const ORDER_0: u8 = 0;
const ORDER_1: u8 = 1;

pub(crate) fn decode(input: &[u8]) -> Result<Vec<u8>, Error> {
	let mut cursor = Cursor::new(input);
	let order = cursor.u8()?;
	let compressed_len = cursor.u32_le()?;
	let raw_len = check_len(cursor.u32_le()?.into())?;
	if usize::try_from(compressed_len) != Ok(cursor.remaining()) {
		return Err(Error::Corrupt(CODEC));
	}
	check_expansion(raw_len, cursor.remaining())?;
	let mut output = vec![0u8; raw_len];
	match order {
		ORDER_0 => decode_order0(&mut cursor, &mut output)?,
		ORDER_1 => decode_order1(&mut cursor, &mut output)?,
		_ => return Err(Error::Corrupt(CODEC)),
	}
	Ok(output)
}

struct Runs {
	symbol: u8,
	remaining: u8,
}

impl Runs {
	fn start(cursor: &mut Cursor<'_>) -> Result<Self, Error> {
		Ok(Self {
			symbol: cursor.u8()?,
			remaining: 0,
		})
	}

	fn advance(&mut self, cursor: &mut Cursor<'_>) -> Result<bool, Error> {
		let next = self.symbol.checked_add(1);
		if self.remaining == 0 && next == Some(cursor.peek()?) {
			self.symbol = cursor.u8()?;
			self.remaining = cursor.u8()?;
		} else if self.remaining > 0 {
			self.remaining -= 1;
			self.symbol = next.ok_or(Error::Corrupt(CODEC))?;
		} else {
			self.symbol = cursor.u8()?;
		}
		Ok(self.symbol != 0)
	}
}

fn read_frequency(cursor: &mut Cursor<'_>) -> Result<u32, Error> {
	let first = cursor.u8()?;
	if first < 0x80 {
		Ok(u32::from(first))
	} else {
		Ok((u32::from(first & 0x7f) << 8) | u32::from(cursor.u8()?))
	}
}

fn read_states(cursor: &mut Cursor<'_>) -> Result<[u32; WAYS], Error> {
	let mut states = [0u32; WAYS];
	for state in &mut states {
		*state = cursor.u32_le()?;
		if *state < LOWER_BOUND {
			return Err(Error::Corrupt(CODEC));
		}
	}
	Ok(states)
}

fn renormalise(mut state: u32, stream: &mut Stream<'_>) -> u32 {
	for _ in 0..2 {
		if state < LOWER_BOUND {
			state = (state << 8) | stream.byte();
		}
	}
	state
}

fn read_order0_table(cursor: &mut Cursor<'_>) -> Result<Order0Table, Error> {
	let mut table = Order0Table::new()?;
	let mut used = 0u32;
	let mut runs = Runs::start(cursor)?;
	loop {
		let frequency = read_frequency(cursor)?;
		used = table
			.fill(runs.symbol, used, frequency)
			.ok_or(Error::Corrupt(CODEC))?;
		if !runs.advance(cursor)? {
			break;
		}
	}
	match TOTAL - used {
		0 => {}
		1 => table.fill_last_slot(),
		_ => return Err(Error::Corrupt(CODEC)),
	}
	Ok(table)
}

fn decode_order0(cursor: &mut Cursor<'_>, output: &mut [u8]) -> Result<(), Error> {
	let table = read_order0_table(cursor)?;
	let mut states = read_states(cursor)?;
	let mut stream = Stream::new(cursor.rest());
	let (groups, remainder) = output.as_chunks_mut::<WAYS>();
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
	Ok(())
}

struct Order1Table {
	symbols: Box<[[u8; SLOTS]; 256]>,
	entries: Box<[[u32; 256]; 256]>,
}

impl Order1Table {
	fn read(cursor: &mut Cursor<'_>) -> Result<Self, Error> {
		let mut table = Self {
			symbols: boxed_array([0; SLOTS])?,
			entries: boxed_array([0; 256])?,
		};
		let mut contexts = Runs::start(cursor)?;
		loop {
			table.read_context(cursor, contexts.symbol)?;
			if !contexts.advance(cursor)? {
				break;
			}
		}
		Ok(table)
	}

	fn read_context(&mut self, cursor: &mut Cursor<'_>, context: u8) -> Result<(), Error> {
		let row = &mut self.symbols[usize::from(context)];
		let entries = &mut self.entries[usize::from(context)];
		let mut used = 0u32;
		let mut runs = Runs::start(cursor)?;
		loop {
			let frequency = match read_frequency(cursor)? {
				0 => TOTAL,
				frequency => frequency,
			};
			let end = used
				.checked_add(frequency)
				.filter(|end| *end <= TOTAL)
				.ok_or(Error::Corrupt(CODEC))?;
			if let Some(slots) = row.get_mut(used as usize..end as usize) {
				slots.fill(runs.symbol);
			}
			entries[usize::from(runs.symbol)] = frequency | (used << 16);
			used = end;
			if !runs.advance(cursor)? {
				break;
			}
		}
		match TOTAL - used {
			0 => Ok(()),
			1 => {
				row[SLOTS - 1] = row[SLOTS - 2];
				Ok(())
			}
			_ => Err(Error::Corrupt(CODEC)),
		}
	}

	fn advance(&self, state: u32, context: u8) -> (u8, u32) {
		let slot = (state & SLOT_MASK) as usize;
		let symbol = self.symbols[usize::from(context)][slot];
		let entry = self.entries[usize::from(context)][usize::from(symbol)];
		let next = (entry & 0xffff)
			.wrapping_mul(state >> TOTAL_BITS)
			.wrapping_add(state & SLOT_MASK)
			.wrapping_sub(entry >> 16);
		(symbol, next)
	}
}

fn decode_order1(cursor: &mut Cursor<'_>, output: &mut [u8]) -> Result<(), Error> {
	let table = Order1Table::read(cursor)?;
	let mut states = read_states(cursor)?;
	let mut stream = Stream::new(cursor.rest());
	let lane_len = output.len() / WAYS;
	let mut lanes = vec![0u8; lane_len * WAYS];
	let mut contexts = [0u8; WAYS];
	for group in lanes.as_chunks_mut::<WAYS>().0 {
		for ((state, context), byte) in states.iter_mut().zip(&mut contexts).zip(group) {
			let (symbol, next) = table.advance(*state, *context);
			*byte = symbol;
			*context = symbol;
			*state = renormalise(next, &mut stream);
		}
	}
	let [.., last_state] = &mut states;
	let [.., last_context] = &mut contexts;
	for byte in output.get_mut(lane_len * WAYS..).unwrap_or_default() {
		let (symbol, next) = table.advance(*last_state, *last_context);
		*byte = symbol;
		*last_context = symbol;
		*last_state = renormalise(next, &mut stream);
	}
	deinterleave::<WAYS>(&lanes, output);
	if stream.is_overrun() {
		return Err(Error::Corrupt(CODEC));
	}
	Ok(())
}
