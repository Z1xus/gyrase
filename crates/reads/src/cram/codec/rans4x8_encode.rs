use super::rans_table::{TOTAL, TOTAL_BITS};

const LOWER_BOUND: u32 = 1 << 23;
const WAYS: usize = 4;
const SYMBOLS: usize = 256;
const HEADER_LEN: usize = 9;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Order {
	Zero,
	One,
}

#[derive(Debug, Clone, Copy, Default)]
struct Symbol {
	start: u32,
	frequency: u32,
}

struct Lane {
	state: u32,
}

impl Lane {
	fn new() -> Self {
		Self { state: LOWER_BOUND }
	}

	fn put(&mut self, symbol: Symbol, bytes: &mut Vec<u8>) {
		let limit = ((LOWER_BOUND >> TOTAL_BITS) << 8) * symbol.frequency;
		while self.state >= limit {
			bytes.push(self.state.to_le_bytes()[0]);
			self.state >>= 8;
		}
		self.state = ((self.state / symbol.frequency) << TOTAL_BITS)
			+ (self.state % symbol.frequency)
			+ symbol.start;
	}

	fn flush(&self, bytes: &mut Vec<u8>) {
		bytes.extend(self.state.to_le_bytes().iter().rev());
	}
}

fn normalise(counts: &[u32; SYMBOLS]) -> [u32; SYMBOLS] {
	let total: u64 = counts.iter().map(|count| u64::from(*count)).sum();
	let mut frequencies = [0u32; SYMBOLS];
	if total == 0 {
		return frequencies;
	}
	let mut assigned = 0u32;
	for (frequency, count) in frequencies.iter_mut().zip(counts) {
		if *count > 0 {
			let scaled = u64::from(*count) * u64::from(TOTAL) / total;
			*frequency = u32::try_from(scaled).unwrap_or(TOTAL).max(1);
			assigned += *frequency;
		}
	}
	while assigned != TOTAL {
		let Some((index, _)) = frequencies
			.iter()
			.enumerate()
			.max_by_key(|(_, frequency)| **frequency)
		else {
			break;
		};
		let Some(frequency) = frequencies.get_mut(index) else {
			break;
		};
		if assigned < TOTAL {
			*frequency += TOTAL - assigned;
			assigned = TOTAL;
		} else {
			let take = (assigned - TOTAL).min(*frequency - 1);
			*frequency -= take;
			assigned -= take;
		}
	}
	frequencies
}

fn symbols(frequencies: &[u32; SYMBOLS]) -> [Symbol; SYMBOLS] {
	let mut table = [Symbol::default(); SYMBOLS];
	let mut start = 0;
	for (symbol, frequency) in table.iter_mut().zip(frequencies) {
		*symbol = Symbol {
			start,
			frequency: *frequency,
		};
		start += frequency;
	}
	table
}

fn write_frequency(out: &mut Vec<u8>, frequency: u32) {
	if frequency < 0x80 {
		out.push(frequency.to_le_bytes()[0]);
	} else {
		let [low, high, ..] = frequency.to_le_bytes();
		out.push(high | 0x80);
		out.push(low);
	}
}

fn write_alphabet(
	out: &mut Vec<u8>,
	present: impl Fn(usize) -> bool,
	mut item: impl FnMut(&mut Vec<u8>, usize),
) {
	let mut run = 0usize;
	for symbol in 0..SYMBOLS {
		if !present(symbol) {
			continue;
		}
		if run > 0 {
			run -= 1;
		} else {
			out.push(symbol.to_le_bytes()[0]);
			if symbol > 0 && present(symbol - 1) {
				let mut end = symbol + 1;
				while end < SYMBOLS && present(end) {
					end += 1;
				}
				run = end - symbol - 1;
				out.push(run.to_le_bytes()[0]);
			}
		}
		item(out, symbol);
	}
	out.push(0);
}

fn finish(mut body: Vec<u8>, lanes: &[Lane; WAYS], mut stream: Vec<u8>) -> Vec<u8> {
	for lane in lanes.iter().rev() {
		lane.flush(&mut stream);
	}
	stream.reverse();
	body.extend_from_slice(&stream);
	body
}

pub(super) fn encode(input: &[u8], order: Order) -> Vec<u8> {
	let (code, body) = match order {
		Order::Zero => (0, encode_order0(input)),
		Order::One => (1, encode_order1(input)),
	};
	let mut out = Vec::with_capacity(HEADER_LEN + body.len());
	out.push(code);
	out.extend_from_slice(&u32::try_from(body.len()).unwrap_or(u32::MAX).to_le_bytes());
	out.extend_from_slice(&u32::try_from(input.len()).unwrap_or(u32::MAX).to_le_bytes());
	out.extend_from_slice(&body);
	out
}

fn encode_order0(input: &[u8]) -> Vec<u8> {
	let mut counts = [0u32; SYMBOLS];
	for byte in input {
		if let Some(count) = counts.get_mut(usize::from(*byte)) {
			*count += 1;
		}
	}
	let frequencies = normalise(&counts);
	let table = symbols(&frequencies);
	let mut body = Vec::new();
	write_alphabet(
		&mut body,
		|symbol| {
			frequencies
				.get(symbol)
				.is_some_and(|frequency| *frequency > 0)
		},
		|out, symbol| write_frequency(out, frequencies.get(symbol).copied().unwrap_or(0)),
	);
	let look = |byte: u8| table.get(usize::from(byte)).copied().unwrap_or_default();
	let mut lanes = [Lane::new(), Lane::new(), Lane::new(), Lane::new()];
	let mut bytes = Vec::with_capacity(input.len() / 2 + 16);
	let (groups, tail) = input.as_chunks::<WAYS>();
	for (lane, byte) in lanes.iter_mut().zip(tail).rev() {
		lane.put(look(*byte), &mut bytes);
	}
	let [first, second, third, fourth] = &mut lanes;
	for [a, b, c, d] in groups.iter().rev() {
		fourth.put(look(*d), &mut bytes);
		third.put(look(*c), &mut bytes);
		second.put(look(*b), &mut bytes);
		first.put(look(*a), &mut bytes);
	}
	finish(body, &lanes, bytes)
}

fn encode_order1(input: &[u8]) -> Vec<u8> {
	let lane_len = input.len() / WAYS;
	let mut counts = vec![[0u32; SYMBOLS]; SYMBOLS];
	for_each_transition(input, lane_len, |context, byte| {
		if let Some(count) = counts
			.get_mut(usize::from(context))
			.and_then(|row| row.get_mut(usize::from(byte)))
		{
			*count += 1;
		}
	});
	let frequencies: Vec<[u32; SYMBOLS]> = counts.iter().map(normalise).collect();
	let tables: Vec<[Symbol; SYMBOLS]> = frequencies.iter().map(symbols).collect();
	let mut body = Vec::new();
	write_alphabet(
		&mut body,
		|context| {
			counts
				.get(context)
				.is_some_and(|row| row.iter().any(|count| *count > 0))
		},
		|out, context| {
			let row = frequencies.get(context);
			write_alphabet(
				out,
				|symbol| row.is_some_and(|row| row.get(symbol).is_some_and(|f| *f > 0)),
				|out, symbol| {
					let frequency = row.and_then(|row| row.get(symbol)).copied().unwrap_or(0);
					write_frequency(out, frequency);
				},
			);
		},
	);
	let look = |context: u8, byte: u8| {
		tables
			.get(usize::from(context))
			.and_then(|row| row.get(usize::from(byte)))
			.copied()
			.unwrap_or_default()
	};
	let mut lanes = [Lane::new(), Lane::new(), Lane::new(), Lane::new()];
	let mut bytes = Vec::with_capacity(input.len() / 2 + 16);
	let tail_start = lane_len * WAYS;
	let [.., last] = &mut lanes;
	for index in (tail_start..input.len()).rev() {
		let context = before(input, index);
		if let Some(byte) = input.get(index) {
			last.put(look(context, *byte), &mut bytes);
		}
	}
	for step in (0..lane_len).rev() {
		for (lane_index, lane) in lanes.iter_mut().enumerate().rev() {
			let index = lane_index * lane_len + step;
			let context = if step == 0 { 0 } else { before(input, index) };
			if let Some(byte) = input.get(index) {
				lane.put(look(context, *byte), &mut bytes);
			}
		}
	}
	finish(body, &lanes, bytes)
}

fn before(input: &[u8], index: usize) -> u8 {
	index
		.checked_sub(1)
		.and_then(|previous| input.get(previous))
		.copied()
		.unwrap_or(0)
}

fn for_each_transition(input: &[u8], lane_len: usize, mut visit: impl FnMut(u8, u8)) {
	let tail_start = lane_len * WAYS;
	for lane in 0..WAYS {
		let mut context = 0u8;
		for index in lane * lane_len..(lane + 1) * lane_len {
			let byte = input.get(index).copied().unwrap_or(0);
			visit(context, byte);
			context = byte;
		}
	}
	let mut context = if lane_len == 0 {
		0
	} else {
		before(input, tail_start)
	};
	for byte in input.get(tail_start..).unwrap_or_default() {
		visit(context, *byte);
		context = *byte;
	}
}
