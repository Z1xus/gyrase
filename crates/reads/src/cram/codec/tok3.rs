use std::rc::Rc;

use super::{arith, check_len, first_reserve, rans_nx16};
use crate::cram::cursor::Cursor;
use crate::cram::{Codec, Error};

const CODEC: Codec = Codec::Tok3;
const MAX_TOKENS: usize = 128;
const KINDS_PER_TOKEN: usize = 16;
const MAX_NAMES: u32 = 10_000_000;
const OUTPUT_SLACK: usize = 1024;
const MAX_FIXED_DIGITS: usize = 9;
const NEW_TOKEN: u8 = 0x80;
const COPY_STREAM: u8 = 0x40;
const KIND_MASK: u8 = 0x0f;

fn corrupt() -> Error {
	Error::Corrupt(CODEC)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
	Type = 0,
	Alpha = 1,
	Char = 2,
	Digits0 = 3,
	DigitsZeroLen = 4,
	Dup = 5,
	Diff = 6,
	Digits = 7,
	Delta = 8,
	Delta0 = 9,
	Match = 10,
	Nop = 11,
	End = 12,
}

impl Kind {
	fn from_byte(byte: u8) -> Option<Self> {
		Some(match byte {
			0 => Self::Type,
			1 => Self::Alpha,
			2 => Self::Char,
			3 => Self::Digits0,
			4 => Self::DigitsZeroLen,
			5 => Self::Dup,
			6 => Self::Diff,
			7 => Self::Digits,
			8 => Self::Delta,
			9 => Self::Delta0,
			10 => Self::Match,
			11 => Self::Nop,
			12 => Self::End,
			_ => return None,
		})
	}
}

#[derive(Debug, Default, Clone)]
struct Stream {
	data: Option<Rc<Vec<u8>>>,
	implicit: usize,
	position: usize,
}

impl Stream {
	fn present(data: Vec<u8>) -> Self {
		Self {
			data: Some(Rc::new(data)),
			implicit: 0,
			position: 0,
		}
	}

	fn next_byte(&mut self) -> Option<u8> {
		let data = self.data.as_ref()?;
		let byte = match data.get(self.position) {
			Some(byte) => *byte,
			None if self.position < data.len() + self.implicit => Kind::Match as u8,
			None => return None,
		};
		self.position += 1;
		Some(byte)
	}

	fn next_u32(&mut self) -> Option<u32> {
		let bytes = self.data.as_ref()?.get(self.position..self.position + 4)?;
		let value = u32::from_le_bytes(bytes.try_into().ok()?);
		self.position += 4;
		Some(value)
	}

	fn next_text(&mut self) -> Option<&[u8]> {
		let data = self.data.as_ref()?;
		let rest = data.get(self.position..)?;
		let length = rest.iter().position(|byte| *byte == 0)?;
		self.position += length + 1;
		rest.get(..length)
	}
}

struct Streams {
	streams: Vec<Stream>,
	columns: usize,
}

impl Streams {
	fn get(&mut self, column: usize, kind: Kind) -> Option<&mut Stream> {
		self.streams
			.get_mut(column * KINDS_PER_TOKEN + kind as usize)
	}
}

fn read_streams(input: &[u8], name_count: usize, use_arith: bool) -> Result<Streams, Error> {
	let mut cursor = Cursor::new(input);
	let mut streams = Streams {
		streams: vec![Stream::default(); MAX_TOKENS * KINDS_PER_TOKEN],
		columns: 0,
	};
	let mut column: Option<usize> = None;
	while !cursor.is_empty() {
		let flags = cursor.u8()?;
		let copy_from = if flags & COPY_STREAM != 0 {
			let high = usize::from(cursor.u8()?);
			Some((high << 4) + usize::from(cursor.u8()?))
		} else {
			None
		};
		if flags & NEW_TOKEN != 0 {
			let next = column.map_or(0, |column| column + 1);
			if next >= MAX_TOKENS {
				return Err(corrupt());
			}
			streams.columns = next + 1;
			column = Some(next);
			let first = next * KINDS_PER_TOKEN;
			for stream in streams
				.streams
				.get_mut(first..first + KINDS_PER_TOKEN)
				.ok_or_else(corrupt)?
			{
				*stream = Stream::default();
			}
			if flags & KIND_MASK != 0 {
				*streams.streams.get_mut(first).ok_or_else(corrupt)? = Stream {
					implicit: name_count - 1,
					..Stream::present(vec![flags & KIND_MASK])
				};
			}
		}
		let column = column.ok_or_else(corrupt)?;
		let index = column * KINDS_PER_TOKEN + usize::from(flags & KIND_MASK);
		let stream = if let Some(source) = copy_from {
			if source >= index {
				return Err(corrupt());
			}
			let source = streams.streams.get(source).ok_or_else(corrupt)?;
			Stream {
				position: 0,
				..source.clone()
			}
		} else {
			let compressed_len = check_len(cursor.uint7()?.into())?;
			let compressed = cursor.take(compressed_len)?;
			Stream::present(if use_arith {
				arith::decode(compressed)?
			} else {
				rans_nx16::decode(compressed)?
			})
		};
		*streams.streams.get_mut(index).ok_or_else(corrupt)? = stream;
	}
	Ok(streams)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TokenKind {
	Nop,
	Char,
	Alpha,
	Digits0,
	Digits,
	End,
}

#[derive(Debug, Clone, Copy)]
struct Token {
	kind: TokenKind,
	number: u32,
	text: usize,
}

const EMPTY_TOKEN: Token = Token {
	kind: TokenKind::Nop,
	number: 0,
	text: 0,
};

#[derive(Debug, Clone, Copy)]
struct Name {
	text_start: usize,
	text_len: usize,
	token_start: usize,
	token_count: usize,
}

pub(crate) fn decode(input: &[u8]) -> Result<Vec<u8>, Error> {
	let mut cursor = Cursor::new(input);
	let declared_len = check_len(cursor.u32_le()?.into())?;
	let name_count = cursor.u32_le()?;
	let use_arith = cursor.u8()? != 0;
	if name_count == 0 || name_count > MAX_NAMES || name_count as usize > declared_len {
		return Err(corrupt());
	}
	let mut streams = read_streams(cursor.rest(), name_count as usize, use_arith)?;
	let mut decoder = NameDecoder {
		limit: declared_len + OUTPUT_SLACK,
		name_count: name_count as usize,
		output: first_reserve(declared_len + OUTPUT_SLACK),
		names: Vec::new(),
		tokens: Vec::new(),
		current: Vec::new(),
	};
	while decoder.next_name(&mut streams)? {}
	Ok(decoder.output)
}

struct NameDecoder {
	limit: usize,
	name_count: usize,
	output: Vec<u8>,
	names: Vec<Name>,
	tokens: Vec<Token>,
	current: Vec<Token>,
}

impl NameDecoder {
	fn next_name(&mut self, streams: &mut Streams) -> Result<bool, Error> {
		let index = self.names.len();
		if index > self.name_count {
			return Err(corrupt());
		}
		let Some(kind_byte) = streams.get(0, Kind::Type).and_then(Stream::next_byte) else {
			return Ok(false);
		};
		if usize::from(kind_byte) >= streams.columns * KINDS_PER_TOKEN {
			return Ok(false);
		}
		let kind = Kind::from_byte(kind_byte).ok_or_else(corrupt)?;
		let distance = streams
			.get(0, kind)
			.and_then(Stream::next_u32)
			.ok_or_else(corrupt)? as usize;
		let previous_index = index.checked_sub(distance).ok_or_else(corrupt)?;
		let previous = self.names.get(previous_index).copied();
		if kind == Kind::Dup {
			self.duplicate(previous.ok_or_else(corrupt)?)?;
		} else {
			self.tokenised(streams, previous)?;
		}
		Ok(true)
	}

	fn reserve(&self, extra: usize) -> Result<(), Error> {
		if self.output.len() + extra > self.limit {
			return Err(corrupt());
		}
		Ok(())
	}

	fn duplicate(&mut self, previous: Name) -> Result<(), Error> {
		let text_start = self.output.len();
		self.reserve(previous.text_len + 1)?;
		self.output
			.extend_from_within(previous.text_start..previous.text_start + previous.text_len);
		self.output.push(0);
		let token_start = self.tokens.len();
		self.tokens
			.extend_from_within(previous.token_start..previous.token_start + previous.token_count);
		self.names.push(Name {
			text_start,
			text_len: previous.text_len,
			token_start,
			token_count: previous.token_count,
		});
		Ok(())
	}

	fn previous_token(&self, previous: Option<Name>, column: usize) -> Result<Token, Error> {
		let previous = previous.ok_or_else(corrupt)?;
		if column >= previous.token_count {
			return Err(corrupt());
		}
		self.tokens
			.get(previous.token_start + column)
			.copied()
			.ok_or_else(corrupt)
	}

	fn tokenised(&mut self, streams: &mut Streams, previous: Option<Name>) -> Result<(), Error> {
		let text_start = self.output.len();
		self.current.clear();
		self.current.push(EMPTY_TOKEN);
		for column in 1..MAX_TOKENS.min(streams.columns) {
			let kind = streams
				.get(column, Kind::Type)
				.and_then(Stream::next_byte)
				.and_then(Kind::from_byte);
			let token = match kind {
				Some(Kind::Char) => self.char_token(streams, column)?,
				Some(Kind::Alpha) => self.alpha_token(streams, column)?,
				Some(Kind::Digits0) => self.digits0_token(streams, column)?,
				Some(Kind::Delta0) => self.delta0_token(streams, column, previous)?,
				Some(Kind::Digits) => self.digits_token(streams, column)?,
				Some(Kind::Delta) => self.delta_token(streams, column, previous)?,
				Some(Kind::Nop) => EMPTY_TOKEN,
				Some(Kind::Match) => self.match_token(column, previous)?,
				_ => {
					self.reserve(1)?;
					self.output.push(0);
					self.current.push(Token {
						kind: TokenKind::End,
						number: 0,
						text: 0,
					});
					self.finish(text_start);
					return Ok(());
				}
			};
			self.current.push(token);
		}
		Err(corrupt())
	}

	fn finish(&mut self, text_start: usize) {
		let token_start = self.tokens.len();
		let token_count = self.current.len() - 1;
		self.tokens.append(&mut self.current);
		self.names.push(Name {
			text_start,
			text_len: self.output.len() - 1 - text_start,
			token_start,
			token_count,
		});
	}

	fn char_token(&mut self, streams: &mut Streams, column: usize) -> Result<Token, Error> {
		let byte = streams
			.get(column, Kind::Char)
			.and_then(Stream::next_byte)
			.ok_or_else(corrupt)?;
		self.reserve(1)?;
		self.output.push(byte);
		Ok(Token {
			kind: TokenKind::Char,
			number: u32::from(byte),
			text: 0,
		})
	}

	fn alpha_token(&mut self, streams: &mut Streams, column: usize) -> Result<Token, Error> {
		let text = streams
			.get(column, Kind::Alpha)
			.and_then(Stream::next_text)
			.ok_or_else(corrupt)?;
		if self.output.len() + text.len() + 1 > self.limit {
			return Err(corrupt());
		}
		let start = self.output.len();
		self.output.extend_from_slice(text);
		Ok(Token {
			kind: TokenKind::Alpha,
			number: u32::try_from(text.len()).map_err(|_| corrupt())?,
			text: start,
		})
	}

	fn digits0_token(&mut self, streams: &mut Streams, column: usize) -> Result<Token, Error> {
		let width = streams
			.get(column, Kind::DigitsZeroLen)
			.and_then(Stream::next_byte)
			.ok_or_else(corrupt)?;
		let value = streams
			.get(column, Kind::Digits0)
			.and_then(Stream::next_u32)
			.ok_or_else(corrupt)?;
		self.push_fixed(value, usize::from(width))?;
		Ok(Token {
			kind: TokenKind::Digits0,
			number: value,
			text: usize::from(width),
		})
	}

	fn delta0_token(
		&mut self,
		streams: &mut Streams,
		column: usize,
		previous: Option<Name>,
	) -> Result<Token, Error> {
		let base = self.previous_token(previous, column)?;
		let delta = streams
			.get(column, Kind::Delta0)
			.and_then(Stream::next_byte)
			.ok_or_else(corrupt)?;
		let value = u32::from(delta).wrapping_add(base.number);
		self.push_fixed(value, base.text)?;
		Ok(Token {
			kind: TokenKind::Digits0,
			number: value,
			text: base.text,
		})
	}

	fn digits_token(&mut self, streams: &mut Streams, column: usize) -> Result<Token, Error> {
		let value = streams
			.get(column, Kind::Digits)
			.and_then(Stream::next_u32)
			.ok_or_else(corrupt)?;
		self.push_variable(value)?;
		Ok(Token {
			kind: TokenKind::Digits,
			number: value,
			text: 0,
		})
	}

	fn delta_token(
		&mut self,
		streams: &mut Streams,
		column: usize,
		previous: Option<Name>,
	) -> Result<Token, Error> {
		let base = self.previous_token(previous, column)?;
		let delta = streams
			.get(column, Kind::Delta)
			.and_then(Stream::next_byte)
			.ok_or_else(corrupt)?;
		let value = u32::from(delta).wrapping_add(base.number);
		self.push_variable(value)?;
		Ok(Token {
			kind: TokenKind::Digits,
			number: value,
			text: 0,
		})
	}

	fn match_token(&mut self, column: usize, previous: Option<Name>) -> Result<Token, Error> {
		let base = self.previous_token(previous, column)?;
		match base.kind {
			TokenKind::Char => {
				self.reserve(1)?;
				self.output
					.push(u8::try_from(base.number).map_err(|_| corrupt())?);
			}
			TokenKind::Alpha => {
				let length = base.number as usize;
				self.reserve(length + 1)?;
				self.output
					.extend_from_within(base.text..base.text + length);
			}
			TokenKind::Digits => self.push_variable(base.number)?,
			TokenKind::Digits0 => self.push_fixed(base.number, base.text)?,
			TokenKind::Nop | TokenKind::End => return Err(corrupt()),
		}
		let text = if base.kind == TokenKind::Alpha {
			self.output.len() - base.number as usize
		} else {
			base.text
		};
		Ok(Token { text, ..base })
	}

	fn push_fixed(&mut self, value: u32, width: usize) -> Result<(), Error> {
		if width > MAX_FIXED_DIGITS
			|| u64::from(value) >= 10u64.pow(u32::try_from(width).map_err(|_| corrupt())?)
		{
			return Err(corrupt());
		}
		self.reserve(width + 20)?;
		let mut remaining = value;
		for exponent in (0..width).rev() {
			let place = 10u32.pow(u32::try_from(exponent).map_err(|_| corrupt())?);
			let digit = remaining / place;
			remaining %= place;
			self.output.push(b'0'.wrapping_add(digit.to_le_bytes()[0]));
		}
		Ok(())
	}

	fn push_variable(&mut self, value: u32) -> Result<(), Error> {
		self.reserve(20)?;
		if value != 0 {
			self.output.extend_from_slice(value.to_string().as_bytes());
		}
		Ok(())
	}
}
