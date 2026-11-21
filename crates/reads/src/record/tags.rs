use std::marker::PhantomData;

use crate::error::RecordError;

pub trait Element: Copy + 'static {
	const CODE: u8;
	const WIDTH: usize;

	fn read(bytes: &[u8]) -> Option<Self>;
	fn write(self, out: &mut Vec<u8>);
}

macro_rules! element {
	($type:ty, $code:expr) => {
		impl Element for $type {
			const CODE: u8 = $code;
			const WIDTH: usize = size_of::<$type>();

			fn read(bytes: &[u8]) -> Option<Self> {
				Some(Self::from_le_bytes(*bytes.first_chunk()?))
			}

			fn write(self, out: &mut Vec<u8>) {
				out.extend_from_slice(&self.to_le_bytes());
			}
		}
	};
}

element!(i8, b'c');
element!(u8, b'C');
element!(i16, b's');
element!(u16, b'S');
element!(i32, b'i');
element!(u32, b'I');
element!(f32, b'f');

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Elements<'a, T: Element> {
	bytes: &'a [u8],
	element: PhantomData<T>,
}

impl<'a, T: Element> Elements<'a, T> {
	fn new(bytes: &'a [u8]) -> Self {
		Self {
			bytes,
			element: PhantomData,
		}
	}

	pub fn len(&self) -> usize {
		self.bytes.len() / T::WIDTH
	}

	pub fn is_empty(&self) -> bool {
		self.bytes.is_empty()
	}

	pub fn get(&self, index: usize) -> Option<T> {
		T::read(self.bytes.get(index.checked_mul(T::WIDTH)?..)?)
	}

	pub fn iter(&self) -> impl Iterator<Item = T> + 'a {
		self.bytes.chunks_exact(T::WIDTH).filter_map(T::read)
	}

	pub(crate) fn as_bytes(&self) -> &'a [u8] {
		self.bytes
	}
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Array<'a> {
	Int8(Elements<'a, i8>),
	UInt8(Elements<'a, u8>),
	Int16(Elements<'a, i16>),
	UInt16(Elements<'a, u16>),
	Int32(Elements<'a, i32>),
	UInt32(Elements<'a, u32>),
	Float(Elements<'a, f32>),
}

impl<'a> Array<'a> {
	fn from_parts(code: u8, bytes: &'a [u8]) -> Option<Self> {
		Some(match code {
			b'c' => Self::Int8(Elements::new(bytes)),
			b'C' => Self::UInt8(Elements::new(bytes)),
			b's' => Self::Int16(Elements::new(bytes)),
			b'S' => Self::UInt16(Elements::new(bytes)),
			b'i' => Self::Int32(Elements::new(bytes)),
			b'I' => Self::UInt32(Elements::new(bytes)),
			b'f' => Self::Float(Elements::new(bytes)),
			_ => return None,
		})
	}

	pub fn len(&self) -> usize {
		match self {
			Self::Int8(items) => items.len(),
			Self::UInt8(items) => items.len(),
			Self::Int16(items) => items.len(),
			Self::UInt16(items) => items.len(),
			Self::Int32(items) => items.len(),
			Self::UInt32(items) => items.len(),
			Self::Float(items) => items.len(),
		}
	}

	pub fn is_empty(&self) -> bool {
		self.len() == 0
	}

	pub(crate) fn code_and_bytes(&self) -> (u8, &'a [u8]) {
		match self {
			Self::Int8(items) => (i8::CODE, items.as_bytes()),
			Self::UInt8(items) => (u8::CODE, items.as_bytes()),
			Self::Int16(items) => (i16::CODE, items.as_bytes()),
			Self::UInt16(items) => (u16::CODE, items.as_bytes()),
			Self::Int32(items) => (i32::CODE, items.as_bytes()),
			Self::UInt32(items) => (u32::CODE, items.as_bytes()),
			Self::Float(items) => (f32::CODE, items.as_bytes()),
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Value<'a> {
	Char(u8),
	Int8(i8),
	UInt8(u8),
	Int16(i16),
	UInt16(u16),
	Int32(i32),
	UInt32(u32),
	Float(f32),
	Double(f64),
	String(&'a [u8]),
	Hex(&'a [u8]),
	Array(Array<'a>),
}

impl Value<'_> {
	pub fn as_integer(&self) -> Option<i64> {
		Some(match *self {
			Self::Int8(value) => i64::from(value),
			Self::UInt8(value) => i64::from(value),
			Self::Int16(value) => i64::from(value),
			Self::UInt16(value) => i64::from(value),
			Self::Int32(value) => i64::from(value),
			Self::UInt32(value) => i64::from(value),
			_ => return None,
		})
	}

	pub fn as_string(&self) -> Option<&[u8]> {
		match *self {
			Self::String(text) => Some(text),
			_ => None,
		}
	}

	pub(crate) fn encode(&self, tag: [u8; 2], out: &mut Vec<u8>) {
		out.extend_from_slice(&tag);
		match *self {
			Self::Char(value) => out.extend_from_slice(&[b'A', value]),
			Self::Int8(value) => encode_fixed(out, value),
			Self::UInt8(value) => encode_fixed(out, value),
			Self::Int16(value) => encode_fixed(out, value),
			Self::UInt16(value) => encode_fixed(out, value),
			Self::Int32(value) => encode_fixed(out, value),
			Self::UInt32(value) => encode_fixed(out, value),
			Self::Float(value) => encode_fixed(out, value),
			Self::Double(value) => {
				out.push(b'd');
				out.extend_from_slice(&value.to_le_bytes());
			}
			Self::String(text) => encode_text(out, b'Z', text),
			Self::Hex(text) => encode_text(out, b'H', text),
			Self::Array(array) => {
				let (code, bytes) = array.code_and_bytes();
				out.extend_from_slice(&[b'B', code]);
				out.extend_from_slice(
					&u32::try_from(array.len()).unwrap_or(u32::MAX).to_le_bytes(),
				);
				out.extend_from_slice(bytes);
			}
		}
	}
}

fn encode_fixed<T: Element>(out: &mut Vec<u8>, value: T) {
	out.push(T::CODE);
	value.write(out);
}

fn encode_text(out: &mut Vec<u8>, kind: u8, text: &[u8]) {
	out.push(kind);
	out.extend_from_slice(text);
	out.push(0);
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Field<'a> {
	pub tag: [u8; 2],
	pub value: Value<'a>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Tags<'a> {
	bytes: &'a [u8],
}

impl<'a> Tags<'a> {
	pub(crate) const fn new(bytes: &'a [u8]) -> Self {
		Self { bytes }
	}

	pub fn is_empty(&self) -> bool {
		self.bytes.is_empty()
	}

	pub fn iter(&self) -> TagIter<'a> {
		TagIter { rest: self.bytes }
	}

	pub fn get(&self, tag: [u8; 2]) -> Result<Option<Value<'a>>, RecordError> {
		for field in self {
			let field = field?;
			if field.tag == tag {
				return Ok(Some(field.value));
			}
		}
		Ok(None)
	}
}

impl<'a> IntoIterator for &Tags<'a> {
	type Item = Result<Field<'a>, RecordError>;
	type IntoIter = TagIter<'a>;

	fn into_iter(self) -> TagIter<'a> {
		self.iter()
	}
}

#[derive(Debug, Clone)]
pub struct TagIter<'a> {
	rest: &'a [u8],
}

impl<'a> Iterator for TagIter<'a> {
	type Item = Result<Field<'a>, RecordError>;

	fn next(&mut self) -> Option<Self::Item> {
		if self.rest.is_empty() {
			return None;
		}
		match parse_field(self.rest) {
			Ok((field, rest)) => {
				self.rest = rest;
				Some(Ok(field))
			}
			Err(error) => {
				self.rest = &[];
				Some(Err(error))
			}
		}
	}
}

fn fixed<const N: usize>(bytes: &[u8]) -> Result<([u8; N], &[u8]), RecordError> {
	let (head, rest) = bytes
		.split_first_chunk::<N>()
		.ok_or(RecordError::MalformedTag)?;
	Ok((*head, rest))
}

fn parse_field(bytes: &[u8]) -> Result<(Field<'_>, &[u8]), RecordError> {
	let ([first, second, kind], rest) = fixed::<3>(bytes)?;
	let tag = [first, second];
	let (value, rest) = match kind {
		b'A' => {
			let ([value], rest) = fixed::<1>(rest)?;
			(Value::Char(value), rest)
		}
		b'c' => {
			let (raw, rest) = fixed::<1>(rest)?;
			(Value::Int8(i8::from_le_bytes(raw)), rest)
		}
		b'C' => {
			let ([value], rest) = fixed::<1>(rest)?;
			(Value::UInt8(value), rest)
		}
		b's' => {
			let (raw, rest) = fixed::<2>(rest)?;
			(Value::Int16(i16::from_le_bytes(raw)), rest)
		}
		b'S' => {
			let (raw, rest) = fixed::<2>(rest)?;
			(Value::UInt16(u16::from_le_bytes(raw)), rest)
		}
		b'i' => {
			let (raw, rest) = fixed::<4>(rest)?;
			(Value::Int32(i32::from_le_bytes(raw)), rest)
		}
		b'I' => {
			let (raw, rest) = fixed::<4>(rest)?;
			(Value::UInt32(u32::from_le_bytes(raw)), rest)
		}
		b'f' => {
			let (raw, rest) = fixed::<4>(rest)?;
			(Value::Float(f32::from_le_bytes(raw)), rest)
		}
		b'd' => {
			let (raw, rest) = fixed::<8>(rest)?;
			(Value::Double(f64::from_le_bytes(raw)), rest)
		}
		b'Z' | b'H' => {
			let end = memchr::memchr(0, rest).ok_or(RecordError::MalformedTag)?;
			let (text, rest) = rest
				.split_at_checked(end)
				.ok_or(RecordError::MalformedTag)?;
			let rest = rest.get(1..).ok_or(RecordError::MalformedTag)?;
			let value = if kind == b'Z' {
				Value::String(text)
			} else {
				Value::Hex(text)
			};
			(value, rest)
		}
		b'B' => parse_array(rest)?,
		_ => return Err(RecordError::MalformedTag),
	};
	Ok((Field { tag, value }, rest))
}

fn parse_array(bytes: &[u8]) -> Result<(Value<'_>, &[u8]), RecordError> {
	let ([code], rest) = fixed::<1>(bytes)?;
	let (count, rest) = fixed::<4>(rest)?;
	let width = match code {
		b'c' | b'C' => 1,
		b's' | b'S' => 2,
		b'i' | b'I' | b'f' => 4,
		_ => return Err(RecordError::MalformedTag),
	};
	let count =
		usize::try_from(u32::from_le_bytes(count)).map_err(|_| RecordError::MalformedTag)?;
	let len = count.checked_mul(width).ok_or(RecordError::MalformedTag)?;
	let (items, rest) = rest
		.split_at_checked(len)
		.ok_or(RecordError::MalformedTag)?;
	let array = Array::from_parts(code, items).ok_or(RecordError::MalformedTag)?;
	Ok((Value::Array(array), rest))
}
