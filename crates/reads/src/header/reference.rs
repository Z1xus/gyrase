use std::fmt::Write;

use super::fields::Fields;
use crate::error::HeaderError;

pub(crate) const MAX_REFERENCE_LEN: u32 = i32::MAX as u32;
const NAME: [u8; 2] = *b"SN";
const LENGTH: [u8; 2] = *b"LN";
const MD5: [u8; 2] = *b"M5";
const ALTERNATE_NAMES: [u8; 2] = *b"AN";

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
	name: String,
	length: u32,
	fields: Fields,
}

impl Reference {
	pub fn new(name: impl Into<String>, length: u32) -> Result<Self, HeaderError> {
		let name = name.into();
		if name.is_empty() || name.contains(['\t', '\n']) || length > MAX_REFERENCE_LEN {
			return Err(HeaderError::InvalidReferenceName);
		}
		let fields = Fields::new()
			.with(NAME, name.clone())
			.with(LENGTH, length.to_string());
		Ok(Self {
			name,
			length,
			fields,
		})
	}

	pub fn with_field(mut self, tag: [u8; 2], value: impl Into<String>) -> Self {
		self.fields = self.fields.with(tag, value);
		self
	}

	pub(crate) fn from_fields(fields: Fields, line: usize) -> Result<Self, HeaderError> {
		let name = fields
			.get(NAME)
			.ok_or(HeaderError::MissingField { line, tag: "SN" })?
			.to_owned();
		let length = fields
			.get(LENGTH)
			.ok_or(HeaderError::MissingField { line, tag: "LN" })?
			.parse::<u32>()
			.ok()
			.filter(|&length| length <= MAX_REFERENCE_LEN)
			.ok_or(HeaderError::InvalidLength { line })?;
		if name.is_empty() {
			return Err(HeaderError::InvalidReferenceName);
		}
		Ok(Self {
			name,
			length,
			fields,
		})
	}

	pub fn name(&self) -> &str {
		&self.name
	}

	pub const fn length(&self) -> u32 {
		self.length
	}

	pub fn fields(&self) -> &Fields {
		&self.fields
	}

	pub fn md5(&self) -> Option<[u8; 16]> {
		parse_md5(self.fields.get(MD5)?)
	}

	pub fn alternate_names(&self) -> impl Iterator<Item = &str> {
		self.fields
			.get(ALTERNATE_NAMES)
			.into_iter()
			.flat_map(|names| names.split(','))
			.filter(|name| !name.is_empty())
	}
}

fn parse_md5(text: &str) -> Option<[u8; 16]> {
	let digits = text.as_bytes();
	if digits.len() != 32 {
		return None;
	}
	let mut checksum = [0u8; 16];
	let (pairs, _) = digits.as_chunks::<2>();
	for (slot, &[high, low]) in checksum.iter_mut().zip(pairs) {
		*slot = hex_value(high)? << 4 | hex_value(low)?;
	}
	Some(checksum)
}

fn hex_value(digit: u8) -> Option<u8> {
	char::from(digit)
		.to_digit(16)
		.and_then(|value| u8::try_from(value).ok())
}

pub(crate) fn md5_hex(checksum: &[u8; 16]) -> String {
	let mut out = String::with_capacity(32);
	for byte in checksum {
		let _ = write!(out, "{byte:02x}");
	}
	out
}
