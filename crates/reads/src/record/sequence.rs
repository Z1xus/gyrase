use gyrase_genome::{Base, Code4};

const fn pair_table() -> [[u8; 2]; 256] {
	const ALPHABET: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
	let mut table = [[0; 2]; 256];
	let mut byte = 0;
	while byte < 256 {
		table[byte] = [ALPHABET[byte >> 4], ALPHABET[byte & 0xf]];
		byte += 1;
	}
	table
}

static PAIRS: [[u8; 2]; 256] = pair_table();

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sequence<'a> {
	packed: &'a [u8],
	len: usize,
}

impl<'a> Sequence<'a> {
	pub(crate) const fn new(packed: &'a [u8], len: usize) -> Self {
		Self { packed, len }
	}

	pub const fn len(&self) -> usize {
		self.len
	}

	pub const fn is_empty(&self) -> bool {
		self.len == 0
	}

	pub fn code(&self, index: usize) -> Option<Code4> {
		if index >= self.len {
			return None;
		}
		let byte = *self.packed.get(index / 2)?;
		let shift = if index.is_multiple_of(2) { 4 } else { 0 };
		Some(Code4::new(byte >> shift))
	}

	pub fn base(&self, index: usize) -> Option<Base> {
		self.code(index).map(Code4::base)
	}

	pub fn codes(&self) -> impl Iterator<Item = Code4> + 'a {
		let sequence = *self;
		(0..sequence.len).filter_map(move |index| sequence.code(index))
	}

	pub fn decode_into(&self, out: &mut Vec<u8>) {
		out.reserve(self.len);
		let whole = self.len / 2;
		for &byte in self.packed.iter().take(whole) {
			out.extend_from_slice(&PAIRS[usize::from(byte)]);
		}
		if self.len % 2 == 1
			&& let Some(&byte) = self.packed.get(whole)
		{
			out.push(PAIRS[usize::from(byte)][0]);
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Qualities<'a> {
	phred: &'a [u8],
}

impl<'a> Qualities<'a> {
	pub(crate) const fn new(phred: &'a [u8]) -> Self {
		Self { phred }
	}

	pub const fn len(&self) -> usize {
		self.phred.len()
	}

	pub const fn is_empty(&self) -> bool {
		self.phred.is_empty()
	}

	pub fn is_missing(&self) -> bool {
		self.phred.first() == Some(&0xff)
	}

	pub const fn as_phred(&self) -> &'a [u8] {
		self.phred
	}

	pub fn decode_ascii_into(&self, out: &mut Vec<u8>) {
		out.extend(self.phred.iter().map(|value| value.wrapping_add(33)));
	}
}
