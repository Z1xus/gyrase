#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Base {
	A,
	C,
	G,
	T,
	N,
}

impl Base {
	pub const ACGT: [Self; 4] = [Self::A, Self::C, Self::G, Self::T];

	pub const fn from_ascii(byte: u8) -> Self {
		match byte {
			b'A' | b'a' => Self::A,
			b'C' | b'c' => Self::C,
			b'G' | b'g' => Self::G,
			b'T' | b't' | b'U' | b'u' => Self::T,
			_ => Self::N,
		}
	}

	pub const fn ascii(self) -> u8 {
		match self {
			Self::A => b'A',
			Self::C => b'C',
			Self::G => b'G',
			Self::T => b'T',
			Self::N => b'N',
		}
	}

	pub const fn index(self) -> Option<usize> {
		match self {
			Self::A => Some(0),
			Self::C => Some(1),
			Self::G => Some(2),
			Self::T => Some(3),
			Self::N => None,
		}
	}

	pub const fn complement(self) -> Self {
		match self {
			Self::A => Self::T,
			Self::C => Self::G,
			Self::G => Self::C,
			Self::T => Self::A,
			Self::N => Self::N,
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Code4(u8);

impl Code4 {
	const ALPHABET: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

	pub const fn new(code: u8) -> Self {
		Self(code & 0xf)
	}

	pub const fn from_ascii(byte: u8) -> Self {
		Self(match byte.to_ascii_uppercase() {
			b'=' => 0,
			b'A' => 1,
			b'C' => 2,
			b'M' => 3,
			b'G' => 4,
			b'R' => 5,
			b'S' => 6,
			b'V' => 7,
			b'T' | b'U' => 8,
			b'W' => 9,
			b'Y' => 10,
			b'H' => 11,
			b'K' => 12,
			b'D' => 13,
			b'B' => 14,
			_ => 15,
		})
	}

	pub const fn code(self) -> u8 {
		self.0
	}

	pub const fn ascii(self) -> u8 {
		Self::ALPHABET[self.0 as usize]
	}

	pub const fn base(self) -> Base {
		match self.0 {
			1 => Base::A,
			2 => Base::C,
			4 => Base::G,
			8 => Base::T,
			_ => Base::N,
		}
	}
}
