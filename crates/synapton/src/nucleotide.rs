pub(crate) const CODE_N: u8 = 4;
pub(crate) const CODE_PAD: u8 = 5;

pub(crate) fn encode(ascii: &[u8], codes: &mut Vec<u8>) {
	codes.extend(ascii.iter().map(|&byte| match byte {
		b'A' | b'a' => 0,
		b'C' | b'c' => 1,
		b'G' | b'g' => 2,
		b'T' | b't' | b'U' | b'u' => 3,
		_ => CODE_N,
	}));
}

pub(crate) fn ascii(code: u8) -> u8 {
	match code {
		0 => b'A',
		1 => b'C',
		2 => b'G',
		3 => b'T',
		_ => b'N',
	}
}
