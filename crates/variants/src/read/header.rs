use std::io::{BufRead, Read};

use crate::Error;

const HEADER_END: &[u8] = b"#CHROM";
const MAX_LINE: u64 = 1 << 26;

pub fn read_header_text(source: &mut impl BufRead) -> Result<Vec<u8>, Error> {
	let mut text = Vec::new();
	loop {
		let line_start = text.len();
		match source.fill_buf()?.first() {
			None => return Err(Error::Header("missing #CHROM line")),
			Some(b'#') => {}
			Some(_) => return Err(Error::Header("record before the #CHROM line")),
		}
		Read::by_ref(source)
			.take(MAX_LINE)
			.read_until(b'\n', &mut text)?;
		let line = text.get(line_start..).unwrap_or_default();
		if !line.ends_with(b"\n")
			&& u64::try_from(line.len()).is_ok_and(|length| length >= MAX_LINE)
		{
			return Err(Error::Header("header line is too long"));
		}
		if line.starts_with(HEADER_END) {
			return Ok(text);
		}
	}
}
