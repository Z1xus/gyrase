use std::io::BufRead;

use crate::Error;

const HEADER_END: &[u8] = b"#CHROM";

pub fn read_header_text(source: &mut impl BufRead) -> Result<Vec<u8>, Error> {
	let mut text = Vec::new();
	loop {
		let line_start = text.len();
		if source.read_until(b'\n', &mut text)? == 0 {
			return Err(Error::Header("missing #CHROM line"));
		}
		let line = text.get(line_start..).unwrap_or_default();
		if !line.starts_with(b"#") {
			return Err(Error::Header("record before the #CHROM line"));
		}
		if line.starts_with(HEADER_END) {
			return Ok(text);
		}
	}
}
