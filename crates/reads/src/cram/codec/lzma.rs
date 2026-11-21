use std::io::Read;

use lzma_rust2::XzReader;

use crate::cram::{Codec, Error};

pub(super) fn decode(input: &[u8], len: usize) -> Result<Vec<u8>, Error> {
	let mut output = super::first_reserve(len);
	XzReader::new(input, false)
		.take(len as u64 + 1)
		.read_to_end(&mut output)
		.map_err(|_| Error::Corrupt(Codec::Lzma))?;
	if output.len() == len {
		Ok(output)
	} else {
		Err(Error::Corrupt(Codec::Lzma))
	}
}
