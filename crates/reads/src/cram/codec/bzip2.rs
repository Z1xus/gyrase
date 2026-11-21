use std::io::Read;

use crate::cram::{Codec, Error};

pub(super) fn decode(input: &[u8], len: usize) -> Result<Vec<u8>, Error> {
	let mut output = Vec::with_capacity(len);
	bzip2::read::BzDecoder::new(input)
		.take(len as u64 + 1)
		.read_to_end(&mut output)
		.map_err(|_| Error::Corrupt(Codec::Bzip2))?;
	if output.len() == len {
		Ok(output)
	} else {
		Err(Error::Corrupt(Codec::Bzip2))
	}
}
