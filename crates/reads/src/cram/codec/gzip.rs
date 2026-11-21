use libdeflater::{CompressionLvl, Compressor, Decompressor};

use crate::cram::{Codec, Error};

const MAX_EXPANSION: usize = 1032;

pub(super) fn decode(input: &[u8], len: usize) -> Result<Vec<u8>, Error> {
	if len / MAX_EXPANSION > input.len() {
		return Err(Error::Corrupt(Codec::Gzip));
	}
	let mut output = vec![0; len];
	let written = Decompressor::new()
		.gzip_decompress(input, &mut output)
		.map_err(|_| Error::Corrupt(Codec::Gzip))?;
	if written == len {
		Ok(output)
	} else {
		Err(Error::Corrupt(Codec::Gzip))
	}
}

pub(super) fn encode(input: &[u8]) -> Option<Vec<u8>> {
	let mut compressor = Compressor::new(CompressionLvl::default());
	let mut output = vec![0; compressor.gzip_compress_bound(input.len())];
	let written = compressor.gzip_compress(input, &mut output).ok()?;
	output.truncate(written);
	Some(output)
}
