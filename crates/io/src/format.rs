use std::fmt;
use std::io::{self, BufRead, BufReader, Read};

use flate2::bufread::MultiGzDecoder;

use crate::bgzf::is_bgzf;
use crate::{BgzfReader, Error, Input, Workers};

const PEEK_LEN: usize = 18;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
	None,
	Gzip,
	Bgzf,
}

type Source = BufReader<Box<dyn Read + Send>>;

pub enum Decoded {
	Plain(Source),
	Gzip(Box<BufReader<MultiGzDecoder<Source>>>),
	Bgzf(BgzfReader),
}

impl fmt::Debug for Decoded {
	fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
		formatter.write_str(match self.compression() {
			Compression::None => "Plain",
			Compression::Gzip => "Gzip",
			Compression::Bgzf => "Bgzf",
		})
	}
}

impl Decoded {
	pub fn compression(&self) -> Compression {
		match self {
			Self::Plain(_) => Compression::None,
			Self::Gzip(_) => Compression::Gzip,
			Self::Bgzf(_) => Compression::Bgzf,
		}
	}
}

fn compression_of(head: &[u8]) -> Compression {
	match head {
		[0x1f, 0x8b, ..] if head.len() >= PEEK_LEN && is_bgzf(head) => Compression::Bgzf,
		[0x1f, 0x8b, ..] => Compression::Gzip,
		_ => Compression::None,
	}
}

pub fn peek(source: impl Read, limit: u64) -> Result<(Compression, Vec<u8>), Error> {
	let mut buffered = BufReader::new(source);
	let compression = compression_of(buffered.fill_buf()?);
	let mut head = Vec::new();
	match compression {
		Compression::None => buffered.take(limit).read_to_end(&mut head)?,
		Compression::Gzip | Compression::Bgzf => MultiGzDecoder::new(buffered)
			.take(limit)
			.read_to_end(&mut head)?,
	};
	Ok((compression, head))
}

pub fn decode(source: Box<dyn Read + Send>, workers: &Workers) -> Result<Decoded, Error> {
	let mut buffered = BufReader::with_capacity(1 << 20, source);
	let compression = compression_of(buffered.fill_buf()?);
	Ok(match compression {
		Compression::None => Decoded::Plain(buffered),
		Compression::Gzip => Decoded::Gzip(Box::new(BufReader::with_capacity(
			1 << 20,
			MultiGzDecoder::new(buffered),
		))),
		Compression::Bgzf => {
			Decoded::Bgzf(BgzfReader::new(Input::Stream(Box::new(buffered)), workers))
		}
	})
}

impl Read for Decoded {
	fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
		match self {
			Self::Plain(reader) => reader.read(buffer),
			Self::Gzip(reader) => reader.read(buffer),
			Self::Bgzf(reader) => reader.read(buffer),
		}
	}
}

impl BufRead for Decoded {
	fn fill_buf(&mut self) -> io::Result<&[u8]> {
		match self {
			Self::Plain(reader) => reader.fill_buf(),
			Self::Gzip(reader) => reader.fill_buf(),
			Self::Bgzf(reader) => reader.fill_buf(),
		}
	}

	fn consume(&mut self, amount: usize) {
		match self {
			Self::Plain(reader) => reader.consume(amount),
			Self::Gzip(reader) => reader.consume(amount),
			Self::Bgzf(reader) => reader.consume(amount),
		}
	}
}
