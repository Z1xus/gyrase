use super::encoded::{Encoded, Format, Span};
use super::sink::Sink;
use crate::Error;
use crate::header::Header;

const BCF_MAGIC: [u8; 5] = *b"BCF\x02\x02";

#[derive(Debug, Clone, Copy)]
pub struct Placed<M> {
	pub begin: M,
	pub end: M,
}

#[derive(Debug)]
pub struct Writer<S: Sink> {
	sink: S,
}

impl<S: Sink> Writer<S> {
	pub fn new(mut sink: S, header: &Header, format: Format) -> Result<Self, Error> {
		let mut text = Vec::new();
		match format {
			Format::Vcf => header.write_vcf(&mut text),
			Format::Bcf => {
				header.write_bcf(&mut text);
				text.push(0);
				let length =
					u32::try_from(text.len()).map_err(|_| Error::Limits("header length"))?;
				sink.write_all(&BCF_MAGIC)?;
				sink.write_all(&length.to_le_bytes())?;
			}
		}
		sink.write_all(&text)?;
		sink.end_block()?;
		Ok(Self { sink })
	}

	pub fn write(
		&mut self,
		encoded: &Encoded,
		mut placed: impl FnMut(Span, Placed<S::Mark>),
	) -> Result<(), Error> {
		for (record, span) in encoded.records() {
			let begin = self.sink.mark();
			self.sink.write_all(record)?;
			placed(
				span,
				Placed {
					begin,
					end: self.sink.mark(),
				},
			);
		}
		Ok(())
	}

	pub fn append(&mut self, encoded: &Encoded) -> Result<(), Error> {
		self.write(encoded, |_, _| {})
	}

	pub fn finish(self) -> Result<S::Done, Error> {
		self.sink.finish()
	}
}
