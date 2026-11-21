use crate::Error;
use crate::bcf::{BcfRecord, Encoder};
use crate::header::Header;
use crate::record::Record;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
	Vcf,
	Bcf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
	pub contig: u32,
	pub start: u32,
	pub end: u32,
}

#[derive(Debug, Default)]
pub struct Encoded {
	bytes: Vec<u8>,
	ends: Vec<usize>,
	spans: Vec<Span>,
}

impl Encoded {
	pub fn clear(&mut self) {
		self.bytes.clear();
		self.ends.clear();
		self.spans.clear();
	}

	pub fn len(&self) -> usize {
		self.ends.len()
	}

	pub fn is_empty(&self) -> bool {
		self.ends.is_empty()
	}

	pub fn records(&self) -> impl Iterator<Item = (&[u8], Span)> {
		let mut start = 0;
		self.ends.iter().zip(&self.spans).map(move |(&end, &span)| {
			let record = self.bytes.get(start..end).unwrap_or_default();
			start = end;
			(record, span)
		})
	}

	pub(crate) fn bytes_mut(&mut self) -> &mut Vec<u8> {
		&mut self.bytes
	}

	pub(crate) fn seal(&mut self, span: Span) {
		self.ends.push(self.bytes.len());
		self.spans.push(span);
	}
}

#[derive(Debug)]
pub struct Formatter {
	format: Format,
	encoder: Encoder,
}

impl Formatter {
	pub fn new(format: Format) -> Self {
		Self {
			format,
			encoder: Encoder::new(),
		}
	}

	pub fn format(&self) -> Format {
		self.format
	}

	pub fn push_record(
		&mut self,
		header: &Header,
		record: &Record,
		out: &mut Encoded,
	) -> Result<(), Error> {
		match self.format {
			Format::Vcf => {
				record.write_vcf(header, out.bytes_mut())?;
				out.bytes_mut().push(b'\n');
			}
			Format::Bcf => {
				let (shared, indiv) = self.encoder.encode(record, header)?;
				push_bcf_parts(out.bytes_mut(), shared, indiv)?;
			}
		}
		out.seal(Span {
			contig: record.contig,
			start: record.start,
			end: record.end(header),
		});
		Ok(())
	}

	pub fn push_line(&self, line: &[u8], span: Span, out: &mut Encoded) {
		out.bytes_mut().extend_from_slice(line);
		out.bytes_mut().push(b'\n');
		out.seal(span);
	}

	pub fn push_bcf_as_vcf(
		&self,
		header: &Header,
		source: &BcfRecord<'_>,
		span: Span,
		out: &mut Encoded,
	) -> Result<(), Error> {
		source.write_vcf(header, out.bytes_mut())?;
		out.bytes_mut().push(b'\n');
		out.seal(span);
		Ok(())
	}

	pub fn push_bcf(
		&self,
		shared: &[u8],
		indiv: &[u8],
		span: Span,
		out: &mut Encoded,
	) -> Result<(), Error> {
		push_bcf_parts(out.bytes_mut(), shared, indiv)?;
		out.seal(span);
		Ok(())
	}
}

fn push_bcf_parts(out: &mut Vec<u8>, shared: &[u8], indiv: &[u8]) -> Result<(), Error> {
	let shared_len = u32::try_from(shared.len()).map_err(|_| Error::Limits("record length"))?;
	let indiv_len = u32::try_from(indiv.len()).map_err(|_| Error::Limits("record length"))?;
	out.extend_from_slice(&shared_len.to_le_bytes());
	out.extend_from_slice(&indiv_len.to_le_bytes());
	out.extend_from_slice(shared);
	out.extend_from_slice(indiv);
	Ok(())
}
