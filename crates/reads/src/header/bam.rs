use super::{Header, Reference};
use crate::error::{BamError, Error, HeaderError};

const MAGIC: [u8; 4] = *b"BAM\x01";
const MAX_TEXT_LEN: usize = 1 << 28;
const MAX_REFERENCES: usize = 1 << 28;
const MAX_NAME_LEN: usize = 1 << 16;
const PREALLOCATED_REFERENCES: usize = 1 << 12;

pub(crate) trait ByteSource {
	fn take(&mut self, len: usize, out: &mut Vec<u8>) -> Result<(), Error>;
}

impl ByteSource for &[u8] {
	fn take(&mut self, len: usize, out: &mut Vec<u8>) -> Result<(), Error> {
		let (head, rest) = self.split_at_checked(len).ok_or(HeaderError::Truncated)?;
		out.clear();
		out.extend_from_slice(head);
		*self = rest;
		Ok(())
	}
}

fn take_i32(source: &mut impl ByteSource, scratch: &mut Vec<u8>) -> Result<i32, Error> {
	source.take(4, scratch)?;
	let raw = scratch.first_chunk::<4>().ok_or(HeaderError::Truncated)?;
	Ok(i32::from_le_bytes(*raw))
}

fn take_length(
	source: &mut impl ByteSource,
	scratch: &mut Vec<u8>,
	limit: usize,
) -> Result<usize, Error> {
	usize::try_from(take_i32(source, scratch)?)
		.ok()
		.filter(|&len| len <= limit)
		.ok_or_else(|| HeaderError::TooLarge { limit }.into())
}

pub(super) fn read(source: &mut impl ByteSource) -> Result<Header, Error> {
	let mut scratch = Vec::new();
	source.take(4, &mut scratch)?;
	if scratch != MAGIC {
		return Err(BamError::NotBam.into());
	}
	let text_len = take_length(source, &mut scratch, MAX_TEXT_LEN)?;
	source.take(text_len, &mut scratch)?;
	let text = scratch.split(|&byte| byte == 0).next().unwrap_or_default();
	let mut header = Header::parse(text)?;
	let count = take_length(source, &mut scratch, MAX_REFERENCES)?;
	let mut binary = Vec::with_capacity(count.min(PREALLOCATED_REFERENCES));
	for _ in 0..count {
		let name_len = take_length(source, &mut scratch, MAX_NAME_LEN)?;
		source.take(name_len, &mut scratch)?;
		let name = scratch.strip_suffix(&[0]).ok_or(HeaderError::Truncated)?;
		let name = String::from_utf8(name.to_vec()).map_err(|_| HeaderError::InvalidUtf8)?;
		let length = u32::try_from(take_i32(source, &mut scratch)?)
			.map_err(|_| HeaderError::InvalidReferenceName)?;
		binary.push(Reference::new(name, length)?);
	}
	reconcile(&mut header, binary)?;
	Ok(header)
}

fn reconcile(header: &mut Header, binary: Vec<Reference>) -> Result<(), HeaderError> {
	if header.references.is_empty() {
		return header.insert_references(binary);
	}
	let consistent = header.references.len() == binary.len();
	let mismatch = header
		.references
		.iter()
		.zip(&binary)
		.position(|(text, binary)| {
			text.name() != binary.name() || text.length() != binary.length()
		});
	match (consistent, mismatch) {
		(true, None) => Ok(()),
		(_, Some(index)) => Err(HeaderError::ReferenceMismatch { index }),
		(false, None) => Err(HeaderError::ReferenceMismatch {
			index: header.references.len().min(binary.len()),
		}),
	}
}

pub(super) fn write(header: &Header, out: &mut Vec<u8>) -> Result<(), HeaderError> {
	let text = header.text();
	let text_len = i32::try_from(text.len()).map_err(|_| HeaderError::TooLarge {
		limit: i32::MAX as usize,
	})?;
	let count = i32::try_from(header.references.len()).map_err(|_| HeaderError::TooLarge {
		limit: i32::MAX as usize,
	})?;
	out.extend_from_slice(&MAGIC);
	out.extend_from_slice(&text_len.to_le_bytes());
	out.extend_from_slice(text.as_bytes());
	out.extend_from_slice(&count.to_le_bytes());
	for reference in &header.references {
		let name_len = i32::try_from(reference.name().len() + 1)
			.map_err(|_| HeaderError::InvalidReferenceName)?;
		out.extend_from_slice(&name_len.to_le_bytes());
		out.extend_from_slice(reference.name().as_bytes());
		out.push(0);
		let length = i32::try_from(reference.length())
			.map_err(|_| HeaderError::ReferenceTooLong(reference.name().to_owned()))?;
		out.extend_from_slice(&length.to_le_bytes());
	}
	Ok(())
}
