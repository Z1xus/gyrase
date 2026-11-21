use std::io::BufRead;

use memchr::memchr;

use super::batch::{FastqBatch, Span};
use crate::error::{Error, FastqError};

const MAX_BATCH_BYTES: usize = 1 << 30;

#[derive(Debug)]
pub struct FastqReader<R: BufRead> {
	source: R,
	records_read: u64,
}

struct Lines {
	header_end: usize,
	bases_end: usize,
	separator_end: usize,
	qualities_end: usize,
	next: usize,
}

impl<R: BufRead> FastqReader<R> {
	pub fn new(source: R) -> Self {
		Self {
			source,
			records_read: 0,
		}
	}

	pub fn read_batch(
		&mut self,
		batch: &mut FastqBatch,
		max_records: usize,
	) -> Result<usize, Error> {
		batch.restart();
		let mut source_done = false;
		while batch.len() < max_records {
			let record = self.records_read + 1;
			let parsed = parse_record(&batch.data, batch.consumed, source_done, record)?;
			match parsed {
				Some((span, next)) => {
					self.records_read += 1;
					batch.spans.push(span);
					batch.consumed = next;
				}
				None if source_done => {
					if has_content(batch) {
						return Err(FastqError::Truncated {
							record: self.records_read + 1,
						}
						.into());
					}
					break;
				}
				None => source_done = self.refill(batch)?,
			}
		}
		Ok(batch.len())
	}

	fn refill(&mut self, batch: &mut FastqBatch) -> Result<bool, Error> {
		let chunk = self.source.fill_buf()?;
		if chunk.is_empty() {
			return Ok(true);
		}
		if batch.data.len() + chunk.len() > MAX_BATCH_BYTES {
			return Err(FastqError::RecordTooLarge {
				record: self.records_read + 1,
			}
			.into());
		}
		batch.data.extend_from_slice(chunk);
		let taken = chunk.len();
		self.source.consume(taken);
		Ok(false)
	}
}

fn has_content(batch: &FastqBatch) -> bool {
	batch
		.data
		.get(batch.consumed..)
		.is_some_and(|rest| rest.iter().any(|byte| !byte.is_ascii_whitespace()))
}

fn parse_record(
	data: &[u8],
	from: usize,
	at_end: bool,
	record: u64,
) -> Result<Option<(Span, usize)>, FastqError> {
	let Some(lines) = locate_lines(data, from, at_end) else {
		return Ok(None);
	};
	let line = |start: usize, end: usize| {
		let text = data.get(start..end).unwrap_or_default();
		start + text.strip_suffix(b"\r").unwrap_or(text).len()
	};
	let header_end = line(from, lines.header_end);
	let bases_start = lines.header_end + 1;
	let bases_end = line(bases_start, lines.bases_end);
	let separator_start = lines.bases_end + 1;
	let qualities_start = lines.separator_end + 1;
	let qualities_end = line(qualities_start, lines.qualities_end);
	if data.get(from) != Some(&b'@') {
		return Err(FastqError::MissingHeader { record });
	}
	if data.get(separator_start) != Some(&b'+') {
		return Err(FastqError::MissingSeparator { record });
	}
	if bases_end - bases_start != qualities_end - qualities_start {
		return Err(FastqError::LengthMismatch { record });
	}
	let name_start = from + 1;
	let name_len = data
		.get(name_start..header_end)
		.unwrap_or_default()
		.iter()
		.position(|byte| matches!(byte, b' ' | b'\t'))
		.unwrap_or(header_end - name_start);
	let too_large = FastqError::RecordTooLarge { record };
	let offset = |value: usize| u32::try_from(value).map_err(|_| too_large.clone());
	let span = Span {
		name_start: offset(name_start)?,
		name_end: offset(name_start + name_len)?,
		header_end: offset(header_end)?,
		bases_start: offset(bases_start)?,
		bases_end: offset(bases_end)?,
		qualities_start: offset(qualities_start)?,
		qualities_end: offset(qualities_end)?,
	};
	Ok(Some((span, lines.next)))
}

fn locate_lines(data: &[u8], from: usize, at_end: bool) -> Option<Lines> {
	let end_of = |start: usize| -> Option<usize> {
		let rest = data.get(start..)?;
		memchr(b'\n', rest).map(|offset| start + offset)
	};
	let header_end = end_of(from)?;
	let bases_end = end_of(header_end + 1)?;
	let separator_end = end_of(bases_end + 1)?;
	let qualities_start = separator_end + 1;
	let (qualities_end, next) = match end_of(qualities_start) {
		Some(end) => (end, end + 1),
		None if at_end && qualities_start <= data.len() => (data.len(), data.len()),
		None => return None,
	};
	Some(Lines {
		header_end,
		bases_end,
		separator_end,
		qualities_end,
		next,
	})
}
