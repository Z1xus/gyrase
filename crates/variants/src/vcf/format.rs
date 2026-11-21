use crate::Error;
use crate::header::{Header, Version};
use crate::number::{push_float, push_int};
use crate::record::{Field, Record};
use crate::value::{Cells, END_INT, MISSING_INT, Store, is_end_float, is_missing_float};

const PREFIXED_PHASING: Version = Version { major: 4, minor: 4 };
const GENOTYPE_KEY: &str = "GT";
const MISSING_STRING_BYTE: u8 = 7;

impl Record {
	pub fn write_vcf(&self, header: &Header, out: &mut Vec<u8>) -> Result<(), Error> {
		let contig = header
			.contig(self.contig)
			.ok_or(Error::Bcf("contig id is not in the header"))?;
		out.extend_from_slice(contig.name.as_bytes());
		out.push(b'\t');
		push_int(out, i64::from(self.start) + 1);
		out.push(b'\t');
		push_or_dot(out, &self.id);
		out.push(b'\t');
		push_or_dot(out, self.allele(0).unwrap_or_default());
		out.push(b'\t');
		self.write_alts(out);
		out.push(b'\t');
		match self.qual {
			Some(qual) if !is_missing_float(qual) => push_float(out, qual),
			_ => out.push(b'.'),
		}
		out.push(b'\t');
		write_filters(self, header, out)?;
		out.push(b'\t');
		write_info(self, header, out)?;
		if header.sample_count() > 0 {
			write_samples(self, header, out)?;
		}
		Ok(())
	}

	fn write_alts(&self, out: &mut Vec<u8>) {
		if self.allele_count() < 2 {
			out.push(b'.');
			return;
		}
		for (index, allele) in self.alleles().skip(1).enumerate() {
			if index > 0 {
				out.push(b',');
			}
			out.extend_from_slice(allele);
		}
	}
}

pub(crate) fn push_or_dot(out: &mut Vec<u8>, bytes: &[u8]) {
	if bytes.is_empty() {
		out.push(b'.');
	} else {
		out.extend_from_slice(bytes);
	}
}

pub(crate) fn tag_name(header: &Header, key: u32) -> Result<&str, Error> {
	header
		.name(key)
		.ok_or(Error::Bcf("tag id is not in the header"))
}

fn write_filters(record: &Record, header: &Header, out: &mut Vec<u8>) -> Result<(), Error> {
	if record.filters.is_empty() {
		out.push(b'.');
		return Ok(());
	}
	for (index, &key) in record.filters.iter().enumerate() {
		if index > 0 {
			out.push(b';');
		}
		out.extend_from_slice(tag_name(header, key)?.as_bytes());
	}
	Ok(())
}

pub(crate) fn write_info(record: &Record, header: &Header, out: &mut Vec<u8>) -> Result<(), Error> {
	if record.info.is_empty() {
		out.push(b'.');
		return Ok(());
	}
	for (index, field) in record.info.iter().enumerate() {
		if index > 0 {
			out.push(b';');
		}
		out.extend_from_slice(tag_name(header, field.key)?.as_bytes());
		if field.width == 0 || matches!(field.cells, Cells::Flag) {
			continue;
		}
		out.push(b'=');
		push_cells(out, &record.store, &field.cells);
	}
	Ok(())
}

fn push_cells(out: &mut Vec<u8>, store: &Store, cells: &Cells) {
	match cells {
		Cells::Flag => {}
		Cells::Ints(_) => push_int_list(out, store.ints(cells).iter().copied()),
		Cells::Floats(_) => push_float_list(out, store.floats(cells).iter().copied()),
		Cells::Text(_) => push_text(out, store.text(cells)),
	}
}

pub(crate) fn push_int_list(out: &mut Vec<u8>, values: impl Iterator<Item = i32>) {
	for (index, value) in values.enumerate() {
		if value == END_INT {
			break;
		}
		if index > 0 {
			out.push(b',');
		}
		if value == MISSING_INT {
			out.push(b'.');
		} else {
			push_int(out, i64::from(value));
		}
	}
}

pub(crate) fn push_float_list(out: &mut Vec<u8>, values: impl Iterator<Item = f32>) {
	for (index, value) in values.enumerate() {
		if is_end_float(value) {
			break;
		}
		if index > 0 {
			out.push(b',');
		}
		if is_missing_float(value) {
			out.push(b'.');
		} else {
			push_float(out, value);
		}
	}
}

pub(crate) fn push_text(out: &mut Vec<u8>, bytes: &[u8]) {
	let end = memchr::memchr(0, bytes).unwrap_or(bytes.len());
	out.extend_from_slice(bytes.get(..end).unwrap_or_default());
}

pub(crate) fn push_sample_text(out: &mut Vec<u8>, bytes: &[u8]) {
	let before = out.len();
	push_text(out, bytes);
	if out.len() == before || out.get(before..) == Some(&[MISSING_STRING_BYTE]) {
		out.truncate(before);
		out.push(b'.');
	}
}

fn write_samples(record: &Record, header: &Header, out: &mut Vec<u8>) -> Result<(), Error> {
	if record.format.is_empty() {
		out.extend_from_slice(b"\t.");
		for _ in 0..header.sample_count() {
			out.extend_from_slice(b"\t.");
		}
		return Ok(());
	}
	for (index, field) in record.format.iter().enumerate() {
		out.push(if index == 0 { b'\t' } else { b':' });
		out.extend_from_slice(tag_name(header, field.key)?.as_bytes());
	}
	for sample in 0..header.sample_count() {
		out.push(b'\t');
		for (index, field) in record.format.iter().enumerate() {
			if index > 0 {
				out.push(b':');
			}
			write_sample_field(record, header, field, sample, out)?;
		}
	}
	Ok(())
}

fn write_sample_field(
	record: &Record,
	header: &Header,
	field: &Field,
	sample: usize,
	out: &mut Vec<u8>,
) -> Result<(), Error> {
	let start = sample
		.checked_mul(field.width)
		.ok_or(Error::Bcf("FORMAT field is too large"))?;
	let end = start + field.width;
	match &field.cells {
		Cells::Ints(_) => {
			let values = record
				.store
				.ints(&field.cells)
				.get(start..end)
				.ok_or(Error::Bcf("FORMAT field is truncated"))?;
			if header.name(field.key) == Some(GENOTYPE_KEY) {
				write_genotype(out, values.iter().copied(), header.version());
			} else {
				push_int_list(out, values.iter().copied());
			}
		}
		Cells::Floats(_) => {
			let values = record
				.store
				.floats(&field.cells)
				.get(start..end)
				.ok_or(Error::Bcf("FORMAT field is truncated"))?;
			push_float_list(out, values.iter().copied());
		}
		Cells::Text(_) => {
			let values = record
				.store
				.text(&field.cells)
				.get(start..end)
				.ok_or(Error::Bcf("FORMAT field is truncated"))?;
			push_sample_text(out, values);
		}
		Cells::Flag => return Err(Error::Bcf("FORMAT field without values")),
	}
	Ok(())
}

pub(crate) fn write_genotype(
	out: &mut Vec<u8>,
	values: impl Iterator<Item = i32>,
	version: Version,
) {
	let begin = out.len();
	let mut first = 0;
	let mut any_unphased = false;
	let mut ploidy = 0usize;
	for (index, value) in values.enumerate() {
		if value == END_INT {
			break;
		}
		if index == 0 {
			first = value;
		} else {
			out.push(if value & 1 == 1 { b'|' } else { b'/' });
			any_unphased |= value & 1 == 0;
		}
		if value == MISSING_INT || value >> 1 == 0 {
			out.push(b'.');
		} else {
			push_int(out, i64::from((value >> 1) - 1));
		}
		ploidy += 1;
	}
	if ploidy == 0 {
		out.push(b'.');
	}
	if version < PREFIXED_PHASING {
		return;
	}
	let prefix = if first & 1 == 1 {
		((ploidy > 1 && any_unphased) || (ploidy <= 1 && first >> 1 == 0)).then_some(b'|')
	} else {
		((ploidy <= 1 && first != 0) || (ploidy > 1 && !any_unphased)).then_some(b'/')
	};
	if let Some(prefix) = prefix {
		out.insert(begin, prefix);
	}
}
