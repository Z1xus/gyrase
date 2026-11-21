use gyrase_genome::ContigId;

use super::number::{push_general, push_short, push_signed, push_unsigned};
use crate::error::RecordError;
use crate::header::Header;
use crate::record::{Array, Cigar, Element, Elements, Field, RecordRef, Value};

const LONG_CIGAR_TAG: [u8; 2] = *b"CG";

pub(super) fn format_record(
	out: &mut Vec<u8>,
	header: &Header,
	record: RecordRef<'_>,
) -> Result<(), RecordError> {
	out.extend_from_slice(record.name());
	out.push(b'\t');
	push_unsigned(out, u64::from(record.flags().bits()));
	out.push(b'\t');
	push_contig(out, header, record.contig())?;
	out.push(b'\t');
	push_unsigned(
		out,
		record
			.position()
			.map_or(0, |position| u64::from(position) + 1),
	);
	out.push(b'\t');
	push_unsigned(out, u64::from(record.mapq()));
	out.push(b'\t');
	push_cigar(out, record.cigar());
	out.push(b'\t');
	push_mate_contig(out, header, record)?;
	out.push(b'\t');
	push_unsigned(
		out,
		record
			.mate_position()
			.map_or(0, |position| u64::from(position) + 1),
	);
	out.push(b'\t');
	push_signed(out, i64::from(record.template_length()));
	out.push(b'\t');
	push_read(out, record);
	push_tags(out, record)
}

fn push_contig(
	out: &mut Vec<u8>,
	header: &Header,
	contig: Option<ContigId>,
) -> Result<(), RecordError> {
	match contig {
		None => out.push(b'*'),
		Some(id) => {
			let reference = header
				.reference(id)
				.ok_or(RecordError::ContigOutOfRange { contig: id.get() })?;
			out.extend_from_slice(reference.name().as_bytes());
		}
	}
	Ok(())
}

fn push_mate_contig(
	out: &mut Vec<u8>,
	header: &Header,
	record: RecordRef<'_>,
) -> Result<(), RecordError> {
	match record.mate_contig() {
		Some(mate) if record.contig() == Some(mate) => out.push(b'='),
		mate => return push_contig(out, header, mate),
	}
	Ok(())
}

fn push_cigar(out: &mut Vec<u8>, cigar: Cigar<'_>) {
	if cigar.is_empty() {
		out.push(b'*');
		return;
	}
	for op in cigar.ops() {
		push_unsigned(out, u64::from(op.len()));
		out.push(op.kind().symbol());
	}
}

fn push_read(out: &mut Vec<u8>, record: RecordRef<'_>) {
	let sequence = record.sequence();
	if sequence.is_empty() {
		out.extend_from_slice(b"*\t*");
		return;
	}
	sequence.decode_into(out);
	out.push(b'\t');
	let qualities = record.qualities();
	if qualities.is_missing() {
		out.push(b'*');
	} else {
		qualities.decode_ascii_into(out);
	}
}

fn push_tags(out: &mut Vec<u8>, record: RecordRef<'_>) -> Result<(), RecordError> {
	let skip_long_cigar = record.has_long_cigar();
	for field in &record.tags() {
		let field = field?;
		if skip_long_cigar && field.tag == LONG_CIGAR_TAG {
			continue;
		}
		out.push(b'\t');
		push_field(out, field)?;
	}
	Ok(())
}

fn push_field(out: &mut Vec<u8>, field: Field<'_>) -> Result<(), RecordError> {
	out.extend_from_slice(&field.tag);
	out.push(b':');
	match field.value {
		Value::Char(value) => {
			out.extend_from_slice(b"A:");
			out.push(value);
		}
		Value::Float(value) => {
			out.extend_from_slice(b"f:");
			push_general(out, f64::from(value));
		}
		Value::Double(value) => {
			out.extend_from_slice(b"d:");
			push_general(out, value);
		}
		Value::String(text) => push_text(out, b'Z', text)?,
		Value::Hex(text) => push_text(out, b'H', text)?,
		Value::Array(array) => push_array(out, array),
		integer => {
			out.extend_from_slice(b"i:");
			push_signed(out, integer.as_integer().ok_or(RecordError::MalformedTag)?);
		}
	}
	Ok(())
}

fn push_text(out: &mut Vec<u8>, kind: u8, text: &[u8]) -> Result<(), RecordError> {
	if text.iter().any(|&byte| byte < b' ') {
		return Err(RecordError::InvalidTagString);
	}
	out.push(kind);
	out.push(b':');
	out.extend_from_slice(text);
	Ok(())
}

fn push_array(out: &mut Vec<u8>, array: Array<'_>) {
	out.extend_from_slice(b"B:");
	match array {
		Array::Int8(items) => push_items(out, b'c', items, |out, value| {
			push_signed(out, i64::from(value));
		}),
		Array::UInt8(items) => push_items(out, b'C', items, |out, value| {
			push_unsigned(out, u64::from(value));
		}),
		Array::Int16(items) => push_items(out, b's', items, |out, value| {
			push_signed(out, i64::from(value));
		}),
		Array::UInt16(items) => push_items(out, b'S', items, |out, value| {
			push_unsigned(out, u64::from(value));
		}),
		Array::Int32(items) => push_items(out, b'i', items, |out, value| {
			push_signed(out, i64::from(value));
		}),
		Array::UInt32(items) => push_items(out, b'I', items, |out, value| {
			push_unsigned(out, u64::from(value));
		}),
		Array::Float(items) => push_items(out, b'f', items, |out, value| {
			push_short(out, f64::from(value));
		}),
	}
}

fn push_items<T: Element>(
	out: &mut Vec<u8>,
	code: u8,
	items: Elements<'_, T>,
	mut push: impl FnMut(&mut Vec<u8>, T),
) {
	out.push(code);
	for value in items.iter() {
		out.push(b',');
		push(out, value);
	}
}
