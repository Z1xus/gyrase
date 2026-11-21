use crate::Error;
use crate::header::Header;
use crate::record::{Field, Record};
use crate::value::{BcfType, Cells, END_INT, MISSING_INT};

const LONG_LENGTH: u8 = 15;
const INT8_MAX: i32 = 0x7f;
const INT8_MIN: i32 = -120;
const INT16_MAX: i32 = 0x7fff;
const INT16_MIN: i32 = -32760;

#[derive(Debug, Default)]
pub struct Encoder {
	shared: Vec<u8>,
	sample_block: Vec<u8>,
}

impl Encoder {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn encode(&mut self, record: &Record, header: &Header) -> Result<(&[u8], &[u8]), Error> {
		self.shared.clear();
		self.sample_block.clear();
		record.check(header)?;
		let alleles = u16::try_from(record.allele_count()).map_err(|_| Error::Limits("alleles"))?;
		let info = u16::try_from(record.info.len()).map_err(|_| Error::Limits("INFO fields"))?;
		let format =
			u8::try_from(record.format.len()).map_err(|_| Error::Limits("FORMAT fields"))?;
		let samples = u32::try_from(header.sample_count()).map_err(|_| Error::Limits("samples"))?;
		if samples > 0x00ff_ffff {
			return Err(Error::Limits("samples"));
		}
		let length = record.end(header) - record.start;
		let shared = &mut self.shared;
		shared.extend_from_slice(
			&i32::try_from(record.contig)
				.map_err(|_| Error::Limits("contigs"))?
				.to_le_bytes(),
		);
		shared.extend_from_slice(
			&i32::try_from(record.start)
				.map_err(|_| Error::Limits("position"))?
				.to_le_bytes(),
		);
		shared.extend_from_slice(
			&i32::try_from(length)
				.map_err(|_| Error::Limits("length"))?
				.to_le_bytes(),
		);
		let qual = record.qual.map_or(crate::MISSING_FLOAT_BITS, f32::to_bits);
		shared.extend_from_slice(&qual.to_le_bytes());
		shared.extend_from_slice(&(u32::from(alleles) << 16 | u32::from(info)).to_le_bytes());
		shared.extend_from_slice(&(u32::from(format) << 24 | samples).to_le_bytes());
		push_string(shared, &record.id);
		for allele in record.alleles() {
			push_string(shared, allele);
		}
		push_ints(
			shared,
			&record
				.filters
				.iter()
				.map(|&key| key.cast_signed())
				.collect::<Vec<_>>(),
			None,
		);
		for field in &record.info {
			push_int1(shared, field.key.cast_signed());
			push_info_value(shared, record, field);
		}
		for field in &record.format {
			push_int1(&mut self.sample_block, field.key.cast_signed());
			push_format_values(&mut self.sample_block, record, field)?;
		}
		Ok((&self.shared, &self.sample_block))
	}
}

fn push_info_value(out: &mut Vec<u8>, record: &Record, field: &Field) {
	match &field.cells {
		Cells::Flag => push_size(out, 0, BcfType::Null),
		Cells::Ints(_) => push_ints(out, record.store.ints(&field.cells), None),
		Cells::Floats(_) => {
			let values = record.store.floats(&field.cells);
			push_size(out, values.len(), BcfType::Float);
			push_float_values(out, values);
		}
		Cells::Text(_) => push_string(out, record.store.text(&field.cells)),
	}
}

fn push_format_values(out: &mut Vec<u8>, record: &Record, field: &Field) -> Result<(), Error> {
	match &field.cells {
		Cells::Flag => return Err(Error::Bcf("FORMAT field without values")),
		Cells::Ints(_) => push_ints(out, record.store.ints(&field.cells), Some(field.width)),
		Cells::Floats(_) => {
			push_size(out, field.width, BcfType::Float);
			push_float_values(out, record.store.floats(&field.cells));
		}
		Cells::Text(_) => {
			push_size(out, field.width, BcfType::Char);
			out.extend_from_slice(record.store.text(&field.cells));
		}
	}
	Ok(())
}

fn push_float_values(out: &mut Vec<u8>, values: &[f32]) {
	for value in values {
		out.extend_from_slice(&value.to_bits().to_le_bytes());
	}
}

pub(crate) fn push_size(out: &mut Vec<u8>, size: usize, bcf_type: BcfType) {
	if size < usize::from(LONG_LENGTH) {
		out.push(u8::try_from(size).unwrap_or(0) << 4 | bcf_type.code());
		return;
	}
	out.push(LONG_LENGTH << 4 | bcf_type.code());
	push_int1(out, i32::try_from(size).unwrap_or(i32::MAX));
}

pub(crate) fn push_string(out: &mut Vec<u8>, bytes: &[u8]) {
	push_size(out, bytes.len(), BcfType::Char);
	out.extend_from_slice(bytes);
}

pub(crate) fn push_int1(out: &mut Vec<u8>, value: i32) {
	match value {
		END_INT => out.extend_from_slice(&[1 << 4 | BcfType::Int8.code(), 0x81]),
		MISSING_INT => out.extend_from_slice(&[1 << 4 | BcfType::Int8.code(), 0x80]),
		_ if (INT8_MIN..=INT8_MAX).contains(&value) => {
			out.extend_from_slice(&[1 << 4 | BcfType::Int8.code(), value.to_le_bytes()[0]]);
		}
		_ if (INT16_MIN..=INT16_MAX).contains(&value) => {
			out.push(1 << 4 | BcfType::Int16.code());
			out.extend_from_slice(&value.to_le_bytes()[..2]);
		}
		_ => {
			out.push(1 << 4 | BcfType::Int32.code());
			out.extend_from_slice(&value.to_le_bytes());
		}
	}
}

pub(crate) fn push_ints(out: &mut Vec<u8>, values: &[i32], width: Option<usize>) {
	match values {
		[] => push_size(out, 0, BcfType::Null),
		[single] => push_int1(out, *single),
		_ => push_int_vector(out, values, width.unwrap_or(values.len())),
	}
}

fn push_int_vector(out: &mut Vec<u8>, values: &[i32], width: usize) {
	let mut max = i32::MIN;
	let mut min = i32::MAX;
	for &value in values {
		max = max.max(value);
		if value > END_INT {
			min = min.min(value);
		}
	}
	if max <= INT8_MAX && min >= INT8_MIN {
		push_size(out, width, BcfType::Int8);
		out.extend(values.iter().map(|&value| match value {
			END_INT => 0x81,
			MISSING_INT => 0x80,
			other => other.to_le_bytes()[0],
		}));
	} else if max <= INT16_MAX && min >= INT16_MIN {
		push_size(out, width, BcfType::Int16);
		for &value in values {
			let narrow: i16 = match value {
				END_INT => i16::MIN + 1,
				MISSING_INT => i16::MIN,
				other => i16::try_from(other).unwrap_or(i16::MIN),
			};
			out.extend_from_slice(&narrow.to_le_bytes());
		}
	} else {
		push_size(out, width, BcfType::Int32);
		for &value in values {
			out.extend_from_slice(&value.to_le_bytes());
		}
	}
}
