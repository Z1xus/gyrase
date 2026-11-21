use crate::cram::Error;
use crate::cram::decode::{Environment, reference_span};
use crate::cram::records::{
	FLAG_FIRST, FLAG_MATE_REVERSE, FLAG_MATE_UNMAPPED, FLAG_PAIRED, FLAG_REVERSE, FLAG_UNMAPPED,
	MATE_REVERSE, MATE_UNMAPPED, Record, RecordTable, Span,
};
use crate::{RecordError, RecordRef, Tags};
use std::io::Write;

const MAX_NAME_PREFIX: usize = 233;
const CORE_LEN: usize = 32;
const NIBBLES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
const LONG_CIGAR_TAG: [u8; 2] = *b"CG";

#[derive(Debug, Default)]
pub(crate) struct Decoded {
	bytes: Vec<u8>,
	ends: Vec<usize>,
}

impl Decoded {
	pub(crate) fn len(&self) -> usize {
		self.ends.len()
	}

	pub(crate) fn record(&self, index: usize) -> Result<RecordRef<'_>, RecordError> {
		let end = *self.ends.get(index).ok_or(RecordError::Truncated)?;
		let start = index
			.checked_sub(1)
			.and_then(|previous| self.ends.get(previous))
			.copied()
			.unwrap_or(0);
		let bytes = self.bytes.get(start..end).ok_or(RecordError::Truncated)?;
		RecordRef::new(bytes)
	}
}

impl Decoded {
	pub(crate) fn retain(&mut self, mut keep: impl FnMut(RecordRef<'_>) -> bool) {
		let mut start = 0;
		let mut written = 0;
		let ends = std::mem::take(&mut self.ends);
		for end in ends {
			let range = start..end;
			start = end;
			let kept = self
				.bytes
				.get(range.clone())
				.and_then(|bytes| RecordRef::new(bytes).ok())
				.is_some_and(&mut keep);
			if kept {
				let length = range.len();
				self.bytes.copy_within(range, written);
				written += length;
				self.ends.push(written);
			}
		}
		self.bytes.truncate(written);
	}
}

fn malformed() -> Error {
	Error::Malformed("record")
}

fn record_at(records: &[Record], index: usize) -> Result<&Record, Error> {
	records.get(index).ok_or_else(malformed)
}

fn record_at_mut(records: &mut [Record], index: usize) -> Result<&mut Record, Error> {
	records.get_mut(index).ok_or_else(malformed)
}

pub(crate) fn finish(
	mut table: RecordTable,
	environment: &Environment,
	record_counter: i64,
) -> Result<Decoded, Error> {
	resolve_mates(&mut table.records)?;
	name_records(&mut table, environment, record_counter)?;
	write(&table, environment)
}

fn resolve_mates(records: &mut [Record]) -> Result<(), Error> {
	let count = records.len();
	for index in 0..count {
		if let Some(mate) = record_at(records, index)?.mate_line {
			link_mate(records, index, mate)?;
		} else {
			detach_mate(record_at_mut(records, index)?);
		}
		let record = record_at_mut(records, index)?;
		if record.template_length.is_none() {
			record.template_length = Some(0);
		}
	}
	for record in records {
		if record.explicit_template_length.is_some() {
			record.template_length = record.explicit_template_length;
		}
	}
	Ok(())
}

fn link_mate(records: &mut [Record], index: usize, mate: usize) -> Result<(), Error> {
	if mate >= records.len() {
		return Err(malformed());
	}
	if record_at(records, index)?.template_length.is_none() {
		resolve_chain(records, index)?;
	}
	let (position, reference, flags) = {
		let mate = record_at(records, mate)?;
		(mate.position, mate.reference, mate.flags)
	};
	let record = record_at_mut(records, index)?;
	record.mate_position = position;
	record.mate_reference = reference;
	record.flags |= FLAG_PAIRED;
	if flags & FLAG_UNMAPPED != 0 {
		record.flags |= FLAG_MATE_UNMAPPED;
		record.template_length = Some(0);
	}
	if record.flags & FLAG_UNMAPPED != 0 {
		record.template_length = Some(0);
	}
	if flags & FLAG_REVERSE != 0 {
		record.flags |= FLAG_MATE_REVERSE;
	}
	Ok(())
}

fn detach_mate(record: &mut Record) {
	if record.mate_flags & MATE_REVERSE != 0 {
		record.flags |= FLAG_PAIRED | FLAG_MATE_REVERSE;
	}
	if record.mate_flags & MATE_UNMAPPED != 0 {
		record.flags |= FLAG_MATE_UNMAPPED;
	}
	if record.flags & FLAG_PAIRED == 0 {
		record.mate_reference = -1;
	}
}

fn resolve_chain(records: &mut [Record], first: usize) -> Result<(), Error> {
	let count = records.len();
	let head = record_at(records, first)?;
	let mut reference = head.reference;
	let mut left = i64::from(head.position);
	let mut right = i64::from(head.end);
	let (mut left_count, mut right_count) = (0u32, 0u32);
	let mut cursor = first;
	for _ in 0..=count {
		let record = record_at(records, cursor)?;
		let (position, end) = (i64::from(record.position), i64::from(record.end));
		if left > position {
			left = position;
			left_count = 1;
		} else if left == position {
			left_count += 1;
		}
		if right < end {
			right = end;
			right_count = 1;
		} else if right == end {
			right_count += 1;
		}
		let Some(next) = record.mate_line else {
			record_at_mut(records, cursor)?.mate_line = Some(first);
			break;
		};
		if next <= cursor || next >= count {
			return Err(malformed());
		}
		cursor = next;
		if record_at(records, cursor)?.reference != reference {
			reference = -1;
		}
		if cursor == first {
			break;
		}
	}
	let (head_value, others) = if reference == -1 {
		(0, 0)
	} else {
		let length = i32::try_from(right - left + 1).map_err(|_| Error::IntegerOverflow)?;
		let head = record_at(records, first)?;
		let leftmost = i64::from(head.position) == left;
		let tie = i64::from(head.end) == right && left_count > 1 && right_count > 1;
		let first_segment = if tie {
			head.flags & FLAG_FIRST != 0
		} else {
			i64::from(head.end) < right || left_count <= 1
		};
		if leftmost && first_segment {
			(length, -length)
		} else {
			(-length, length)
		}
	};
	record_at_mut(records, first)?.template_length = Some(head_value);
	let mut cursor = record_at(records, first)?.mate_line;
	for _ in 0..count {
		let Some(index) = cursor.filter(|index| *index != first) else {
			break;
		};
		let record = record_at_mut(records, index)?;
		record.template_length = Some(others);
		cursor = record.mate_line;
	}
	Ok(())
}

fn name_records(
	table: &mut RecordTable,
	environment: &Environment,
	record_counter: i64,
) -> Result<(), Error> {
	let stored = table.names.len();
	let prefix = environment
		.prefix
		.get(..environment.prefix.len().min(MAX_NAME_PREFIX))
		.unwrap_or_default();
	for index in 0..table.records.len() {
		let record = record_at(&table.records, index)?;
		if !record.name.is_empty() {
			continue;
		}
		let mate = record.mate_line;
		let mate_name = mate
			.and_then(|mate| table.records.get(mate))
			.map(|mate| mate.name)
			.filter(|name| !name.is_empty() && name.end <= stored);
		let start = table.names.len();
		if let Some(name) = mate_name {
			table.names.extend_from_within(name.start..name.end);
		} else {
			let number = mate.filter(|mate| *mate < index).unwrap_or(index);
			let number = record_counter
				.checked_add(i64::try_from(number).map_err(|_| Error::IntegerOverflow)?)
				.and_then(|number| number.checked_add(1))
				.ok_or(Error::IntegerOverflow)?;
			table.names.extend_from_slice(prefix);
			write!(table.names, ":{number}")?;
		}
		record_at_mut(&mut table.records, index)?.name = Span::since(start, table.names.len());
	}
	Ok(())
}

fn long_placeholder(length: usize, cigar: &[u32]) -> Result<[u32; 2], Error> {
	let op = |len: u64, kind: u32| {
		u32::try_from(len)
			.ok()
			.filter(|len| len >> 28 == 0)
			.map(|len| len << 4 | kind)
			.ok_or(Error::IntegerOverflow)
	};
	Ok([
		op(
			u64::try_from(length).map_err(|_| Error::IntegerOverflow)?,
			4,
		)?,
		op(reference_span(cigar), 3)?,
	])
}

const fn nibble_table() -> [u8; 256] {
	let mut table = [15u8; 256];
	let mut code = 0u8;
	while code < 16 {
		let base = NIBBLES[code as usize];
		table[base as usize] = code;
		table[base.to_ascii_lowercase() as usize] = code;
		code += 1;
	}
	table
}

const NIBBLE_TABLE: [u8; 256] = nibble_table();

fn nibble(base: u8) -> u8 {
	NIBBLE_TABLE.get(usize::from(base)).copied().unwrap_or(15)
}

fn region_bin(start: i64, end: i64) -> u16 {
	let last = end - 1;
	[(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)]
		.into_iter()
		.find(|(shift, _)| start >> shift == last >> shift)
		.and_then(|(shift, offset)| u16::try_from(offset + (start >> shift)).ok())
		.unwrap_or(0)
}

fn write(table: &RecordTable, environment: &Environment) -> Result<Decoded, Error> {
	let mut out = Decoded {
		bytes: Vec::with_capacity(
			table.records.len() * (CORE_LEN + 64)
				+ table.bases.len() * 3 / 2
				+ table.tags.len()
				+ table.names.len(),
		),
		ends: Vec::with_capacity(table.records.len()),
	};
	for record in &table.records {
		let name = record.name.of(&table.names);
		let cigar = record.cigar.of(&table.cigars);
		let bytes = &mut out.bytes;
		bytes.extend_from_slice(&record.reference.to_le_bytes());
		bytes.extend_from_slice(
			&record
				.position
				.checked_sub(1)
				.ok_or(Error::IntegerOverflow)?
				.to_le_bytes(),
		);
		bytes.push(u8::try_from(name.len() + 1).map_err(|_| malformed())?);
		bytes.push(record.mapq);
		let bin = region_bin(i64::from(record.position) - 1, i64::from(record.end));
		bytes.extend_from_slice(&bin.to_le_bytes());
		let long = u16::try_from(cigar.len()).is_err();
		let placeholder;
		let stored = if long {
			placeholder = long_placeholder(record.length, cigar)?;
			&placeholder[..]
		} else {
			cigar
		};
		let operations = u16::try_from(stored.len()).map_err(|_| Error::IntegerOverflow)?;
		bytes.extend_from_slice(&operations.to_le_bytes());
		bytes.extend_from_slice(&record.flags.to_le_bytes());
		let length = u32::try_from(record.length).map_err(|_| Error::IntegerOverflow)?;
		bytes.extend_from_slice(&length.to_le_bytes());
		bytes.extend_from_slice(&record.mate_reference.to_le_bytes());
		bytes.extend_from_slice(
			&record
				.mate_position
				.checked_sub(1)
				.ok_or(Error::IntegerOverflow)?
				.to_le_bytes(),
		);
		bytes.extend_from_slice(&record.template_length.unwrap_or(0).to_le_bytes());
		bytes.extend_from_slice(name);
		bytes.push(0);
		for op in stored {
			bytes.extend_from_slice(&op.to_le_bytes());
		}
		let bases = record.bases.of(&table.bases);
		let (pairs, odd) = bases.as_chunks::<2>();
		bytes.extend(
			pairs
				.iter()
				.map(|[high, low]| (nibble(*high) << 4) | nibble(*low)),
		);
		if let Some(last) = odd.first() {
			bytes.push(nibble(*last) << 4);
		}
		bytes.extend_from_slice(record.qualities.of(&table.qualities));
		let tags = record.tags.of(&table.tags);
		if long {
			for field in &Tags::new(tags) {
				let field = field?;
				if field.tag != LONG_CIGAR_TAG {
					field.value.encode(field.tag, bytes);
				}
			}
		} else {
			bytes.extend_from_slice(tags);
		}
		if record.read_group >= 0 {
			let group = usize::try_from(record.read_group)
				.ok()
				.and_then(|index| environment.read_groups.get(index))
				.ok_or_else(malformed)?;
			bytes.extend_from_slice(b"RGZ");
			bytes.extend_from_slice(group);
			bytes.push(0);
		}
		if long {
			bytes.extend_from_slice(b"CGBI");
			let count = u32::try_from(cigar.len()).map_err(|_| Error::IntegerOverflow)?;
			bytes.extend_from_slice(&count.to_le_bytes());
			for op in cigar {
				bytes.extend_from_slice(&op.to_le_bytes());
			}
		}
		out.ends.push(bytes.len());
	}
	Ok(out)
}
