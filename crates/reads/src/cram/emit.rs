use crate::cram::Error;
use crate::cram::decode::Environment;
use crate::cram::records::{
	FLAG_FIRST, FLAG_MATE_REVERSE, FLAG_MATE_UNMAPPED, FLAG_PAIRED, FLAG_REVERSE, FLAG_UNMAPPED,
	MATE_REVERSE, MATE_UNMAPPED, Rec, Records, Span,
};
use crate::{RecordError, RecordRef};
use std::io::Write;

const MAX_NAME_PREFIX: usize = 233;
const CORE_LEN: usize = 32;
const NIBBLES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

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

fn record_at(recs: &[Rec], index: usize) -> Result<&Rec, Error> {
	recs.get(index).ok_or_else(malformed)
}

fn record_at_mut(recs: &mut [Rec], index: usize) -> Result<&mut Rec, Error> {
	recs.get_mut(index).ok_or_else(malformed)
}

pub(crate) fn finish(
	mut records: Records,
	environment: &Environment,
	record_counter: i64,
) -> Result<Decoded, Error> {
	resolve_mates(&mut records.recs)?;
	name_records(&mut records, environment, record_counter)?;
	write(&records, environment)
}

fn resolve_mates(recs: &mut [Rec]) -> Result<(), Error> {
	let count = recs.len();
	for index in 0..count {
		if let Some(mate) = record_at(recs, index)?.mate_line {
			link_mate(recs, index, mate)?;
		} else {
			detach_mate(record_at_mut(recs, index)?);
		}
		let rec = record_at_mut(recs, index)?;
		if rec.template_length.is_none() {
			rec.template_length = Some(0);
		}
	}
	for rec in recs {
		if rec.explicit_template_length.is_some() {
			rec.template_length = rec.explicit_template_length;
		}
	}
	Ok(())
}

fn link_mate(recs: &mut [Rec], index: usize, mate: usize) -> Result<(), Error> {
	if mate >= recs.len() {
		return Err(malformed());
	}
	if record_at(recs, index)?.template_length.is_none() {
		resolve_chain(recs, index)?;
	}
	let (position, reference, flags) = {
		let mate = record_at(recs, mate)?;
		(mate.position, mate.reference, mate.flags)
	};
	let rec = record_at_mut(recs, index)?;
	rec.mate_position = position;
	rec.mate_reference = reference;
	rec.flags |= FLAG_PAIRED;
	if flags & FLAG_UNMAPPED != 0 {
		rec.flags |= FLAG_MATE_UNMAPPED;
		rec.template_length = Some(0);
	}
	if rec.flags & FLAG_UNMAPPED != 0 {
		rec.template_length = Some(0);
	}
	if flags & FLAG_REVERSE != 0 {
		rec.flags |= FLAG_MATE_REVERSE;
	}
	Ok(())
}

fn detach_mate(rec: &mut Rec) {
	if rec.mate_flags & MATE_REVERSE != 0 {
		rec.flags |= FLAG_PAIRED | FLAG_MATE_REVERSE;
	}
	if rec.mate_flags & MATE_UNMAPPED != 0 {
		rec.flags |= FLAG_MATE_UNMAPPED;
	}
	if rec.flags & FLAG_PAIRED == 0 {
		rec.mate_reference = -1;
	}
}

fn resolve_chain(recs: &mut [Rec], first: usize) -> Result<(), Error> {
	let count = recs.len();
	let head = record_at(recs, first)?;
	let mut reference = head.reference;
	let mut left = i64::from(head.position);
	let mut right = i64::from(head.end);
	let (mut left_count, mut right_count) = (0u32, 0u32);
	let mut cursor = first;
	for _ in 0..=count {
		let rec = record_at(recs, cursor)?;
		let (position, end) = (i64::from(rec.position), i64::from(rec.end));
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
		let Some(next) = rec.mate_line else {
			record_at_mut(recs, cursor)?.mate_line = Some(first);
			break;
		};
		if next <= cursor || next >= count {
			return Err(malformed());
		}
		cursor = next;
		if record_at(recs, cursor)?.reference != reference {
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
		let head = record_at(recs, first)?;
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
	record_at_mut(recs, first)?.template_length = Some(head_value);
	let mut cursor = record_at(recs, first)?.mate_line;
	for _ in 0..count {
		let Some(index) = cursor.filter(|index| *index != first) else {
			break;
		};
		let rec = record_at_mut(recs, index)?;
		rec.template_length = Some(others);
		cursor = rec.mate_line;
	}
	Ok(())
}

fn name_records(
	records: &mut Records,
	environment: &Environment,
	record_counter: i64,
) -> Result<(), Error> {
	let stored = records.names.len();
	let prefix = environment
		.prefix
		.get(..environment.prefix.len().min(MAX_NAME_PREFIX))
		.unwrap_or_default();
	for index in 0..records.recs.len() {
		let rec = record_at(&records.recs, index)?;
		if !rec.name.is_empty() {
			continue;
		}
		let mate = rec.mate_line;
		let mate_name = mate
			.and_then(|mate| records.recs.get(mate))
			.map(|mate| mate.name)
			.filter(|name| !name.is_empty() && name.end <= stored);
		let start = records.names.len();
		if let Some(name) = mate_name {
			records.names.extend_from_within(name.start..name.end);
		} else {
			let number = mate.filter(|mate| *mate < index).unwrap_or(index);
			let number = record_counter
				.checked_add(i64::try_from(number).map_err(|_| Error::IntegerOverflow)?)
				.and_then(|number| number.checked_add(1))
				.ok_or(Error::IntegerOverflow)?;
			records.names.extend_from_slice(prefix);
			write!(records.names, ":{number}")?;
		}
		record_at_mut(&mut records.recs, index)?.name = Span::since(start, records.names.len());
	}
	Ok(())
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

fn write(records: &Records, environment: &Environment) -> Result<Decoded, Error> {
	let mut out = Decoded {
		bytes: Vec::with_capacity(
			records.recs.len() * (CORE_LEN + 64)
				+ records.bases.len() * 3 / 2
				+ records.tags.len()
				+ records.names.len(),
		),
		ends: Vec::with_capacity(records.recs.len()),
	};
	for rec in &records.recs {
		let name = rec.name.of(&records.names);
		let cigar = rec.cigar.of(&records.cigars);
		let bytes = &mut out.bytes;
		bytes.extend_from_slice(&rec.reference.to_le_bytes());
		bytes.extend_from_slice(&(rec.position - 1).to_le_bytes());
		bytes.push(u8::try_from(name.len() + 1).map_err(|_| malformed())?);
		bytes.push(rec.mapq);
		let bin = region_bin(i64::from(rec.position) - 1, i64::from(rec.end));
		bytes.extend_from_slice(&bin.to_le_bytes());
		let operations = u16::try_from(cigar.len()).map_err(|_| Error::IntegerOverflow)?;
		bytes.extend_from_slice(&operations.to_le_bytes());
		bytes.extend_from_slice(&rec.flags.to_le_bytes());
		let length = u32::try_from(rec.length).map_err(|_| Error::IntegerOverflow)?;
		bytes.extend_from_slice(&length.to_le_bytes());
		bytes.extend_from_slice(&rec.mate_reference.to_le_bytes());
		bytes.extend_from_slice(&(rec.mate_position - 1).to_le_bytes());
		bytes.extend_from_slice(&rec.template_length.unwrap_or(0).to_le_bytes());
		bytes.extend_from_slice(name);
		bytes.push(0);
		for op in cigar {
			bytes.extend_from_slice(&op.to_le_bytes());
		}
		let bases = rec.bases.of(&records.bases);
		let (pairs, odd) = bases.as_chunks::<2>();
		bytes.extend(
			pairs
				.iter()
				.map(|[high, low]| (nibble(*high) << 4) | nibble(*low)),
		);
		if let Some(last) = odd.first() {
			bytes.push(nibble(*last) << 4);
		}
		bytes.extend_from_slice(rec.qualities.of(&records.qualities));
		bytes.extend_from_slice(rec.tags.of(&records.tags));
		if rec.read_group >= 0 {
			let group = usize::try_from(rec.read_group)
				.ok()
				.and_then(|index| environment.read_groups.get(index))
				.ok_or_else(malformed)?;
			bytes.extend_from_slice(b"RGZ");
			bytes.extend_from_slice(group);
			bytes.push(0);
		}
		out.ends.push(bytes.len());
	}
	Ok(out)
}
