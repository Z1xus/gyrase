use std::marker::PhantomData;

use gyrase_genome::{Code4, ContigId};

use super::bin::region_bin;
use super::cigar::{CigarKind, CigarOp, MAX_OP_LEN};
use super::fixed::{self, CORE_LEN, Core, Layout};
use super::flags::Flags;
use super::tags::{Element, Value};
use super::view::RecordRef;
use crate::error::RecordError;

const MAX_NAME_LEN: usize = 254;
const LONG_CIGAR_TAG: [u8; 2] = *b"CG";
const MISSING_QUALITY: u8 = 0xff;
const ASCII_QUALITY_OFFSET: u8 = 33;
const ASCII_QUALITY_MAX: u8 = 126;

#[derive(Debug, Clone, Copy)]
pub enum QualityInput<'a> {
	Missing,
	Phred(&'a [u8]),
	Ascii(&'a [u8]),
}

#[derive(Debug, Clone)]
pub struct Record {
	core: Core,
	tail: Vec<u8>,
	layout: Layout,
}

impl Default for Record {
	fn default() -> Self {
		Self::new()
	}
}

impl Record {
	pub fn new() -> Self {
		let mut record = Self {
			core: [0; CORE_LEN],
			tail: Vec::new(),
			layout: Layout {
				name_end: 0,
				cigar_end: 0,
				seq_end: 0,
				qual_end: 0,
				seq_len: 0,
			},
		};
		record.clear();
		record
	}

	pub fn clear(&mut self) {
		self.core = [0; CORE_LEN];
		fixed::set_contig(&mut self.core, -1);
		fixed::set_position(&mut self.core, -1);
		fixed::set_mate_contig(&mut self.core, -1);
		fixed::set_mate_position(&mut self.core, -1);
		fixed::set_flags(&mut self.core, Flags::UNMAPPED.bits());
		fixed::set_name_len(&mut self.core, 2);
		self.tail.clear();
		self.tail.extend_from_slice(b"*\0");
		self.layout = Layout {
			name_end: 2,
			cigar_end: 2,
			seq_end: 2,
			qual_end: 2,
			seq_len: 0,
		};
		self.refresh_bin();
	}

	pub fn assign(&mut self, record: RecordRef<'_>) {
		let (core, rest) = record.parts();
		self.core = *core;
		self.tail.clear();
		self.tail.extend_from_slice(rest);
		self.layout = Layout::from_counts(
			fixed::name_len(core),
			fixed::cigar_ops(core),
			fixed::seq_len(core),
		)
		.unwrap_or(self.layout);
	}

	pub fn view(&self) -> RecordRef<'_> {
		RecordRef::from_parts(&self.core, &self.tail, self.layout)
	}

	pub fn set_name(&mut self, name: &[u8]) -> Result<(), RecordError> {
		if name.is_empty() {
			return Err(RecordError::InvalidName);
		}
		if name.len() > MAX_NAME_LEN || name.contains(&0) {
			return Err(RecordError::NameTooLong);
		}
		let name_len = name.len() + 1;
		let encoded = u8::try_from(name_len).map_err(|_| RecordError::NameTooLong)?;
		self.tail
			.splice(
				..self.layout.name_end,
				name.iter().copied().chain(std::iter::once(0)),
			)
			.for_each(drop);
		fixed::set_name_len(&mut self.core, encoded);
		self.relayout(name_len, fixed::cigar_ops(&self.core), self.layout.seq_len)
	}

	pub fn set_flags(&mut self, flags: Flags) {
		fixed::set_flags(&mut self.core, flags.bits());
		self.refresh_bin();
	}

	pub fn set_contig(&mut self, contig: Option<ContigId>) -> Result<(), RecordError> {
		fixed::set_contig(&mut self.core, contig_raw(contig)?);
		Ok(())
	}

	pub fn set_position(&mut self, position: Option<u32>) -> Result<(), RecordError> {
		fixed::set_position(&mut self.core, position_raw(position)?);
		self.refresh_bin();
		Ok(())
	}

	pub fn set_mapq(&mut self, mapq: u8) {
		fixed::set_mapq(&mut self.core, mapq);
	}

	pub fn set_mate_contig(&mut self, contig: Option<ContigId>) -> Result<(), RecordError> {
		fixed::set_mate_contig(&mut self.core, contig_raw(contig)?);
		Ok(())
	}

	pub fn set_mate_position(&mut self, position: Option<u32>) -> Result<(), RecordError> {
		fixed::set_mate_position(&mut self.core, position_raw(position)?);
		Ok(())
	}

	pub fn set_template_length(&mut self, length: i32) {
		fixed::set_template_len(&mut self.core, length);
	}

	pub fn set_cigar(&mut self, ops: &[CigarOp]) -> Result<(), RecordError> {
		let Ok(count) = u16::try_from(ops.len()) else {
			return self.set_long_cigar(ops);
		};
		if ops.iter().any(|op| op.len() > MAX_OP_LEN) {
			return Err(RecordError::CigarLengthTooLarge);
		}
		if self.view().has_long_cigar() {
			self.remove_tags(LONG_CIGAR_TAG);
		}
		let encoded = ops.iter().flat_map(|op| op.raw().to_le_bytes());
		self.tail
			.splice(self.layout.name_end..self.layout.cigar_end, encoded)
			.for_each(drop);
		fixed::set_cigar_ops(&mut self.core, count);
		self.relayout(
			self.layout.name_end,
			usize::from(count),
			self.layout.seq_len,
		)?;
		self.refresh_bin();
		Ok(())
	}

	fn set_long_cigar(&mut self, ops: &[CigarOp]) -> Result<(), RecordError> {
		let reference: u64 = ops
			.iter()
			.filter(|op| op.kind().consumes_reference())
			.map(|op| u64::from(op.len()))
			.sum();
		let placeholder = [
			u32::try_from(self.layout.seq_len)
				.ok()
				.and_then(|len| CigarOp::new(len, CigarKind::SoftClip)),
			u32::try_from(reference)
				.ok()
				.and_then(|len| CigarOp::new(len, CigarKind::Skip)),
		];
		let [Some(clip), Some(skip)] = placeholder else {
			return Err(RecordError::CigarLengthTooLarge);
		};
		self.remove_tags(LONG_CIGAR_TAG);
		self.set_cigar(&[clip, skip])?;
		let mut array = self.begin_array::<u32>(LONG_CIGAR_TAG);
		for op in ops {
			array.push(op.raw())?;
		}
		Ok(())
	}

	pub fn set_sequence(
		&mut self,
		bases: &[u8],
		qualities: QualityInput<'_>,
	) -> Result<(), RecordError> {
		let seq_len = u32::try_from(bases.len())
			.ok()
			.filter(|&len| i32::try_from(len).is_ok())
			.ok_or(RecordError::SequenceTooLong)?;
		check_qualities(bases.len(), qualities)?;
		let long_cigar = self.view().has_long_cigar();
		let clip = CigarOp::new(seq_len, CigarKind::SoftClip);
		if long_cigar && clip.is_none() {
			return Err(RecordError::CigarLengthTooLarge);
		}
		let packed = bases.chunks(2).map(|pair| match *pair {
			[high, low] => Code4::from_ascii(high).code() << 4 | Code4::from_ascii(low).code(),
			[high] => Code4::from_ascii(high).code() << 4,
			_ => 0,
		});
		let phred = (0..bases.len()).map(|index| quality_byte(qualities, index));
		self.tail
			.splice(
				self.layout.cigar_end..self.layout.qual_end,
				packed.chain(phred),
			)
			.for_each(drop);
		fixed::set_seq_len(&mut self.core, seq_len);
		self.relayout(
			self.layout.name_end,
			fixed::cigar_ops(&self.core),
			bases.len(),
		)?;
		if let (true, Some(clip)) = (long_cigar, clip) {
			let at = self.layout.name_end;
			if let Some(slot) = self.tail.get_mut(at..at + 4) {
				slot.copy_from_slice(&clip.raw().to_le_bytes());
			}
		}
		Ok(())
	}

	fn remove_tags(&mut self, tag: [u8; 2]) {
		let start = self.layout.qual_end;
		let total = self.tail.len() - start;
		let mut found = Vec::new();
		let mut fields = self.view().tags().iter();
		loop {
			let from = total - fields.remaining();
			match fields.next() {
				Some(Ok(field)) if field.tag == tag => {
					found.push(start + from..start + total - fields.remaining());
				}
				Some(Ok(_)) => {}
				_ => break,
			}
		}
		for range in found.into_iter().rev() {
			self.tail.drain(range);
		}
	}

	pub fn clear_tags(&mut self) {
		self.tail.truncate(self.layout.qual_end);
	}

	pub fn push_tag(&mut self, tag: [u8; 2], value: Value<'_>) -> Result<(), RecordError> {
		if let Value::String(text) | Value::Hex(text) = value
			&& text.iter().any(|&byte| byte < b' ')
		{
			return Err(RecordError::InvalidTagString);
		}
		value.encode(tag, &mut self.tail);
		Ok(())
	}

	pub fn push_integer(&mut self, tag: [u8; 2], value: i64) -> Result<(), RecordError> {
		let encoded = if value < 0 {
			if let Ok(small) = i8::try_from(value) {
				Value::Int8(small)
			} else if let Ok(small) = i16::try_from(value) {
				Value::Int16(small)
			} else {
				Value::Int32(i32::try_from(value).map_err(|_| RecordError::IntegerOutOfRange)?)
			}
		} else if let Ok(small) = u8::try_from(value) {
			Value::UInt8(small)
		} else if let Ok(small) = u16::try_from(value) {
			Value::UInt16(small)
		} else {
			Value::UInt32(u32::try_from(value).map_err(|_| RecordError::IntegerOutOfRange)?)
		};
		self.push_tag(tag, encoded)
	}

	pub fn begin_array<T: Element>(&mut self, tag: [u8; 2]) -> ArrayWriter<'_, T> {
		self.tail.extend_from_slice(&tag);
		self.tail.extend_from_slice(&[b'B', T::CODE]);
		let count_at = self.tail.len();
		self.tail.extend_from_slice(&0u32.to_le_bytes());
		ArrayWriter {
			tail: &mut self.tail,
			count_at,
			count: 0,
			element: PhantomData,
		}
	}

	fn relayout(
		&mut self,
		name_len: usize,
		cigar_ops: usize,
		seq_len: usize,
	) -> Result<(), RecordError> {
		self.layout =
			Layout::from_counts(name_len, cigar_ops, seq_len).ok_or(RecordError::Truncated)?;
		Ok(())
	}

	fn refresh_bin(&mut self) {
		let start = i64::from(fixed::position(&self.core));
		let span = if Flags::new(fixed::flags(&self.core)).is_unmapped() {
			1
		} else {
			i64::try_from(self.view().cigar().reference_len())
				.unwrap_or(1)
				.max(1)
		};
		fixed::set_bin(&mut self.core, region_bin(start, start + span));
	}
}

#[derive(Debug)]
pub struct ArrayWriter<'a, T: Element> {
	tail: &'a mut Vec<u8>,
	count_at: usize,
	count: u32,
	element: PhantomData<T>,
}

impl<T: Element> ArrayWriter<'_, T> {
	pub fn push(&mut self, value: T) -> Result<(), RecordError> {
		self.count = self.count.checked_add(1).ok_or(RecordError::ArrayTooLong)?;
		value.write(self.tail);
		if let Some(slot) = self.tail.get_mut(self.count_at..self.count_at + 4) {
			slot.copy_from_slice(&self.count.to_le_bytes());
		}
		Ok(())
	}
}

fn contig_raw(contig: Option<ContigId>) -> Result<i32, RecordError> {
	contig.map_or(Ok(-1), |id| {
		i32::try_from(id.get()).map_err(|_| RecordError::InvalidCoordinate)
	})
}

fn position_raw(position: Option<u32>) -> Result<i32, RecordError> {
	position.map_or(Ok(-1), |value| {
		i32::try_from(value).map_err(|_| RecordError::InvalidCoordinate)
	})
}

fn check_qualities(len: usize, qualities: QualityInput<'_>) -> Result<(), RecordError> {
	match qualities {
		QualityInput::Missing => Ok(()),
		QualityInput::Phred(values) if values.len() == len => Ok(()),
		QualityInput::Ascii(values) if values.len() == len => values
			.iter()
			.all(|&value| (ASCII_QUALITY_OFFSET..=ASCII_QUALITY_MAX).contains(&value))
			.then_some(())
			.ok_or(RecordError::InvalidQuality),
		_ => Err(RecordError::LengthMismatch),
	}
}

fn quality_byte(qualities: QualityInput<'_>, index: usize) -> u8 {
	match qualities {
		QualityInput::Missing => MISSING_QUALITY,
		QualityInput::Phred(values) => values.get(index).copied().unwrap_or(MISSING_QUALITY),
		QualityInput::Ascii(values) => values.get(index).map_or(MISSING_QUALITY, |value| {
			value.wrapping_sub(ASCII_QUALITY_OFFSET)
		}),
	}
}
