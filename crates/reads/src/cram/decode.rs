use std::sync::Arc;

use gyrase_genome::Reference;

use crate::cram::Error;
use crate::cram::compression::{Bases, CompressionHeader, Positions, Qualities, ReadNames};
use crate::cram::emit::{self, Decoded};
use crate::cram::encoding::{Data, Encoding};
use crate::cram::records::{
	CF_DETACHED, CF_EXPLICIT_LENGTH, CF_MATE_DOWNSTREAM, CF_NO_SEQUENCE, CF_QUALITIES,
	CIGAR_DELETION, CIGAR_HARD_CLIP, CIGAR_INSERTION, CIGAR_MATCH, CIGAR_PADDING, CIGAR_SKIP,
	CIGAR_SOFT_CLIP, FLAG_REVERSE, FLAG_UNMAPPED, Rec, Records, Span,
};
use crate::cram::reference::{References, Window};
use crate::cram::regenerate;
use crate::cram::slice::{SliceBlocks, SliceHeader};
use crate::{Header, Line};

const MAX_NAME_LEN: usize = 254;
const MAX_BASES: usize = 1 << 28;
const FIRST_RESERVE: usize = 1 << 16;
const WINDOW_CHUNK: usize = 1 << 16;
const MAX_RECORDS: usize = 1 << 24;
const MISSING_QUALITY: u8 = 255;
const FILLED_QUALITY: u8 = 30;
const DEFAULT_MAPPING_QUALITY: u8 = 40;

fn malformed() -> Error {
	Error::Malformed("record")
}

#[derive(Debug)]
pub(crate) struct Environment {
	pub(crate) references: References,
	pub(crate) read_groups: Vec<Vec<u8>>,
	pub(crate) unknown_read_group: Option<usize>,
	pub(crate) prefix: Vec<u8>,
}

impl Environment {
	pub(crate) fn new(header: &Header, reference: Option<Arc<Reference>>, prefix: &str) -> Self {
		Self::with_references(header, References::new(header, reference), prefix)
	}

	pub(crate) fn unbacked(header: &Header) -> Self {
		Self::with_references(header, References::unbacked(header), "")
	}

	fn with_references(header: &Header, references: References, prefix: &str) -> Self {
		let read_groups: Vec<Vec<u8>> = header
			.lines()
			.iter()
			.filter_map(|line| match line {
				Line::ReadGroup(fields) => {
					Some(fields.get(*b"ID").unwrap_or_default().as_bytes().to_vec())
				}
				_ => None,
			})
			.collect();
		let unknown_read_group = read_groups.len().checked_sub(1).filter(|last| {
			read_groups
				.get(*last)
				.is_some_and(|name| name == b"UNKNOWN")
		});
		Self {
			references,
			read_groups,
			unknown_read_group,
			prefix: prefix.as_bytes().to_vec(),
		}
	}
}

#[derive(Debug, Default, Clone, Copy)]
struct TagFlags {
	md: bool,
	nm: bool,
}

#[derive(Debug, Default)]
struct CigarBuilder {
	pending: Option<(u32, usize)>,
}

impl CigarBuilder {
	fn add(&mut self, cigars: &mut Vec<u32>, op: u32, len: usize) -> Result<(), Error> {
		if len == 0 {
			return Ok(());
		}
		if let Some((pending, total)) = &mut self.pending
			&& *pending == op
		{
			*total = total.checked_add(len).ok_or(Error::IntegerOverflow)?;
			return Ok(());
		}
		self.flush(cigars)?;
		self.pending = Some((op, len));
		Ok(())
	}

	fn soft_clip(&mut self, cigars: &mut Vec<u32>, len: usize) -> Result<(), Error> {
		self.flush(cigars)?;
		if len > 0 {
			cigars.push(encode_op(CIGAR_SOFT_CLIP, len)?);
		}
		Ok(())
	}

	fn flush(&mut self, cigars: &mut Vec<u32>) -> Result<(), Error> {
		if let Some((op, len)) = self.pending.take() {
			cigars.push(encode_op(op, len)?);
		}
		Ok(())
	}
}

fn encode_op(op: u32, len: usize) -> Result<u32, Error> {
	let len = u32::try_from(len)
		.ok()
		.filter(|len| *len < 1 << 28)
		.ok_or(Error::IntegerOverflow)?;
	Ok((len << 4) | op)
}

#[derive(Debug)]
struct Walk {
	read: usize,
	reference: u64,
	builder: CigarBuilder,
}

fn optional(encoding: &Encoding, data: &mut Data<'_>) -> Result<Option<i32>, Error> {
	if encoding.is_null() {
		Ok(None)
	} else {
		encoding.int(data).map(Some)
	}
}

fn natural(value: i32) -> Result<usize, Error> {
	usize::try_from(value).map_err(|_| malformed())
}

fn prepare_window(
	environment: &Environment,
	header: &CompressionHeader,
	slice: &SliceHeader,
	embedded: Option<&[u8]>,
) -> Result<Option<Window>, Error> {
	if slice.reference < 0 || (embedded.is_none() && header.bases == Bases::Explicit) {
		return Ok(None);
	}
	let index = natural(slice.reference)?;
	let length = environment
		.references
		.length(index)
		.ok_or(Error::Malformed("reference id"))?;
	let start = u64::try_from(slice.start.saturating_sub(1)).map_err(|_| malformed())?;
	let end = start
		.checked_add(u64::try_from(slice.span).map_err(|_| malformed())?)
		.ok_or(Error::IntegerOverflow)?;
	if let Some(bases) = embedded {
		return Ok(Some(Window::embedded(length, start, bases.to_vec())));
	}
	if slice.checksum != [0; 16] && environment.references.has_source() {
		let actual = environment.references.checksum(index, start, end)?;
		if actual != slice.checksum {
			return Err(Error::ReferenceChecksum(index.to_string()));
		}
	}
	environment.references.window(index, start, end).map(Some)
}

pub(crate) fn decode_slice(
	environment: &Environment,
	header: &CompressionHeader,
	blocks: &SliceBlocks<'_>,
) -> Result<Decoded, Error> {
	let records = decode_records(environment, header, blocks)?;
	emit::finish(records, environment, blocks.header.record_counter)
}

pub(crate) fn decode_extents(
	environment: &Environment,
	header: &CompressionHeader,
	blocks: &SliceBlocks<'_>,
) -> Result<Vec<(i32, i32, i32)>, Error> {
	let records = decode_records(environment, header, blocks)?;
	Ok(records
		.recs
		.iter()
		.map(|rec| (rec.reference, rec.position, rec.end))
		.collect())
}

fn decode_records(
	environment: &Environment,
	header: &CompressionHeader,
	blocks: &SliceBlocks<'_>,
) -> Result<Records, Error> {
	let slice = &blocks.header;
	let records = usize::try_from(slice.records)
		.ok()
		.filter(|count| *count <= MAX_RECORDS)
		.ok_or(Error::TooLarge {
			len: u64::from(slice.records),
		})?;
	let payload = blocks.payload(&header.slots)?;
	let window = prepare_window(environment, header, slice, payload.embedded.as_deref())?;
	let mut decoder = Decoder {
		environment,
		header,
		slice,
		data: payload.data(),
		window_reference: slice.reference,
		window,
		last_position: slice.start,
		out: Records::default(),
		bases: Vec::new(),
		qualities: Vec::new(),
		array: Vec::new(),
		reference_bases: Vec::new(),
		md: Vec::new(),
		length_budget: 0,
	};
	decoder.out.recs.reserve(records.min(FIRST_RESERVE));
	for _ in 0..records {
		decoder.record()?;
	}
	Ok(decoder.out)
}

struct Decoder<'a> {
	environment: &'a Environment,
	header: &'a CompressionHeader,
	slice: &'a SliceHeader,
	data: Data<'a>,
	window: Option<Window>,
	window_reference: i32,
	last_position: i32,
	out: Records,
	bases: Vec<u8>,
	qualities: Vec<u8>,
	array: Vec<u8>,
	reference_bases: Vec<u8>,
	md: Vec<u8>,
	length_budget: usize,
}

impl Decoder<'_> {
	fn record(&mut self) -> Result<(), Error> {
		let header = self.header;
		let series = &header.series;
		let index = self.out.recs.len();
		let mut rec = Rec {
			flags: match optional(&series.flags, &mut self.data)? {
				Some(flags) => u16::try_from(flags).map_err(|_| malformed())?,
				None => FLAG_UNMAPPED,
			},
			cram_flags: match optional(&series.compression_flags, &mut self.data)? {
				Some(flags) => u32::try_from(flags).map_err(|_| malformed())?,
				None => 0,
			},
			..Rec::default()
		};
		rec.reference = if self.slice.reference == -2 {
			optional(&series.reference_id, &mut self.data)?.unwrap_or(-1)
		} else {
			self.slice.reference
		};
		self.select_window(rec.reference)?;
		if let Some(length) = optional(&series.read_length, &mut self.data)? {
			rec.length = natural(length)?;
			self.length_budget = self.length_budget.saturating_add(rec.length);
			if self.length_budget > MAX_BASES {
				return Err(Error::TooLarge {
					len: self.length_budget as u64,
				});
			}
		}
		rec.position = self.position()?;
		rec.read_group = self.read_group()?;
		if header.read_names == ReadNames::Stored {
			rec.name = self.name()?;
		}
		self.mate(&mut rec, index)?;
		let flags = self.tags(&mut rec)?;
		let side = self.out.names.len().max(self.out.tags.len());
		if side > MAX_BASES {
			return Err(Error::TooLarge { len: side as u64 });
		}
		if rec.flags & FLAG_UNMAPPED == 0 {
			self.mapped(&mut rec, flags)?;
		} else {
			self.unmapped(&mut rec)?;
		}
		if header.qualities == Qualities::Aligned
			&& !series.quality.is_null()
			&& rec.flags & FLAG_REVERSE != 0
		{
			self.qualities.reverse();
		}
		self.store(&mut rec)?;
		self.out.recs.push(rec);
		Ok(())
	}

	fn select_window(&mut self, reference: i32) -> Result<(), Error> {
		if self.slice.reference != -2 || reference == self.window_reference {
			return Ok(());
		}
		self.window_reference = reference;
		self.window = None;
		if reference < 0 || self.header.bases == Bases::Explicit {
			return Ok(());
		}
		usize::try_from(reference)
			.ok()
			.filter(|index| *index < self.environment.references.count())
			.ok_or(Error::Malformed("reference id"))?;
		Ok(())
	}

	fn cover(&mut self, position: u64, count: usize) -> Result<(), Error> {
		if self.slice.reference != -2
			|| self.window_reference < 0
			|| self.header.bases == Bases::Explicit
		{
			return Ok(());
		}
		let end = advance(position, count.max(1))?;
		if self
			.window
			.as_ref()
			.is_some_and(|window| window.covers(position, end))
		{
			return Ok(());
		}
		let index = natural(self.window_reference)?;
		let to = end.max(advance(position, WINDOW_CHUNK)?);
		self.window = Some(self.environment.references.window(index, position, to)?);
		Ok(())
	}

	fn position(&mut self) -> Result<i32, Error> {
		let Some(value) = optional(&self.header.series.position, &mut self.data)? else {
			return Ok(self.slice.start);
		};
		let position = match self.header.positions {
			Positions::Delta => value
				.checked_add(self.last_position)
				.ok_or(Error::IntegerOverflow)?,
			Positions::Absolute => value,
		};
		self.last_position = position;
		if self.slice.reference >= 0 && position < self.slice.start {
			return Err(malformed());
		}
		Ok(position)
	}

	fn read_group(&mut self) -> Result<i32, Error> {
		let Some(group) = optional(&self.header.series.read_group, &mut self.data)? else {
			return Ok(-1);
		};
		if group == -1 || usize::try_from(group).ok() == self.environment.unknown_read_group {
			return Ok(-1);
		}
		let index = natural(group)?;
		if index >= self.environment.read_groups.len() {
			return Err(malformed());
		}
		Ok(group)
	}

	fn name(&mut self) -> Result<Span, Error> {
		let series = &self.header.series;
		let start = self.out.names.len();
		if series.name.is_null() {
			return Ok(Span::since(start, start));
		}
		series.name.array(&mut self.data, &mut self.out.names)?;
		let name = self.out.names.get(start..).unwrap_or_default();
		if name.len() > MAX_NAME_LEN || name.contains(&0) {
			return Err(malformed());
		}
		Ok(Span::since(start, self.out.names.len()))
	}

	fn mate_reference(&mut self) -> Result<i32, Error> {
		let reference = optional(&self.header.series.mate_reference, &mut self.data)?.unwrap_or(-1);
		let count = self.environment.references.count();
		if reference < -1 || usize::try_from(reference).is_ok_and(|index| index >= count) {
			return Err(malformed());
		}
		Ok(reference)
	}

	fn mate(&mut self, rec: &mut Rec, index: usize) -> Result<(), Error> {
		let series = &self.header.series;
		rec.mate_reference = -1;
		if rec.cram_flags & CF_DETACHED != 0 {
			rec.mate_flags = match optional(&series.mate_flags, &mut self.data)? {
				Some(flags) => u8::try_from(flags).map_err(|_| malformed())?,
				None => 0,
			};
			if self.header.read_names == ReadNames::Generated {
				rec.name = self.name()?;
			}
			rec.mate_reference = self.mate_reference()?;
			rec.mate_position = optional(&series.mate_position, &mut self.data)?.unwrap_or(0);
			rec.template_length = optional(&series.template_length, &mut self.data)?;
		} else if rec.cram_flags & CF_MATE_DOWNSTREAM != 0 {
			if let Some(distance) = optional(&series.next_fragment, &mut self.data)? {
				let mate = index
					.checked_add(natural(distance)?)
					.and_then(|mate| mate.checked_add(1))
					.ok_or(Error::IntegerOverflow)?;
				rec.mate_line = Some(mate);
			}
			if rec.cram_flags & CF_EXPLICIT_LENGTH != 0 {
				rec.explicit_template_length = optional(&series.template_length, &mut self.data)?;
			}
		} else if rec.cram_flags & CF_EXPLICIT_LENGTH != 0 {
			rec.explicit_template_length = optional(&series.template_length, &mut self.data)?;
		}
		Ok(())
	}

	fn tags(&mut self, rec: &mut Rec) -> Result<TagFlags, Error> {
		let header = self.header;
		let start = self.out.tags.len();
		let mut flags = TagFlags::default();
		if let Some(line) = optional(&header.series.tag_line, &mut self.data)? {
			let bindings = usize::try_from(line)
				.ok()
				.and_then(|line| header.tag_lines.get(line))
				.ok_or_else(malformed)?;
			for binding in bindings {
				let encoding = header
					.tag_encodings
					.get(binding.encoding)
					.ok_or_else(malformed)?;
				let entry = self.out.tags.len();
				self.out.tags.extend_from_slice(&binding.name);
				self.out.tags.push(binding.kind);
				let value = self.out.tags.len();
				encoding.tag_value(binding.kind, &mut self.data, &mut self.out.tags)?;
				finish_value(binding.kind, &mut self.out.tags, value)?;
				match (binding.name, binding.kind) {
					([b'c', b'F'], b'C') => {
						let bits = self.out.tags.get(value).copied().unwrap_or(0);
						self.out.tags.truncate(entry);
						flags.md |= bits & 1 != 0;
						flags.nm |= bits & 2 != 0;
					}
					([b'M', b'D'], _) => flags.md = true,
					([b'N', b'M'], _) => flags.nm = true,
					_ => {}
				}
			}
		}
		rec.tags = Span::since(start, self.out.tags.len());
		Ok(flags)
	}

	fn unmapped(&mut self, rec: &mut Rec) -> Result<(), Error> {
		let series = &self.header.series;
		rec.end = rec.position;
		self.bases.clear();
		if rec.length > 0 && rec.cram_flags & CF_NO_SEQUENCE == 0 && !series.base.is_null() {
			series
				.base
				.bytes(rec.length, &mut self.data, &mut self.bases)?;
		} else {
			self.bases.resize(rec.length, b'=');
		}
		self.qualities.clear();
		if rec.cram_flags & CF_QUALITIES != 0 && !series.quality.is_null() {
			series
				.quality
				.bytes(rec.length, &mut self.data, &mut self.qualities)?;
		} else {
			self.qualities.resize(rec.length, MISSING_QUALITY);
		}
		Ok(())
	}

	fn mapped(&mut self, rec: &mut Rec, flags: TagFlags) -> Result<(), Error> {
		let header = self.header;
		let series = &header.series;
		if rec.position <= 0 {
			return Err(malformed());
		}
		let preserved = rec.cram_flags & CF_QUALITIES != 0;
		self.bases.clear();
		self.bases.resize(rec.length, b'=');
		self.qualities.clear();
		self.qualities.resize(rec.length, MISSING_QUALITY);
		let cigar_start = self.out.cigars.len();
		let mut walk = Walk {
			read: 0,
			reference: u64::try_from(rec.position - 1).map_err(|_| malformed())?,
			builder: CigarBuilder::default(),
		};
		let count = natural(series.feature_count.int(&mut self.data)?)?;
		let mut previous = 0i32;
		for _ in 0..count {
			let code = series.feature_code.byte(&mut self.data)?;
			let delta = series.feature_position.int(&mut self.data)?;
			previous = previous.checked_add(delta).ok_or(Error::IntegerOverflow)?;
			let position = usize::try_from(previous)
				.ok()
				.filter(|position| *position > 0)
				.ok_or_else(malformed)?;
			if rec.length != 0 && position > rec.length {
				let allowed = if matches!(code, b'N' | b'P' | b'H' | b'D') {
					rec.length + 1
				} else {
					rec.length
				};
				if position > allowed {
					return Err(malformed());
				}
			}
			let index = position - 1;
			if matches!(code, b'S' | b'X' | b'I' | b'i' | b'b' | b'B') && index < walk.read {
				return Err(malformed());
			}
			let gap = index.saturating_sub(walk.read);
			self.matches(&mut walk, gap)?;
			self.feature(code, index, preserved, &mut walk)?;
		}
		let tail = rec.length.saturating_sub(walk.read);
		self.matches(&mut walk, tail)?;
		walk.builder.flush(&mut self.out.cigars)?;
		rec.cigar = Span::since(cigar_start, self.out.cigars.len());
		let span = reference_span(rec.cigar.of(&self.out.cigars));
		let last = i64::from(rec.position)
			+ i64::try_from(span.max(1)).map_err(|_| Error::IntegerOverflow)?
			- 1;
		rec.end = i32::try_from(last).map_err(|_| Error::IntegerOverflow)?;
		rec.mapq = match optional(&series.mapping_quality, &mut self.data)? {
			Some(mapq) => u8::try_from(mapq).map_err(|_| malformed())?,
			None => DEFAULT_MAPPING_QUALITY,
		};
		if preserved && !series.quality.is_null() {
			self.qualities.clear();
			series
				.quality
				.bytes(rec.length, &mut self.data, &mut self.qualities)?;
		}
		if rec.cram_flags & CF_NO_SEQUENCE == 0 && !(flags.md && flags.nm) {
			self.regenerate(rec, flags)?;
		}
		Ok(())
	}

	fn matches(&mut self, walk: &mut Walk, count: usize) -> Result<(), Error> {
		if count == 0 {
			return Ok(());
		}
		let end = walk.read.checked_add(count).ok_or(Error::IntegerOverflow)?;
		self.cover(walk.reference, count)?;
		let target = self.bases.get_mut(walk.read..end).ok_or_else(malformed)?;
		if let Some(window) = &self.window {
			window.copy(walk.reference, target)?;
		}
		walk.builder.add(&mut self.out.cigars, CIGAR_MATCH, count)?;
		walk.read = end;
		walk.reference = advance(walk.reference, count)?;
		Ok(())
	}

	fn array_into(
		encoding: &Encoding,
		data: &mut Data<'_>,
		scratch: &mut Vec<u8>,
		target: &mut [u8],
		at: usize,
	) -> Result<usize, Error> {
		scratch.clear();
		encoding.array(data, scratch)?;
		let end = at
			.checked_add(scratch.len())
			.ok_or(Error::IntegerOverflow)?;
		target
			.get_mut(at..end)
			.ok_or_else(malformed)?
			.copy_from_slice(scratch);
		Ok(scratch.len())
	}

	fn missing_qualities(&mut self, preserved: bool) {
		if !preserved && self.qualities.first() == Some(&MISSING_QUALITY) {
			self.qualities.fill(FILLED_QUALITY);
		}
	}

	fn feature(
		&mut self,
		code: u8,
		index: usize,
		preserved: bool,
		walk: &mut Walk,
	) -> Result<(), Error> {
		let series = &self.header.series;
		match code {
			b'S' => {
				let length = Self::array_into(
					&series.soft_clip,
					&mut self.data,
					&mut self.array,
					&mut self.bases,
					index,
				)?;
				walk.builder.soft_clip(&mut self.out.cigars, length)?;
				walk.read = index + length;
			}
			b'X' => {
				let substitution = u8::try_from(series.substitution.int(&mut self.data)?)
					.map_err(|_| malformed())?;
				self.cover(walk.reference, 1)?;
				let reference = match &self.window {
					Some(window) => window.at(walk.reference)?,
					None => b'N',
				};
				let base = self.header.substitutions.base(reference, substitution);
				*self.bases.get_mut(index).ok_or_else(malformed)? = base;
				walk.builder.add(&mut self.out.cigars, CIGAR_MATCH, 1)?;
				walk.read = index + 1;
				walk.reference = advance(walk.reference, 1)?;
			}
			b'D' => {
				let length = natural(series.deletion_length.int(&mut self.data)?)?;
				walk.builder
					.add(&mut self.out.cigars, CIGAR_DELETION, length)?;
				walk.reference = advance(walk.reference, length)?;
			}
			b'I' => {
				let length = Self::array_into(
					&series.insertion,
					&mut self.data,
					&mut self.array,
					&mut self.bases,
					index,
				)?;
				walk.builder
					.add(&mut self.out.cigars, CIGAR_INSERTION, length)?;
				walk.read = index + length;
			}
			b'i' => {
				let base = series.base.byte(&mut self.data)?;
				*self.bases.get_mut(index).ok_or_else(malformed)? = base;
				walk.builder.add(&mut self.out.cigars, CIGAR_INSERTION, 1)?;
				walk.read = index + 1;
			}
			b'b' => {
				let length = Self::array_into(
					&series.bases,
					&mut self.data,
					&mut self.array,
					&mut self.bases,
					index,
				)?;
				walk.builder
					.add(&mut self.out.cigars, CIGAR_MATCH, length)?;
				walk.read = index + length;
				walk.reference = advance(walk.reference, length)?;
			}
			_ => self.other_feature(code, index, preserved, walk)?,
		}
		Ok(())
	}

	fn other_feature(
		&mut self,
		code: u8,
		index: usize,
		preserved: bool,
		walk: &mut Walk,
	) -> Result<(), Error> {
		let series = &self.header.series;
		match code {
			b'q' => {
				self.missing_qualities(preserved);
				Self::array_into(
					&series.quality_stretch,
					&mut self.data,
					&mut self.array,
					&mut self.qualities,
					index,
				)?;
			}
			b'B' => {
				let base = series.base.byte(&mut self.data)?;
				*self.bases.get_mut(index).ok_or_else(malformed)? = base;
				self.missing_qualities(preserved);
				let quality = series.quality.byte(&mut self.data)?;
				*self.qualities.get_mut(index).ok_or_else(malformed)? = quality;
				walk.builder.add(&mut self.out.cigars, CIGAR_MATCH, 1)?;
				walk.read = index + 1;
				walk.reference = advance(walk.reference, 1)?;
			}
			b'Q' => {
				self.missing_qualities(preserved);
				let quality = series.quality.byte(&mut self.data)?;
				*self.qualities.get_mut(index).ok_or_else(malformed)? = quality;
			}
			b'H' => {
				let length = natural(series.hard_clip.int(&mut self.data)?)?;
				walk.builder
					.add(&mut self.out.cigars, CIGAR_HARD_CLIP, length)?;
			}
			b'P' => {
				let length = natural(series.padding.int(&mut self.data)?)?;
				walk.builder
					.add(&mut self.out.cigars, CIGAR_PADDING, length)?;
			}
			b'N' => {
				let length = natural(series.reference_skip.int(&mut self.data)?)?;
				walk.builder.add(&mut self.out.cigars, CIGAR_SKIP, length)?;
				walk.reference = advance(walk.reference, length)?;
			}
			_ => return Err(malformed()),
		}
		Ok(())
	}

	fn regenerate(&mut self, rec: &mut Rec, flags: TagFlags) -> Result<(), Error> {
		if rec.reference < 0 || rec.cigar.is_empty() {
			return Ok(());
		}
		let start = u64::try_from(rec.position - 1).map_err(|_| malformed())?;
		let span = usize::try_from(reference_span(rec.cigar.of(&self.out.cigars)))
			.map_err(|_| Error::IntegerOverflow)?;
		self.cover(start, span)?;
		let Some(window) = &self.window else {
			return Ok(());
		};
		let distance = regenerate::describe(
			window,
			start,
			rec.cigar.of(&self.out.cigars).iter().copied(),
			&self.bases,
			&mut self.md,
			&mut self.reference_bases,
		)?;
		let tags = &mut self.out.tags;
		if !flags.md {
			tags.extend_from_slice(b"MDZ");
			tags.extend_from_slice(&self.md);
			tags.push(0);
		}
		if !flags.nm {
			tags.extend_from_slice(b"NM");
			if let Ok(small) = u8::try_from(distance) {
				tags.push(b'C');
				tags.push(small);
			} else if let Ok(medium) = u16::try_from(distance) {
				tags.push(b'S');
				tags.extend_from_slice(&medium.to_le_bytes());
			} else {
				tags.push(b'I');
				tags.extend_from_slice(&distance.to_le_bytes());
			}
		}
		rec.tags.end = tags.len();
		Ok(())
	}

	fn store(&mut self, rec: &mut Rec) -> Result<(), Error> {
		let start = self.out.bases.len();
		if rec.cram_flags & CF_NO_SEQUENCE != 0 {
			rec.length = 0;
		} else {
			if start + self.bases.len() > MAX_BASES {
				return Err(Error::TooLarge {
					len: (start + self.bases.len()) as u64,
				});
			}
			self.out.bases.extend_from_slice(&self.bases);
			self.out.qualities.extend_from_slice(&self.qualities);
		}
		rec.bases = Span::since(start, self.out.bases.len());
		rec.qualities = Span::since(start, self.out.qualities.len());
		Ok(())
	}
}

pub(crate) fn advance(position: u64, count: usize) -> Result<u64, Error> {
	position
		.checked_add(u64::try_from(count).map_err(|_| Error::IntegerOverflow)?)
		.ok_or(Error::IntegerOverflow)
}

pub(crate) fn reference_span(cigar: &[u32]) -> u64 {
	cigar
		.iter()
		.filter(|op| matches!(*op & 0xf, 0 | 2 | 3 | 7 | 8))
		.map(|op| u64::from(op >> 4))
		.sum()
}

fn finish_value(kind: u8, tags: &mut Vec<u8>, start: usize) -> Result<(), Error> {
	if matches!(kind, b'Z' | b'H') {
		if tags.len() == start || tags.last() != Some(&0) {
			tags.push(0);
		}
		let value = tags.get(start..).unwrap_or_default();
		return if memchr::memchr(0, value) == Some(value.len() - 1) {
			Ok(())
		} else {
			Err(malformed())
		};
	}
	let value = tags.get(start..).unwrap_or_default();
	let size = match kind {
		b'A' | b'c' | b'C' => 1,
		b's' | b'S' => 2,
		b'i' | b'I' | b'f' => 4,
		b'd' => 8,
		b'B' => return validate_array(value),
		_ => return Err(malformed()),
	};
	if value.len() == size {
		Ok(())
	} else {
		Err(malformed())
	}
}

fn validate_array(value: &[u8]) -> Result<(), Error> {
	let (&[kind, c0, c1, c2, c3], data) = value.split_first_chunk::<5>().ok_or_else(malformed)?;
	let size = match kind {
		b'c' | b'C' => 1,
		b's' | b'S' => 2,
		b'i' | b'I' | b'f' => 4,
		_ => return Err(malformed()),
	};
	let count = u32::from_le_bytes([c0, c1, c2, c3]) as usize;
	if count.checked_mul(size) == Some(data.len()) {
		Ok(())
	} else {
		Err(malformed())
	}
}
