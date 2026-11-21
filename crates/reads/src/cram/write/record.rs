use super::header::{TEXT_STOP, TagKey};
use super::layout::{int, put_itf8};
use super::series::{Series, Streams};
use crate::cram::Error;
use crate::cram::decode::{Environment, advance};
use crate::cram::records::{
	CF_DETACHED, CF_NO_SEQUENCE, CF_QUALITIES, FLAG_MATE_REVERSE, FLAG_MATE_UNMAPPED,
	FLAG_UNMAPPED, MATE_REVERSE, MATE_UNMAPPED,
};
use crate::cram::reference::Window;
use crate::cram::regenerate;
use crate::{CigarKind, CigarOp, RecordRef, Value};

const SUBSTITUTION_BASES: [u8; 5] = *b"ACGTN";
const NO_READ_GROUP: i32 = -1;
const NO_MATE: i32 = -1;
const PSEUDO_TAG: [u8; 2] = *b"cF";
const LONG_CIGAR_TAG: [u8; 2] = *b"CG";

fn one_based(position: u64) -> Result<i32, Error> {
	i32::try_from(position)
		.ok()
		.and_then(|position| position.checked_add(1))
		.ok_or(Error::IntegerOverflow)
}

fn substitution(expected: u8, base: u8) -> Option<u8> {
	if !SUBSTITUTION_BASES.contains(&expected) {
		return None;
	}
	let rank = SUBSTITUTION_BASES
		.iter()
		.filter(|candidate| **candidate != expected)
		.position(|candidate| *candidate == base)?;
	u8::try_from(rank).ok()
}

#[derive(Debug, Default)]
pub(super) struct Dictionary {
	pub(super) lines: Vec<Vec<TagKey>>,
	pub(super) keys: Vec<TagKey>,
}

impl Dictionary {
	fn line(&mut self, line: &[TagKey]) -> usize {
		for key in line {
			if !self.keys.contains(key) {
				self.keys.push(*key);
			}
		}
		if let Some(index) = self.lines.iter().position(|known| known == line) {
			return index;
		}
		self.lines.push(line.to_vec());
		self.lines.len() - 1
	}
}

#[derive(Debug)]
struct TagEntry {
	key: TagKey,
	start: usize,
	end: usize,
}

#[derive(Debug, Default)]
struct Walk {
	previous: i32,
	count: i32,
}

impl Walk {
	fn feature(&mut self, streams: &mut Streams, read: usize, code: u8) -> Result<(), Error> {
		let position = int(read)?.checked_add(1).ok_or(Error::IntegerOverflow)?;
		streams.byte(Series::FeatureCode, code);
		streams.int(Series::FeaturePosition, position - self.previous);
		self.previous = position;
		self.count += 1;
		Ok(())
	}
}

#[derive(Debug)]
pub(super) struct RecordEncoder<'a> {
	environment: &'a Environment,
	window: Option<&'a Window>,
	pub(super) streams: Streams,
	pub(super) dictionary: Dictionary,
	last_position: i32,
	sequenced: bool,
	bases: Vec<u8>,
	md: Vec<u8>,
	reference_bases: Vec<u8>,
	scratch: Vec<u8>,
	entries: Vec<TagEntry>,
	tag_bytes: Vec<u8>,
	line: Vec<TagKey>,
}

impl<'a> RecordEncoder<'a> {
	pub(super) fn new(
		environment: &'a Environment,
		window: Option<&'a Window>,
		start: i32,
	) -> Self {
		Self {
			environment,
			window,
			streams: Streams::new(),
			dictionary: Dictionary::default(),
			last_position: start,
			sequenced: true,
			bases: Vec::new(),
			md: Vec::new(),
			reference_bases: Vec::new(),
			scratch: Vec::new(),
			entries: Vec::new(),
			tag_bytes: Vec::new(),
			line: Vec::new(),
		}
	}

	pub(super) fn encode(&mut self, record: RecordRef<'_>) -> Result<(), Error> {
		let flags = record.flags().bits();
		let mapped = flags & FLAG_UNMAPPED == 0;
		let length = record.sequence().len();
		self.sequenced = length > 0;
		let read_length = if mapped && !self.sequenced {
			usize::try_from(record.cigar().query_len()).map_err(|_| Error::IntegerOverflow)?
		} else {
			length
		};
		let qualities = record.qualities();
		let has_qualities = !qualities.is_empty() && !qualities.is_missing();
		let mut compression = CF_DETACHED;
		if has_qualities {
			compression |= CF_QUALITIES;
		}
		if !self.sequenced {
			compression |= CF_NO_SEQUENCE;
		}
		self.streams.int(Series::Flags, i32::from(flags));
		self.streams.int(
			Series::CompressionFlags,
			i32::try_from(compression).map_err(|_| Error::IntegerOverflow)?,
		);
		self.streams.int(Series::ReadLength, int(read_length)?);
		let position = match record.position() {
			Some(position) => one_based(u64::from(position))?,
			None if record.contig().is_none() => 0,
			None => return Err(Error::Unsupported("placed record without a position")),
		};
		self.streams
			.int(Series::Position, position - self.last_position);
		self.last_position = position;
		self.bases.clear();
		if self.sequenced {
			record.sequence().decode_into(&mut self.bases);
		} else {
			self.bases.resize(read_length, b'N');
		}
		let distance = if mapped && self.sequenced {
			Some(self.regenerate(record, position, length)?)
		} else {
			None
		};
		let group = self.tags(record, distance)?;
		self.streams.int(Series::ReadGroup, group);
		if record.name().contains(&0) {
			return Err(Error::Unsupported("read name with a NUL byte"));
		}
		self.streams.terminated(Series::Name, record.name());
		self.mate(record, flags)?;
		if mapped {
			self.mapped(record, position)?;
		} else {
			self.streams.bytes(Series::Base, &self.bases);
		}
		if has_qualities {
			self.streams.bytes(Series::Quality, qualities.as_phred());
		}
		Ok(())
	}

	fn regenerate(
		&mut self,
		record: RecordRef<'_>,
		position: i32,
		length: usize,
	) -> Result<u32, Error> {
		let window = self.window.ok_or(Error::ReferenceRequired)?;
		let cigar = record.cigar();
		if cigar.is_empty() || position <= 0 {
			return Err(Error::Unsupported(
				"mapped record without a CIGAR or position",
			));
		}
		if cigar.query_len() != length as u64 {
			return Err(Error::Unsupported("CIGAR length differs from the sequence"));
		}
		let start = u64::try_from(position - 1).map_err(|_| Error::IntegerOverflow)?;
		regenerate::describe(
			&mut |position, target| window.copy(position, target),
			start,
			cigar.ops().map(CigarOp::raw),
			&self.bases,
			&mut self.md,
			&mut self.reference_bases,
		)
	}

	fn read_group(&self, id: &[u8]) -> Option<i32> {
		let index = self
			.environment
			.read_groups
			.iter()
			.position(|group| group == id)
			.filter(|index| Some(*index) != self.environment.unknown_read_group)?;
		i32::try_from(index).ok()
	}

	fn push_tag(&mut self, name: [u8; 2], value: &Value<'_>) -> Result<(), Error> {
		self.scratch.clear();
		value.encode(name, &mut self.scratch);
		let [first, second, kind, rest @ ..] = self.scratch.as_slice() else {
			return Err(Error::Unsupported("empty tag"));
		};
		let text = matches!(kind, b'Z' | b'H');
		let body = if text {
			rest.strip_suffix(&[0]).unwrap_or(rest)
		} else {
			rest
		};
		if text && (body.contains(&0) || body.contains(&TEXT_STOP)) {
			return Err(Error::Unsupported("tag text with a tab or NUL byte"));
		}
		let start = self.tag_bytes.len();
		self.tag_bytes.extend_from_slice(body);
		self.entries.push(TagEntry {
			key: TagKey {
				name: [*first, *second],
				kind: *kind,
			},
			start,
			end: self.tag_bytes.len(),
		});
		Ok(())
	}

	fn tags(&mut self, record: RecordRef<'_>, distance: Option<u32>) -> Result<i32, Error> {
		self.entries.clear();
		self.tag_bytes.clear();
		let mut group = NO_READ_GROUP;
		let (mut md_present, mut nm_present) = (false, false);
		let long_cigar = record.has_long_cigar();
		for field in &record.tags() {
			let field = field?;
			if long_cigar && field.tag == LONG_CIGAR_TAG {
				continue;
			}
			match (field.tag, field.value) {
				(name, Value::String(id)) if name == *b"RG" && group == NO_READ_GROUP => {
					if let Some(index) = self.read_group(id) {
						group = index;
						continue;
					}
				}
				(name, Value::String(text)) if name == *b"MD" && distance.is_some() => {
					md_present = true;
					if text == self.md.as_slice() {
						continue;
					}
				}
				(name, value) if name == *b"NM" && distance.is_some() => {
					nm_present = true;
					if value.as_integer() == distance.map(i64::from) {
						continue;
					}
				}
				_ => {}
			}
			self.push_tag(field.tag, &field.value)?;
		}
		if distance.is_some() {
			let bits = u8::from(!md_present) | (u8::from(!nm_present) << 1);
			if bits != 0 {
				let start = self.tag_bytes.len();
				self.tag_bytes.push(bits);
				self.entries.push(TagEntry {
					key: TagKey {
						name: PSEUDO_TAG,
						kind: b'C',
					},
					start,
					end: start + 1,
				});
			}
		}
		self.line.clear();
		self.line.extend(self.entries.iter().map(|entry| entry.key));
		let line = self.dictionary.line(&self.line);
		self.streams.int(Series::TagLine, int(line)?);
		for entry in &self.entries {
			let value = self
				.tag_bytes
				.get(entry.start..entry.end)
				.unwrap_or_default();
			let stream = self.streams.tag(entry.key.id());
			match entry.key.kind {
				b'Z' | b'H' => {
					stream.extend_from_slice(value);
					stream.extend_from_slice(&[0, TEXT_STOP]);
				}
				b'B' => {
					put_itf8(stream, int(value.len())?);
					stream.extend_from_slice(value);
				}
				_ => stream.extend_from_slice(value),
			}
		}
		Ok(group)
	}

	fn mate(&mut self, record: RecordRef<'_>, flags: u16) -> Result<(), Error> {
		let mut mate_flags = 0u8;
		if flags & FLAG_MATE_REVERSE != 0 {
			mate_flags |= MATE_REVERSE;
		}
		if flags & FLAG_MATE_UNMAPPED != 0 {
			mate_flags |= MATE_UNMAPPED;
		}
		self.streams.int(Series::MateFlags, i32::from(mate_flags));
		let reference = match record.mate_contig() {
			Some(contig) => i32::try_from(contig.get()).map_err(|_| Error::IntegerOverflow)?,
			None => NO_MATE,
		};
		self.streams.int(Series::MateReference, reference);
		let position = match record.mate_position() {
			Some(position) => one_based(u64::from(position))?,
			None => 0,
		};
		self.streams.int(Series::MatePosition, position);
		self.streams
			.int(Series::TemplateLength, record.template_length());
		Ok(())
	}

	fn mapped(&mut self, record: RecordRef<'_>, position: i32) -> Result<(), Error> {
		let mut walk = Walk::default();
		let mut read = 0usize;
		let mut reference = u64::try_from(position - 1).map_err(|_| Error::IntegerOverflow)?;
		for op in record.cigar().ops() {
			let length = op.len() as usize;
			match op.kind() {
				CigarKind::Match | CigarKind::SeqMatch | CigarKind::SeqMismatch => {
					self.matches(&mut walk, read, reference, length)?;
					read += length;
					reference = advance(reference, length)?;
				}
				CigarKind::Insertion | CigarKind::SoftClip => {
					let (code, series) = if op.kind() == CigarKind::Insertion {
						(b'I', Series::Insertion)
					} else {
						(b'S', Series::SoftClip)
					};
					walk.feature(&mut self.streams, read, code)?;
					let bases = self
						.bases
						.get(read..read + length)
						.ok_or(Error::Unsupported("CIGAR length differs from the sequence"))?;
					self.streams.terminated(series, bases);
					read += length;
				}
				CigarKind::Deletion | CigarKind::Skip => {
					let (code, series) = if op.kind() == CigarKind::Deletion {
						(b'D', Series::DeletionLength)
					} else {
						(b'N', Series::ReferenceSkip)
					};
					walk.feature(&mut self.streams, read, code)?;
					self.streams.int(series, int(length)?);
					reference = advance(reference, length)?;
				}
				CigarKind::HardClip | CigarKind::Padding => {
					let (code, series) = if op.kind() == CigarKind::HardClip {
						(b'H', Series::HardClip)
					} else {
						(b'P', Series::Padding)
					};
					walk.feature(&mut self.streams, read, code)?;
					self.streams.int(series, int(length)?);
				}
			}
		}
		self.streams.int(Series::FeatureCount, walk.count);
		self.streams
			.int(Series::MappingQuality, i32::from(record.mapq()));
		Ok(())
	}

	fn matches(
		&mut self,
		walk: &mut Walk,
		read: usize,
		reference: u64,
		length: usize,
	) -> Result<(), Error> {
		if !self.sequenced {
			return Ok(());
		}
		let window = self.window.ok_or(Error::ReferenceRequired)?;
		self.reference_bases.resize(length, 0);
		window.copy(reference, &mut self.reference_bases)?;
		let aligned = self
			.bases
			.get(read..read + length)
			.ok_or(Error::Unsupported("CIGAR length differs from the sequence"))?;
		if aligned == self.reference_bases.as_slice() {
			return Ok(());
		}
		for (offset, (base, expected)) in aligned.iter().zip(&self.reference_bases).enumerate() {
			if base == expected {
				continue;
			}
			if let Some(code) = substitution(*expected, *base) {
				walk.feature(&mut self.streams, read + offset, b'X')?;
				self.streams.int(Series::Substitution, i32::from(code));
			} else {
				walk.feature(&mut self.streams, read + offset, b'b')?;
				self.streams.terminated(Series::Bases, &[*base]);
			}
		}
		Ok(())
	}
}
