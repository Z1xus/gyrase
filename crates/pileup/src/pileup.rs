use std::mem;
use std::num::NonZeroU32;

use gyrase_genome::{Base, Code4, Interval};
use gyrase_reads::{CigarKind, CigarOp, RecordRef};

use crate::column::merge_mates;
use crate::indel::{left_align_deletion, left_align_insertion};
use crate::{Change, Column, Filter, Indel, Observation, Strand};

const MISSING_QUALITY: u8 = 0xff;
const ASSUMED_QUALITY: u8 = 30;
const FNV_OFFSET: u32 = 0x811c_9dc5;
const FNV_PRIME: u32 = 0x0100_0193;
const FREE_MISMATCHES: usize = 1;
const COUNTED_MISMATCH_QUALITY: u8 = 20;
const MAX_MAPPING_QUALITY: u8 = 60;

#[derive(Debug, Clone, Copy)]
pub struct Window<'a> {
	pub start: u64,
	pub bases: &'a [u8],
}

#[derive(Debug, Default)]
pub struct Scratch {
	offsets: Vec<usize>,
	cursors: Vec<usize>,
	coverage: Vec<i64>,
	observations: Vec<Observation>,
	indels: Vec<Indel>,
	inserted: Vec<u8>,
	read_bases: Vec<u8>,
	read_qualities: Vec<u8>,
	reads: Vec<Stored>,
	segments: Vec<Segment>,
	sequence: Vec<u8>,
	longest: u64,
}

#[derive(Debug, Clone, Copy)]
struct Segment {
	reference: u64,
	query: usize,
	len: usize,
}

#[derive(Debug, Clone, Copy)]
struct Stored {
	start: u64,
	end: u64,
	clipped_start: u64,
	offset: usize,
	len: usize,
	segments: (usize, usize),
	mapping_quality: u8,
	fragment: Option<NonZeroU32>,
}

#[derive(Debug, Clone, Copy)]
pub struct Read<'a> {
	pub start: u64,
	pub end: u64,
	pub clipped_start: u64,
	pub bases: &'a [u8],
	pub qualities: &'a [u8],
	pub mapping_quality: u8,
	pub fragment: Option<NonZeroU32>,
	segments: &'a [Segment],
}

#[derive(Debug, Clone, Copy)]
pub struct Pileup<'a> {
	pub interval: Interval,
	pub window: Window<'a>,
	pub filter: Filter,
	offsets: &'a [usize],
	observations: &'a [Observation],
	indels: &'a [Indel],
	inserted: &'a [u8],
	reads: &'a [Stored],
	segments: &'a [Segment],
	sequence: &'a [u8],
	longest: u64,
}

struct Shard<'a> {
	interval: Interval,
	window: Window<'a>,
	filter: Filter,
}

struct Alignment<'a> {
	start: u64,
	end: u64,
	bases: &'a [u8],
	qualities: &'a [u8],
	mapping_quality: u8,
	strand: Strand,
	fragment: Option<NonZeroU32>,
}

impl<'a> Read<'a> {
	pub fn slice(&self, from: u64, to: u64) -> Option<(u64, &'a [u8], &'a [u8])> {
		let first = self
			.segments
			.iter()
			.find(|segment| segment.reference + segment.len as u64 > from)?;
		let last = self
			.segments
			.iter()
			.rev()
			.find(|segment| segment.reference < to)?;
		let start = from.max(first.reference);
		let begin = first.query + usize::try_from(start - first.reference).ok()?;
		let end = last.query + usize::try_from(to - last.reference).ok()?.min(last.len);
		let range = begin..end.min(self.bases.len());
		(!range.is_empty()).then(|| (start, &self.bases[range.clone()], &self.qualities[range]))
	}
}

impl Scratch {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn pileup<'a, 'r, E>(
		&'a mut self,
		interval: Interval,
		window: Window<'a>,
		filter: Filter,
		records: impl Iterator<Item = Result<RecordRef<'r>, E>>,
	) -> Result<Pileup<'a>, E> {
		let mut records = records.collect::<Result<Vec<_>, E>>()?;
		records.sort_by_cached_key(|record| {
			(
				record.position(),
				fragment_id(record.name()),
				record.flags().bits(),
			)
		});
		for tied in records.chunk_by_mut(|first, second| {
			first.position() == second.position()
				&& first.flags() == second.flags()
				&& fragment_id(first.name()) == fragment_id(second.name())
		}) {
			tied.sort_by(|first, second| {
				first
					.name()
					.cmp(second.name())
					.then_with(|| {
						first
							.cigar()
							.ops()
							.map(CigarOp::raw)
							.cmp(second.cigar().ops().map(CigarOp::raw))
					})
					.then_with(|| first.mapq().cmp(&second.mapq()))
					.then_with(|| first.mate_position().cmp(&second.mate_position()))
					.then_with(|| {
						first
							.sequence()
							.codes()
							.map(Code4::code)
							.cmp(second.sequence().codes().map(Code4::code))
					})
					.then_with(|| {
						first
							.qualities()
							.as_phred()
							.cmp(second.qualities().as_phred())
					})
			});
		}
		let span = usize::try_from(interval.len()).unwrap_or(0);
		self.indels.clear();
		self.inserted.clear();
		self.reads.clear();
		self.segments.clear();
		self.sequence.clear();
		self.longest = 0;
		let shard = Shard {
			interval,
			window,
			filter,
		};
		self.coverage.clear();
		self.coverage.resize(span + 1, 0);
		for record in &records {
			self.count(record, &shard);
		}
		self.offsets.clear();
		self.offsets.push(0);
		self.cursors.clear();
		let mut depth = 0;
		let mut total = 0;
		for &change in &self.coverage[..span] {
			depth += change;
			self.cursors.push(total);
			total += usize::try_from(depth)
				.unwrap_or(0)
				.min(filter.max_depth as usize);
			self.offsets.push(total);
		}
		self.observations.clear();
		self.observations.resize(total, placeholder());
		for record in &records {
			self.place(record, &shard);
		}
		self.merge_mates(span);
		self.indels.sort_by_key(|indel| indel.anchor);
		Ok(Pileup {
			interval,
			window,
			filter,
			offsets: &self.offsets,
			observations: &self.observations,
			indels: &self.indels,
			inserted: &self.inserted,
			reads: &self.reads,
			segments: &self.segments,
			sequence: &self.sequence,
			longest: self.longest,
		})
	}

	fn count(&mut self, record: &RecordRef<'_>, shard: &Shard<'_>) {
		let Some((start, _)) = overlap(record, shard) else {
			return;
		};
		let interval = shard.interval;
		let mut reference_position = start;
		for operation in record.cigar().ops() {
			let len = u64::from(operation.len());
			if matches!(
				operation.kind(),
				CigarKind::Match | CigarKind::SeqMatch | CigarKind::SeqMismatch
			) {
				let first = reference_position.clamp(interval.start, interval.end) - interval.start;
				let last =
					(reference_position + len).clamp(interval.start, interval.end) - interval.start;
				if let (Ok(first), Ok(last)) = (usize::try_from(first), usize::try_from(last)) {
					self.coverage[first] += 1;
					self.coverage[last] -= 1;
				}
			}
			if operation.kind().consumes_reference() {
				reference_position += len;
			}
		}
	}

	fn place(&mut self, record: &RecordRef<'_>, shard: &Shard<'_>) {
		let Some((start, end)) = overlap(record, shard) else {
			return;
		};
		let flags = record.flags();
		let mut bases = mem::take(&mut self.read_bases);
		bases.clear();
		record.sequence().decode_into(&mut bases);
		let mut qualities = mem::take(&mut self.read_qualities);
		qualities.clear();
		let recorded = record.qualities();
		if !recorded.is_missing() {
			qualities.extend(
				recorded
					.as_phred()
					.iter()
					.take(bases.len())
					.map(|&quality| match quality {
						MISSING_QUALITY => ASSUMED_QUALITY,
						quality => quality,
					}),
			);
		}
		qualities.resize(bases.len(), ASSUMED_QUALITY);
		let read = Alignment {
			start,
			end,
			bases: &bases,
			qualities: &qualities,
			mapping_quality: burdened_mapping_quality(record, start, &bases, &qualities, shard),
			strand: if flags.is_reverse() {
				Strand::Reverse
			} else {
				Strand::Forward
			},
			fragment: overlaps_mate(record, start, end).then(|| fragment_id(record.name())),
		};
		self.store(&read, record);
		let mut reference_position = start;
		let mut query = 0;
		for operation in record.cigar().ops() {
			let len = operation.len() as usize;
			match operation.kind() {
				CigarKind::Match | CigarKind::SeqMatch | CigarKind::SeqMismatch => {
					self.observe(&read, reference_position, query, len, shard);
				}
				CigarKind::Insertion if query > 0 && reference_position > 0 => {
					let inserted = bases.get(query..query + len).unwrap_or_default();
					self.insertion(&read, reference_position - 1, inserted, shard);
				}
				CigarKind::Deletion if query > 0 && reference_position > 0 => {
					self.deletion(&read, reference_position - 1, operation.len(), shard);
				}
				_ => {}
			}
			if operation.kind().consumes_reference() {
				reference_position += u64::from(operation.len());
			}
			if operation.kind().consumes_query() {
				query += len;
			}
		}
		self.read_bases = bases;
		self.read_qualities = qualities;
	}

	fn store(&mut self, read: &Alignment<'_>, record: &RecordRef<'_>) {
		let first = self.segments.len();
		let mut reference = read.start;
		let mut query = 0;
		let mut leading_clip = 0;
		for operation in record.cigar().ops() {
			let len = operation.len() as usize;
			match operation.kind() {
				CigarKind::SoftClip if query == 0 => {
					leading_clip = operation.len().into();
					self.segments.push(Segment {
						reference: reference.saturating_sub(leading_clip),
						query,
						len,
					});
				}
				CigarKind::Match
				| CigarKind::SeqMatch
				| CigarKind::SeqMismatch
				| CigarKind::SoftClip => {
					self.segments.push(Segment {
						reference,
						query,
						len,
					});
				}
				_ => {}
			}
			if operation.kind().consumes_reference() {
				reference += u64::from(operation.len());
			}
			if operation.kind().consumes_query() {
				query += len;
			}
		}
		let offset = self.sequence.len();
		self.sequence.extend_from_slice(read.bases);
		self.sequence.extend_from_slice(read.qualities);
		self.longest = self.longest.max(read.end - read.start + leading_clip);
		self.reads.push(Stored {
			start: read.start,
			end: read.end,
			clipped_start: read.start.saturating_sub(leading_clip),
			offset,
			len: read.bases.len(),
			segments: (first, self.segments.len()),
			mapping_quality: read.mapping_quality,
			fragment: read.fragment,
		});
	}

	fn observe(
		&mut self,
		read: &Alignment<'_>,
		reference_position: u64,
		query: usize,
		len: usize,
		shard: &Shard<'_>,
	) {
		let interval = shard.interval;
		let first = reference_position.max(interval.start);
		let last = (reference_position + len as u64).min(interval.end);
		let (Ok(skip), Ok(offset), Ok(count)) = (
			usize::try_from(first - reference_position),
			usize::try_from(first - interval.start),
			usize::try_from(last.saturating_sub(first)),
		) else {
			return;
		};
		let query = query + skip;
		let (Some(bases), Some(qualities), Some(cursors), Some(ends)) = (
			read.bases.get(query..query + count),
			read.qualities.get(query..query + count),
			self.cursors.get_mut(offset..offset + count),
			self.offsets.get(offset + 1..offset + 1 + count),
		) else {
			return;
		};
		let columns = bases.iter().zip(qualities).zip(cursors).zip(ends);
		for (position, (((&byte, &quality), cursor), &end)) in (first..last).zip(columns) {
			let base = Base::from_ascii(byte);
			if base == Base::N || quality < shard.filter.min_base_quality || *cursor >= end {
				continue;
			}
			let flank = (position - read.start).min(read.end - 1 - position);
			self.observations[*cursor] = Observation {
				base,
				quality,
				mapping_quality: read.mapping_quality,
				strand: read.strand,
				flank: u8::try_from(flank).unwrap_or(u8::MAX),
				fragment: read.fragment,
			};
			*cursor += 1;
		}
	}

	fn insertion(&mut self, read: &Alignment<'_>, anchor: u64, bases: &[u8], shard: &Shard<'_>) {
		let Some(offset) = window_offset(anchor, shard.window) else {
			return;
		};
		let (Ok(start), Ok(len)) = (
			u32::try_from(self.inserted.len()),
			u32::try_from(bases.len()),
		) else {
			return;
		};
		let begin = self.inserted.len();
		self.inserted.extend_from_slice(bases);
		let aligned = shard.window.start
			+ left_align_insertion(shard.window.bases, offset, &mut self.inserted[begin..]) as u64;
		if shard.interval.contains(shard.interval.contig, aligned) {
			self.push_indel(read, aligned, Change::Insertion { start, len });
		} else {
			self.inserted.truncate(begin);
		}
	}

	fn deletion(&mut self, read: &Alignment<'_>, anchor: u64, len: u32, shard: &Shard<'_>) {
		let Some(offset) = window_offset(anchor, shard.window) else {
			return;
		};
		let aligned = shard.window.start
			+ left_align_deletion(shard.window.bases, offset, len as usize) as u64;
		if shard.interval.contains(shard.interval.contig, aligned) {
			self.push_indel(read, aligned, Change::Deletion { len });
		}
	}

	fn push_indel(&mut self, read: &Alignment<'_>, anchor: u64, change: Change) {
		let flank = anchor
			.saturating_sub(read.start)
			.min(read.end.saturating_sub(anchor + 1));
		self.indels.push(Indel {
			anchor,
			change,
			mapping_quality: read.mapping_quality,
			strand: read.strand,
			flank: u8::try_from(flank).unwrap_or(u8::MAX),
			fragment: read.fragment,
		});
	}

	fn merge_mates(&mut self, span: usize) {
		let mut write = 0;
		for column in 0..span {
			let (begin, end) = (self.offsets[column], self.cursors[column]);
			let kept = merge_mates(&mut self.observations[begin..end]);
			self.observations.copy_within(begin..begin + kept, write);
			self.offsets[column] = write;
			write += kept;
		}
		self.offsets[span] = write;
		self.observations.truncate(write);
	}
}

impl<'a> Pileup<'a> {
	pub fn columns(&self) -> impl Iterator<Item = Column<'a>> + '_ {
		self.offsets
			.windows(2)
			.zip(self.interval.start..)
			.map(|(range, position)| Column {
				position,
				reference: self.reference_base(position),
				observations: &self.observations[range[0]..range[1]],
			})
	}

	pub fn column(&self, position: u64) -> Option<Column<'a>> {
		let offset = usize::try_from(position.checked_sub(self.interval.start)?).ok()?;
		let range = self.offsets.get(offset..offset + 2)?;
		Some(Column {
			position,
			reference: self.reference_base(position),
			observations: self.observations.get(range[0]..range[1])?,
		})
	}

	pub fn indels(&self) -> &'a [Indel] {
		self.indels
	}

	pub fn indels_at(&self, anchor: u64) -> &'a [Indel] {
		let start = self.indels.partition_point(|indel| indel.anchor < anchor);
		let end = self.indels.partition_point(|indel| indel.anchor <= anchor);
		&self.indels[start..end]
	}

	pub fn reads_covering(&self, position: u64) -> impl Iterator<Item = Read<'a>> + '_ {
		let end = self.reads.partition_point(|read| read.start <= position);
		let begin = self
			.reads
			.partition_point(|read| read.start + self.longest < position);
		self.reads[begin.min(end)..end]
			.iter()
			.filter(move |read| read.end > position)
			.map(|read| Read {
				start: read.start,
				end: read.end,
				clipped_start: read.clipped_start,
				bases: &self.sequence[read.offset..read.offset + read.len],
				qualities: &self.sequence[read.offset + read.len..read.offset + 2 * read.len],
				mapping_quality: read.mapping_quality,
				fragment: read.fragment,
				segments: &self.segments[read.segments.0..read.segments.1],
			})
	}

	pub fn inserted(&self, start: u32, len: u32) -> &'a [u8] {
		let start = start as usize;
		self.inserted
			.get(start..start + len as usize)
			.unwrap_or_default()
	}

	pub fn reference_base(&self, position: u64) -> Base {
		window_offset(position, self.window)
			.and_then(|offset| self.window.bases.get(offset))
			.map_or(Base::N, |&byte| Base::from_ascii(byte))
	}
}

fn overlap(record: &RecordRef<'_>, shard: &Shard<'_>) -> Option<(u64, u64)> {
	let interval = shard.interval;
	if record.contig() != Some(interval.contig)
		|| !shard.filter.accepts(record.flags(), record.mapq())
	{
		return None;
	}
	let start = u64::from(record.position()?);
	let end = record.reference_end()?;
	(end > interval.start && start < interval.end).then_some((start, end))
}

fn burdened_mapping_quality(
	record: &RecordRef<'_>,
	start: u64,
	bases: &[u8],
	qualities: &[u8],
	shard: &Shard<'_>,
) -> u8 {
	let excess =
		mismatches(record, start, bases, qualities, shard.window).saturating_sub(FREE_MISMATCHES);
	let penalty = u8::try_from(excess)
		.unwrap_or(u8::MAX)
		.saturating_mul(shard.filter.mismatch_penalty);
	record
		.mapq()
		.min(MAX_MAPPING_QUALITY.saturating_sub(penalty))
}

fn mismatches(
	record: &RecordRef<'_>,
	start: u64,
	bases: &[u8],
	qualities: &[u8],
	window: Window<'_>,
) -> usize {
	let mut reference_position = start;
	let mut query = 0;
	let mut count = 0;
	for operation in record.cigar().ops() {
		let len = operation.len() as usize;
		if matches!(
			operation.kind(),
			CigarKind::Match | CigarKind::SeqMatch | CigarKind::SeqMismatch
		) {
			let expected = window_offset(reference_position, window)
				.and_then(|offset| window.bases.get(offset..))
				.unwrap_or_default();
			count += expected
				.iter()
				.zip(bases.get(query..query + len).unwrap_or_default())
				.zip(qualities.get(query..query + len).unwrap_or_default())
				.filter(|&((&expected, &base), &quality)| {
					base != expected
						&& base != b'N' && expected != b'N'
						&& quality >= COUNTED_MISMATCH_QUALITY
				})
				.count();
		}
		if operation.kind().consumes_reference() {
			reference_position += u64::from(operation.len());
		}
		if operation.kind().consumes_query() {
			query += len;
		}
	}
	count
}

fn window_offset(position: u64, window: Window<'_>) -> Option<usize> {
	let offset = usize::try_from(position.checked_sub(window.start)?).ok()?;
	(offset < window.bases.len()).then_some(offset)
}

fn overlaps_mate(record: &RecordRef<'_>, start: u64, end: u64) -> bool {
	let flags = record.flags();
	let Some(mate_start) = record.mate_position().map(u64::from) else {
		return false;
	};
	flags.is_paired()
		&& !flags.is_mate_unmapped()
		&& record.mate_contig() == record.contig()
		&& mate_start < end
		&& start < mate_start + (end - start)
}

fn fragment_id(name: &[u8]) -> NonZeroU32 {
	let hash = name.iter().fold(FNV_OFFSET, |hash, &byte| {
		(hash ^ u32::from(byte)).wrapping_mul(FNV_PRIME)
	});
	NonZeroU32::new(hash).unwrap_or(NonZeroU32::MIN)
}

const fn placeholder() -> Observation {
	Observation {
		base: Base::N,
		quality: 0,
		mapping_quality: 0,
		strand: Strand::Forward,
		flank: 0,
		fragment: None,
	}
}
