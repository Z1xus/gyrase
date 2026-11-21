use gyrase_genome::{ContigId, Contigs};
use gyrase_reads::{
	CigarKind as BamKind, CigarOp as BamOp, Flags, QualityInput, Record as BamRecord, RecordError,
	RecordRef, Value,
};

use crate::Error;
use crate::aligner::{Read, ReadPair};
use crate::alignment::{
	AlignedBatch, Alignment, CigarKind, CigarOp, End, Flag, Mate, Placement, ReadId, Record, Strand,
};

const PHRED_OFFSET: u8 = 33;
const MIN_MATE_SCORE_QUALITY: u8 = 15;

#[derive(Debug, thiserror::Error)]
pub enum EmitError<E: std::error::Error + 'static> {
	#[error(transparent)]
	Convert(Error),
	#[error(transparent)]
	Sink(E),
}

#[derive(Debug, Clone, Copy)]
pub enum Source<'a> {
	Reads(&'a [Read<'a>]),
	Pairs(&'a [ReadPair<'a>]),
}

impl<'a> Source<'a> {
	fn read(&self, id: ReadId) -> Option<Read<'a>> {
		let index = id.index as usize;
		match (self, id.end) {
			(Self::Reads(reads), End::Single) => reads.get(index).copied(),
			(Self::Pairs(pairs), End::First) => pairs.get(index).map(|pair| pair.first),
			(Self::Pairs(pairs), End::Second) => pairs.get(index).map(|pair| pair.second),
			_ => None,
		}
	}

	fn mate(&self, id: ReadId) -> Option<Read<'a>> {
		let end = match id.end {
			End::First => End::Second,
			End::Second => End::First,
			End::Single => return None,
		};
		self.read(ReadId {
			index: id.index,
			end,
		})
	}
}

#[derive(Debug, Clone, Copy)]
struct Location {
	contig: ContigId,
	position: u32,
}

#[derive(Debug)]
pub struct BamEmitter<'a> {
	contigs: &'a Contigs,
	read_group: Option<&'a [u8]>,
	record: BamRecord,
	operations: Vec<BamOp>,
	bases: Vec<u8>,
	qualities: Vec<u8>,
	text: Vec<u8>,
}

impl<'a> BamEmitter<'a> {
	pub fn new(contigs: &'a Contigs, read_group: Option<&'a [u8]>) -> Self {
		Self {
			contigs,
			read_group,
			record: BamRecord::new(),
			operations: Vec::new(),
			bases: Vec::new(),
			qualities: Vec::new(),
			text: Vec::new(),
		}
	}

	pub fn emit<E: std::error::Error + 'static>(
		&mut self,
		batch: &AlignedBatch,
		source: Source<'_>,
		mut sink: impl FnMut(RecordRef<'_>) -> Result<(), E>,
	) -> Result<(), EmitError<E>> {
		for template in batch
			.records()
			.chunk_by(|left, right| left.read.index == right.read.index)
		{
			for record in template {
				self.build(batch, source, template, record)
					.map_err(EmitError::Convert)?;
				sink(self.record.view()).map_err(EmitError::Sink)?;
			}
		}
		Ok(())
	}

	fn build(
		&mut self,
		batch: &AlignedBatch,
		source: Source<'_>,
		template: &[Record],
		record: &Record,
	) -> Result<(), Error> {
		let read = source
			.read(record.read)
			.ok_or(Error::ForeignBatch(record.read.index))?;
		self.record.clear();
		self.record.set_name(read.name)?;
		let strand = match &record.placement {
			Placement::Mapped(alignment) => {
				self.place(batch, alignment, home(template, record, alignment))?;
				alignment.strand
			}
			Placement::Unmapped { flags, mate } => {
				self.record.set_flags(Flags::new(flags.bits()));
				if let Mate::Mapped(found) = mate {
					let location = Location {
						contig: found.contig,
						position: found.position,
					};
					self.set_location(location)?;
					self.set_mate_location(Some(location))?;
				}
				Strand::Forward
			}
		};
		self.set_sequence(read, strand)?;
		if let Placement::Mapped(alignment) = &record.placement {
			self.push_alignment_tags(batch, alignment)?;
			if !alignment.flags.contains(Flag::Secondary) {
				self.push_supplementary_list(batch, template, record)?;
			}
		}
		if !matches!(record.read.end, End::Single) {
			self.push_mate_tags(batch, source, template, record)?;
		}
		if let Some(group) = self.read_group {
			self.record.push_tag(*b"RG", Value::String(group))?;
		}
		Ok(())
	}

	fn place(
		&mut self,
		batch: &AlignedBatch,
		alignment: &Alignment,
		home: Location,
	) -> Result<(), RecordError> {
		self.record.set_flags(Flags::new(alignment.flags.bits()));
		self.set_location(Location {
			contig: alignment.contig,
			position: alignment.position,
		})?;
		self.record.set_mapq(alignment.mapq);
		self.record.set_template_length(alignment.template_length);
		match alignment.mate {
			Mate::Mapped(found) => self.set_mate_location(Some(Location {
				contig: found.contig,
				position: found.position,
			}))?,
			Mate::Unmapped => self.set_mate_location(Some(home))?,
			Mate::Unpaired => {}
		}
		self.operations.clear();
		for operation in batch.cigar(alignment) {
			self.operations.push(bam_operation(*operation)?);
		}
		self.record.set_cigar(&self.operations)
	}

	fn set_location(&mut self, location: Location) -> Result<(), RecordError> {
		self.record.set_contig(Some(location.contig))?;
		self.record.set_position(Some(location.position))
	}

	fn set_mate_location(&mut self, location: Option<Location>) -> Result<(), RecordError> {
		self.record
			.set_mate_contig(location.map(|found| found.contig))?;
		self.record
			.set_mate_position(location.map(|found| found.position))
	}

	fn set_sequence(&mut self, read: Read<'_>, strand: Strand) -> Result<(), RecordError> {
		self.bases.clear();
		self.qualities.clear();
		match strand {
			Strand::Forward => {
				self.bases.extend_from_slice(read.sequence);
				self.qualities.extend_from_slice(read.qualities);
			}
			Strand::Reverse => {
				self.bases
					.extend(read.sequence.iter().rev().map(|&base| complement(base)));
				self.qualities.extend(read.qualities.iter().rev());
			}
		}
		let qualities = if self.qualities.is_empty() {
			QualityInput::Missing
		} else {
			QualityInput::Ascii(&self.qualities)
		};
		self.record.set_sequence(&self.bases, qualities)
	}

	fn push_alignment_tags(
		&mut self,
		batch: &AlignedBatch,
		alignment: &Alignment,
	) -> Result<(), RecordError> {
		self.record
			.push_integer(*b"NM", i64::from(alignment.edit_distance))?;
		let md = batch.md(alignment);
		if !md.is_empty() {
			self.record.push_tag(*b"MD", Value::String(md))?;
		}
		self.record
			.push_integer(*b"AS", i64::from(alignment.score))?;
		if let Some(suboptimal) = alignment.suboptimal_score {
			self.record.push_integer(*b"XS", i64::from(suboptimal))?;
		}
		Ok(())
	}

	fn push_supplementary_list(
		&mut self,
		batch: &AlignedBatch,
		template: &[Record],
		own: &Record,
	) -> Result<(), RecordError> {
		self.text.clear();
		for other in template {
			let Placement::Mapped(alignment) = &other.placement else {
				continue;
			};
			let chimeric = !alignment.flags.contains(Flag::Secondary);
			if other.read.end != own.read.end || std::ptr::eq(other, own) || !chimeric {
				continue;
			}
			self.write_alignment_text(batch, alignment);
		}
		if self.text.is_empty() {
			return Ok(());
		}
		self.record.push_tag(*b"SA", Value::String(&self.text))
	}

	fn write_alignment_text(&mut self, batch: &AlignedBatch, alignment: &Alignment) {
		let name = self
			.contigs
			.get(alignment.contig)
			.map_or(&b"*"[..], |contig| contig.name.as_bytes());
		self.text.extend_from_slice(name);
		self.text.push(b',');
		push_decimal(&mut self.text, u64::from(alignment.position) + 1);
		self.text.push(b',');
		self.text.push(match alignment.strand {
			Strand::Forward => b'+',
			Strand::Reverse => b'-',
		});
		self.text.push(b',');
		push_cigar(&mut self.text, batch.cigar(alignment));
		self.text.push(b',');
		push_decimal(&mut self.text, u64::from(alignment.mapq));
		self.text.push(b',');
		push_decimal(&mut self.text, u64::from(alignment.edit_distance));
		self.text.push(b';');
	}

	fn push_mate_tags(
		&mut self,
		batch: &AlignedBatch,
		source: Source<'_>,
		template: &[Record],
		own: &Record,
	) -> Result<(), RecordError> {
		let mate_primary = template.iter().find(|candidate| {
			candidate.read.end != own.read.end && candidate.read.index == own.read.index
		});
		if let Some(Placement::Mapped(alignment)) = mate_primary.map(|found| &found.placement) {
			self.text.clear();
			push_cigar(&mut self.text, batch.cigar(alignment));
			self.record.push_tag(*b"MC", Value::String(&self.text))?;
		}
		if let Some(mate) = source.mate(own.read)
			&& !mate.qualities.is_empty()
		{
			self.record
				.push_integer(*b"ms", mate_score(mate.qualities))?;
		}
		Ok(())
	}
}

fn home(template: &[Record], own: &Record, alignment: &Alignment) -> Location {
	let primary = template
		.iter()
		.find(|candidate| candidate.read.end == own.read.end)
		.and_then(|candidate| match &candidate.placement {
			Placement::Mapped(found) => Some(found),
			Placement::Unmapped { .. } => None,
		})
		.unwrap_or(alignment);
	Location {
		contig: primary.contig,
		position: primary.position,
	}
}

fn mate_score(qualities: &[u8]) -> i64 {
	qualities
		.iter()
		.map(|&quality| quality.saturating_sub(PHRED_OFFSET))
		.filter(|&quality| quality >= MIN_MATE_SCORE_QUALITY)
		.map(i64::from)
		.sum()
}

fn bam_operation(operation: CigarOp) -> Result<BamOp, RecordError> {
	let kind = match operation.kind() {
		CigarKind::Match => BamKind::Match,
		CigarKind::Insertion => BamKind::Insertion,
		CigarKind::Deletion => BamKind::Deletion,
		CigarKind::SoftClip => BamKind::SoftClip,
	};
	BamOp::new(operation.length(), kind).ok_or(RecordError::CigarLengthTooLarge)
}

fn complement(base: u8) -> u8 {
	match base.to_ascii_uppercase() {
		b'A' => b'T',
		b'C' => b'G',
		b'G' => b'C',
		b'T' | b'U' => b'A',
		b'R' => b'Y',
		b'Y' => b'R',
		b'K' => b'M',
		b'M' => b'K',
		b'B' => b'V',
		b'V' => b'B',
		b'D' => b'H',
		b'H' => b'D',
		other => other,
	}
}

fn push_cigar(out: &mut Vec<u8>, operations: &[CigarOp]) {
	if operations.is_empty() {
		out.push(b'*');
	}
	for operation in operations {
		push_decimal(out, u64::from(operation.length()));
		out.push(match operation.kind() {
			CigarKind::Match => b'M',
			CigarKind::Insertion => b'I',
			CigarKind::Deletion => b'D',
			CigarKind::SoftClip => b'S',
		});
	}
}

fn push_decimal(out: &mut Vec<u8>, mut value: u64) {
	let mut digits = [0u8; 20];
	let mut start = digits.len();
	loop {
		start -= 1;
		if let Some(slot) = digits.get_mut(start) {
			*slot = b'0' + u8::try_from(value % 10).unwrap_or(0);
		}
		value /= 10;
		if value == 0 {
			break;
		}
	}
	out.extend_from_slice(digits.get(start..).unwrap_or_default());
}
