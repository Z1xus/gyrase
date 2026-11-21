use gyrase_genome::{ContigId, Interval};
use gyrase_index::{Builder, Chunk, Entry, Geometry, Index, Placement, TabixHeader};
use gyrase_io::{BlockPosition, Finished};

use crate::Error;
use crate::header::Header;
use crate::write::{Placed, Span};

const MIN_SHIFT: u32 = 14;
const DEPTH: u32 = 5;
const TABIX_VCF_PRESET: u32 = 2;
const TABIX_COMMENT: u32 = b'#' as u32;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexKind {
	Vcf,
	Bcf,
}

#[derive(Debug, Default)]
pub struct IndexCollector {
	marks: Vec<(Span, Placed<BlockPosition>)>,
}

impl IndexCollector {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn push(&mut self, span: Span, placed: Placed<BlockPosition>) {
		self.marks.push((span, placed));
	}

	pub fn build<W>(
		&self,
		finished: &Finished<W>,
		header: &Header,
		kind: IndexKind,
	) -> Result<Index, Error> {
		let mut builder = Builder::new(Geometry::new(MIN_SHIFT, DEPTH)?);
		let mut seen: Vec<u32> = Vec::new();
		for (span, placed) in &self.marks {
			let offsets = Chunk {
				begin: finished
					.resolve(placed.begin)
					.ok_or(Error::Chunk("record start is not in the written blocks"))?,
				end: finished
					.resolve(placed.end)
					.ok_or(Error::Chunk("record end is not in the written blocks"))?,
			};
			let slot = match kind {
				IndexKind::Bcf => span.contig,
				IndexKind::Vcf => tabix_slot(&mut seen, span.contig)?,
			};
			builder.push(&Entry {
				placement: Placement::Placed(Interval {
					contig: ContigId::new(slot),
					start: u64::from(span.start),
					end: u64::from(span.end),
				}),
				offsets,
				mapped: true,
			})?;
		}
		let references = match kind {
			IndexKind::Bcf => header.contig_count(),
			IndexKind::Vcf => seen.len(),
		};
		let index =
			builder.finish(u32::try_from(references).map_err(|_| Error::Limits("contigs"))?)?;
		match kind {
			IndexKind::Bcf => Ok(index),
			IndexKind::Vcf => {
				let names = seen
					.iter()
					.map(|&contig| {
						header
							.contig(contig)
							.map(|entry| entry.name.clone())
							.ok_or(Error::Header("contig is not in the header"))
					})
					.collect::<Result<Vec<_>, _>>()?;
				Ok(index.with_tabix(&TabixHeader {
					preset: TABIX_VCF_PRESET,
					sequence_column: 1,
					begin_column: 2,
					end_column: 0,
					comment: TABIX_COMMENT,
					skip_lines: 0,
					names,
				}))
			}
		}
	}
}

fn tabix_slot(seen: &mut Vec<u32>, contig: u32) -> Result<u32, Error> {
	let position = seen
		.iter()
		.position(|&known| known == contig)
		.unwrap_or_else(|| {
			seen.push(contig);
			seen.len() - 1
		});
	u32::try_from(position).map_err(|_| Error::Limits("contigs"))
}
