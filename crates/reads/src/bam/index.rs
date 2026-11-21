use gyrase_genome::Interval;
use gyrase_index::{Builder, Entry, Geometry, Index, Placement};
use gyrase_io::Chunk;

use super::{BamReader, RecordBatch};
use crate::error::Error;
use crate::header::Header;
use crate::record::RecordRef;

const LEVEL_BITS: u32 = 3;
const SLACK: u64 = 256;

pub fn index_bam(reader: &mut BamReader, geometry: Geometry) -> Result<Index, Error> {
	let references = u32::try_from(reader.header().references().len()).unwrap_or(u32::MAX);
	let mut builder = Builder::new(geometry);
	let mut batch = RecordBatch::new();
	while reader.read_batch(&mut batch)? {
		for located in batch.iter() {
			let located = located?;
			builder.push(&entry(
				located.record,
				Chunk {
					begin: located.start,
					end: located.end,
				},
			))?;
		}
	}
	Ok(builder.finish(references)?)
}

pub fn csi_geometry(header: &Header, min_shift: u32) -> Result<Geometry, Error> {
	let longest = longest_reference(header);
	let mut depth = 0;
	let mut size = 1u64 << min_shift;
	while longest > size {
		depth += 1;
		size <<= LEVEL_BITS;
	}
	Ok(Geometry::new(min_shift, depth)?)
}

pub(crate) fn fits_bai(header: &Header) -> bool {
	longest_reference(header) <= Geometry::BAI.max_position()
}

pub(crate) fn entry(record: RecordRef<'_>, offsets: Chunk) -> Entry {
	let placement = match (record.contig(), record.position(), record.reference_end()) {
		(Some(contig), Some(start), Some(end)) => Placement::Placed(Interval {
			contig,
			start: u64::from(start),
			end,
		}),
		_ => Placement::Unplaced,
	};
	Entry {
		placement,
		offsets,
		mapped: record.flags().is_mapped(),
	}
}

fn longest_reference(header: &Header) -> u64 {
	header
		.references()
		.iter()
		.map(|reference| u64::from(reference.length()))
		.max()
		.unwrap_or(0)
		+ SLACK
}
