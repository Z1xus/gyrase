use gyrase_genome::Interval;
use gyrase_index::{Builder, Entry, Format, Geometry, Index, Placement};
use gyrase_io::Chunk;

use super::{BamReader, RecordBatch};
use crate::error::Error;
use crate::header::Header;
use crate::record::RecordRef;

const SLACK: u64 = 256;
const CSI_MIN_SHIFT: u32 = 14;

pub fn index_bam<E: From<Error>>(
	reader: &mut BamReader,
	geometry: Geometry,
	mut progress: impl FnMut(u64) -> Result<(), E>,
) -> Result<Index, E> {
	let references = u32::try_from(reader.header().references().len()).unwrap_or(u32::MAX);
	let mut builder = Builder::new(geometry);
	let mut batch = RecordBatch::new();
	let mut position = 0;
	while reader.read_batch(&mut batch)? {
		for located in batch.iter() {
			let located = located.map_err(Error::from)?;
			position = located.end.block();
			builder
				.push(&entry(
					located.record,
					Chunk {
						begin: located.start,
						end: located.end,
					},
				))
				.map_err(Error::from)?;
		}
		progress(position)?;
	}
	Ok(builder.finish(references).map_err(Error::from)?)
}

pub fn index_layout(header: &Header) -> Result<(Format, Geometry), Error> {
	if longest_reference(header) <= Geometry::BAI.max_position() {
		Ok((Format::Bai, Geometry::BAI))
	} else {
		Ok((Format::Csi, csi_geometry(header, CSI_MIN_SHIFT)?))
	}
}

pub fn csi_geometry(header: &Header, min_shift: u32) -> Result<Geometry, Error> {
	Ok(Geometry::covering(min_shift, longest_reference(header))?)
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
