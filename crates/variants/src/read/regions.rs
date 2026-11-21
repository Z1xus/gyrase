use std::fs::File;
use std::path::Path;

use gyrase_genome::Interval;
use gyrase_index::{DEFAULT_MAX_GAP, Index};
use gyrase_io::Workers;

use super::any::Reader;
use super::bcf::BcfReader;
use super::chunks::ChunkReader;
use super::vcf::VcfReader;
use crate::Error;
use crate::region::Regions;

pub fn open_regions(
	path: &Path,
	index: &Index,
	texts: &[&str],
	workers: &Workers,
) -> Result<(Reader<ChunkReader>, Regions), Error> {
	let opened = Reader::open(Box::new(File::open(path)?), workers)?;
	let header = opened.header().clone();
	let is_bcf = matches!(opened, Reader::Bcf(_));
	drop(opened);
	let contigs = header.contigs()?;
	let intervals = texts
		.iter()
		.map(|text| Interval::parse(text, &contigs))
		.collect::<Result<Vec<_>, _>>()?;
	let tabix = index.tabix().ok();
	let indexed: Vec<Interval> = intervals
		.iter()
		.filter_map(|interval| {
			let contig = match &tabix {
				Some(names) => {
					let entry = header.contig(interval.contig.get())?;
					names.contig(&entry.name)?
				}
				None => interval.contig,
			};
			Some(Interval {
				contig,
				..*interval
			})
		})
		.collect();
	let plan = index.plan(&indexed, DEFAULT_MAX_GAP);
	let chunks = ChunkReader::new(File::open(path)?, plan.chunks(), workers);
	let reader = if is_bcf {
		Reader::Bcf(BcfReader::with_header(chunks, header))
	} else {
		Reader::Vcf(VcfReader::with_header(chunks, header))
	};
	Ok((reader, Regions::new(&intervals)))
}
