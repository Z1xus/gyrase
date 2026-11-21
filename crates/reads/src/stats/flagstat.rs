use std::fmt::Write;

use gyrase_io::Workers;
use rayon::prelude::*;

use crate::bam::{BamReader, RecordBatch};
use crate::error::{Error, RecordError};
use crate::record::RecordRef;

const MIN_DIFFERENT_CONTIG_MAPQ: u8 = 5;

type Row<'a> = (&'a [u64; 2], &'a str, Option<&'a [u64; 2]>);

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Flagstat {
	reads: [u64; 2],
	primary: [u64; 2],
	secondary: [u64; 2],
	supplementary: [u64; 2],
	duplicates: [u64; 2],
	primary_duplicates: [u64; 2],
	mapped: [u64; 2],
	primary_mapped: [u64; 2],
	paired: [u64; 2],
	read1: [u64; 2],
	read2: [u64; 2],
	properly_paired: [u64; 2],
	both_mapped: [u64; 2],
	singletons: [u64; 2],
	different_contig: [u64; 2],
	different_contig_confident: [u64; 2],
}

impl Flagstat {
	pub fn add(&mut self, record: RecordRef<'_>) {
		let flags = record.flags();
		let class = usize::from(flags.is_qc_fail());
		bump(&mut self.reads, class);
		if flags.is_secondary() {
			bump(&mut self.secondary, class);
		} else if flags.is_supplementary() {
			bump(&mut self.supplementary, class);
		} else {
			bump(&mut self.primary, class);
			if flags.is_paired() {
				bump(&mut self.paired, class);
				if flags.is_proper_pair() && flags.is_mapped() {
					bump(&mut self.properly_paired, class);
				}
				if flags.is_first_in_pair() {
					bump(&mut self.read1, class);
				}
				if flags.is_second_in_pair() {
					bump(&mut self.read2, class);
				}
				if flags.is_mate_unmapped() && flags.is_mapped() {
					bump(&mut self.singletons, class);
				}
				if flags.is_mapped() && !flags.is_mate_unmapped() {
					bump(&mut self.both_mapped, class);
					if record.mate_contig() != record.contig() {
						bump(&mut self.different_contig, class);
						if record.mapq() >= MIN_DIFFERENT_CONTIG_MAPQ {
							bump(&mut self.different_contig_confident, class);
						}
					}
				}
			}
			if flags.is_mapped() {
				bump(&mut self.primary_mapped, class);
			}
			if flags.is_duplicate() {
				bump(&mut self.primary_duplicates, class);
			}
		}
		if flags.is_mapped() {
			bump(&mut self.mapped, class);
		}
		if flags.is_duplicate() {
			bump(&mut self.duplicates, class);
		}
	}

	pub fn merge(&mut self, other: &Self) {
		for (mine, theirs) in self.counters_mut().into_iter().zip(other.counters()) {
			mine[0] += theirs[0];
			mine[1] += theirs[1];
		}
	}

	pub fn total(&self) -> u64 {
		self.reads[0] + self.reads[1]
	}

	pub fn render(&self) -> String {
		let mut out = String::new();
		let lines: [Row<'_>; 16] = [
			(
				&self.reads,
				"in total (QC-passed reads + QC-failed reads)",
				None,
			),
			(&self.primary, "primary", None),
			(&self.secondary, "secondary", None),
			(&self.supplementary, "supplementary", None),
			(&self.duplicates, "duplicates", None),
			(&self.primary_duplicates, "primary duplicates", None),
			(&self.mapped, "mapped", Some(&self.reads)),
			(&self.primary_mapped, "primary mapped", Some(&self.primary)),
			(&self.paired, "paired in sequencing", None),
			(&self.read1, "read1", None),
			(&self.read2, "read2", None),
			(&self.properly_paired, "properly paired", Some(&self.paired)),
			(&self.both_mapped, "with itself and mate mapped", None),
			(&self.singletons, "singletons", Some(&self.paired)),
			(
				&self.different_contig,
				"with mate mapped to a different chr",
				None,
			),
			(
				&self.different_contig_confident,
				"with mate mapped to a different chr (mapQ>=5)",
				None,
			),
		];
		for (counts, label, total) in lines {
			let _ = write!(out, "{} + {} {label}", counts[0], counts[1]);
			if let Some(total) = total {
				let _ = write!(
					out,
					" ({} : {})",
					percent(counts[0], total[0]),
					percent(counts[1], total[1])
				);
			}
			out.push('\n');
		}
		out
	}

	fn counters(&self) -> [&[u64; 2]; 16] {
		[
			&self.reads,
			&self.primary,
			&self.secondary,
			&self.supplementary,
			&self.duplicates,
			&self.primary_duplicates,
			&self.mapped,
			&self.primary_mapped,
			&self.paired,
			&self.read1,
			&self.read2,
			&self.properly_paired,
			&self.both_mapped,
			&self.singletons,
			&self.different_contig,
			&self.different_contig_confident,
		]
	}

	fn counters_mut(&mut self) -> [&mut [u64; 2]; 16] {
		[
			&mut self.reads,
			&mut self.primary,
			&mut self.secondary,
			&mut self.supplementary,
			&mut self.duplicates,
			&mut self.primary_duplicates,
			&mut self.mapped,
			&mut self.primary_mapped,
			&mut self.paired,
			&mut self.read1,
			&mut self.read2,
			&mut self.properly_paired,
			&mut self.both_mapped,
			&mut self.singletons,
			&mut self.different_contig,
			&mut self.different_contig_confident,
		]
	}
}

fn bump(counter: &mut [u64; 2], class: usize) {
	if let Some(slot) = counter.get_mut(class) {
		*slot += 1;
	}
}

#[allow(clippy::cast_precision_loss)]
fn percent(count: u64, total: u64) -> String {
	if total == 0 {
		return "N/A".to_owned();
	}
	let ratio = count as f32 / total as f32;
	format!("{:.2}%", f64::from(ratio) * 100.0)
}

pub fn flagstat(reader: &mut BamReader, workers: &Workers) -> Result<Flagstat, Error> {
	let mut total = Flagstat::default();
	let mut batch = RecordBatch::new();
	while reader.read_batch(&mut batch)? {
		total.merge(&count_batch(&batch, workers)?);
	}
	Ok(total)
}

fn count_batch(batch: &RecordBatch, workers: &Workers) -> Result<Flagstat, RecordError> {
	workers.install(|| {
		(0..batch.len())
			.into_par_iter()
			.try_fold(Flagstat::default, |mut stats, index| {
				stats.add(batch.at(index)?.record);
				Ok(stats)
			})
			.try_reduce(Flagstat::default, |mut left, right| {
				left.merge(&right);
				Ok(left)
			})
	})
}
