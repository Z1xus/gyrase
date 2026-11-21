use std::collections::HashMap;
use std::mem;

use gyrase_genome::{ContigId, Interval};
use gyrase_io::VirtualOffset;

use crate::finalize::{Draft, finalize};
use crate::wire::MAX_REFERENCES;
use crate::{Chunk, Error, Geometry, Index, RefIndex};

const UNSET: u64 = u64::MAX;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Placement {
	Unplaced,
	Placed(Interval),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Entry {
	pub placement: Placement,
	pub offsets: Chunk,
	pub mapped: bool,
}

#[derive(Debug)]
struct Run {
	contig: ContigId,
	bins: Vec<(u32, Vec<Chunk>)>,
	puts: Vec<u32>,
	positions: HashMap<u32, usize>,
	linear: Vec<u64>,
	open_bin: u32,
	open_begin: VirtualOffset,
	region_begin: VirtualOffset,
	last_start: u64,
	last_end: VirtualOffset,
	mapped: u64,
	unmapped: u64,
}

#[derive(Debug)]
enum State {
	Empty,
	Placed(Run),
	Unplaced,
}

#[derive(Debug)]
pub struct Builder {
	geometry: Geometry,
	completed: Vec<RefIndex>,
	state: State,
	unplaced: u64,
}

impl Run {
	fn new(geometry: Geometry, contig: ContigId, start: u64, end: u64, entry: &Entry) -> Self {
		let mut run = Self {
			contig,
			bins: Vec::new(),
			puts: Vec::new(),
			positions: HashMap::new(),
			linear: Vec::new(),
			open_bin: geometry.bin(start, end),
			open_begin: entry.offsets.begin,
			region_begin: entry.offsets.begin,
			last_start: start,
			last_end: entry.offsets.end,
			mapped: 0,
			unmapped: 0,
		};
		run.cover(geometry, start, end, entry);
		run
	}

	fn add(
		&mut self,
		geometry: Geometry,
		start: u64,
		end: u64,
		entry: &Entry,
	) -> Result<(), Error> {
		if start < self.last_start {
			return Err(Error::Unsorted {
				contig: self.contig.get(),
				previous: self.last_start,
				start,
			});
		}
		let bin = geometry.bin(start, end);
		if bin != self.open_bin {
			let closed = Chunk {
				begin: self.open_begin,
				end: entry.offsets.begin,
			};
			self.add_chunk(self.open_bin, closed);
			self.open_bin = bin;
			self.open_begin = entry.offsets.begin;
		}
		self.cover(geometry, start, end, entry);
		Ok(())
	}

	fn cover(&mut self, geometry: Geometry, start: u64, end: u64, entry: &Entry) {
		let shift = geometry.min_shift();
		let first_window = start >> shift;
		let last_window = end.saturating_sub(1) >> shift;
		let windows =
			usize::try_from(last_window).map_or(usize::MAX, |last| last.saturating_add(1));
		if self.linear.len() < windows {
			self.linear.resize(windows, UNSET);
		}
		for window in first_window..=last_window {
			if let Some(slot) = usize::try_from(window)
				.ok()
				.and_then(|window| self.linear.get_mut(window))
				&& *slot == UNSET
			{
				*slot = entry.offsets.begin.raw();
			}
		}
		if entry.mapped {
			self.mapped += 1;
		} else {
			self.unmapped += 1;
		}
		self.last_start = start;
		self.last_end = entry.offsets.end;
	}

	fn add_chunk(&mut self, bin: u32, chunk: Chunk) {
		self.puts.push(bin);
		let position = *self.positions.entry(bin).or_insert(self.bins.len());
		if let Some((_, chunks)) = self.bins.get_mut(position) {
			chunks.push(chunk);
		} else {
			self.bins.push((bin, vec![chunk]));
		}
	}

	fn close(mut self, geometry: Geometry, closing: VirtualOffset) -> Result<RefIndex, Error> {
		self.add_chunk(
			self.open_bin,
			Chunk {
				begin: self.open_begin,
				end: closing,
			},
		);
		let metadata_bin = geometry.metadata_bin();
		self.add_chunk(
			metadata_bin,
			Chunk {
				begin: self.region_begin,
				end: closing,
			},
		);
		self.add_chunk(
			metadata_bin,
			Chunk {
				begin: VirtualOffset::from_raw(self.mapped),
				end: VirtualOffset::from_raw(self.unmapped),
			},
		);
		finalize(
			geometry,
			Draft {
				bins: self.bins,
				puts: self.puts,
				linear: self.linear,
			},
		)
	}
}

impl Builder {
	pub fn new(geometry: Geometry) -> Self {
		Self {
			geometry,
			completed: Vec::new(),
			state: State::Empty,
			unplaced: 0,
		}
	}

	pub fn push(&mut self, entry: &Entry) -> Result<(), Error> {
		match entry.placement {
			Placement::Unplaced => self.push_unplaced(entry),
			Placement::Placed(interval) => self.push_placed(interval, entry),
		}
	}

	pub fn finish(mut self, references: u32) -> Result<Index, Error> {
		if let State::Placed(run) = mem::replace(&mut self.state, State::Empty) {
			let closing = run.last_end;
			self.complete(run, closing)?;
		}
		let references = usize::try_from(references)
			.ok()
			.filter(|&count| count <= MAX_REFERENCES)
			.ok_or(Error::ReferenceOutOfRange {
				contig: references,
				limit: MAX_REFERENCES,
			})?;
		if self.completed.len() > references {
			return Err(Error::ReferenceOutOfRange {
				contig: u32::try_from(self.completed.len()).unwrap_or(u32::MAX),
				limit: references,
			});
		}
		self.completed.resize_with(references, RefIndex::default);
		Ok(Index::assemble(
			self.geometry,
			self.completed,
			self.unplaced,
			Vec::new(),
		))
	}

	fn push_placed(&mut self, interval: Interval, entry: &Entry) -> Result<(), Error> {
		let end = self.validate(interval)?;
		match &mut self.state {
			State::Unplaced => return Err(Error::PlacedAfterUnplaced),
			State::Placed(run) if run.contig == interval.contig => {
				return run.add(self.geometry, interval.start, end, entry);
			}
			State::Placed(run) if interval.contig < run.contig => {
				return Err(Error::ContigOrder {
					contig: interval.contig.get(),
					previous: run.contig.get(),
				});
			}
			State::Placed(_) | State::Empty => {}
		}
		if let State::Placed(previous) = mem::replace(&mut self.state, State::Empty) {
			self.complete(previous, entry.offsets.begin)?;
		}
		self.state = State::Placed(Run::new(
			self.geometry,
			interval.contig,
			interval.start,
			end,
			entry,
		));
		Ok(())
	}

	fn push_unplaced(&mut self, entry: &Entry) -> Result<(), Error> {
		if let State::Placed(previous) = mem::replace(&mut self.state, State::Unplaced) {
			self.complete(previous, entry.offsets.begin)?;
		}
		self.unplaced += 1;
		Ok(())
	}

	fn validate(&self, interval: Interval) -> Result<u64, Error> {
		if interval.end < interval.start {
			return Err(Error::InvalidInterval {
				start: interval.start,
				end: interval.end,
			});
		}
		let limit = self.geometry.max_position();
		if interval.end > limit {
			return Err(Error::PositionOutOfRange {
				position: interval.end,
				limit,
			});
		}
		if interval.contig.index() >= MAX_REFERENCES {
			return Err(Error::ReferenceOutOfRange {
				contig: interval.contig.get(),
				limit: MAX_REFERENCES,
			});
		}
		Ok(interval.end.max(1))
	}

	fn complete(&mut self, run: Run, closing: VirtualOffset) -> Result<(), Error> {
		let position = run.contig.index();
		let reference = run.close(self.geometry, closing)?;
		self.completed.resize_with(position, RefIndex::default);
		self.completed.push(reference);
		Ok(())
	}
}
