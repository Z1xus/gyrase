use std::collections::VecDeque;
use std::io::Write;

use gyrase_index::{Builder, Entry, Geometry, Index};
use gyrase_io::{BlockPosition, Chunk, Level, VirtualOffset, Workers};

use super::BamWriter;
use super::index::entry;
use crate::error::{BamError, Error};
use crate::header::Header;
use crate::record::RecordRef;

#[derive(Debug)]
struct Pending {
	start: BlockPosition,
	end: BlockPosition,
	entry: Entry,
}

#[derive(Debug)]
pub(crate) struct IndexedWriter<W: Write + Send + 'static> {
	writer: BamWriter<W>,
	builder: Builder,
	pending: VecDeque<Pending>,
	references: u32,
}

impl<W: Write + Send + 'static> IndexedWriter<W> {
	pub(crate) fn new(
		inner: W,
		header: &Header,
		geometry: Geometry,
		level: Level,
		workers: &Workers,
	) -> Result<Self, Error> {
		Ok(Self {
			writer: BamWriter::new(inner, header, level, workers)?,
			builder: Builder::new(geometry),
			pending: VecDeque::new(),
			references: u32::try_from(header.references().len()).unwrap_or(u32::MAX),
		})
	}

	pub(crate) fn write(&mut self, record: RecordRef<'_>) -> Result<(), Error> {
		let written = self.writer.write(record)?;
		let unknown = VirtualOffset::new(0, 0);
		self.pending.push_back(Pending {
			start: written.start,
			end: written.end,
			entry: entry(
				record,
				Chunk {
					begin: unknown,
					end: unknown,
				},
			),
		});
		while let Some(front) = self.pending.front() {
			let (Some(begin), Some(end)) = (
				self.writer.resolve(front.start),
				self.writer.resolve(front.end),
			) else {
				break;
			};
			let mut entry = front.entry;
			entry.offsets = Chunk { begin, end };
			self.builder.push(&entry)?;
			self.pending.pop_front();
		}
		Ok(())
	}

	pub(crate) fn finish(self) -> Result<(W, Index), Error> {
		let Self {
			writer,
			mut builder,
			pending,
			references,
		} = self;
		let last = writer.position();
		let finished = writer.finish()?;
		for mut item in pending {
			if item.end == last && last.within != 0 {
				item.end = BlockPosition {
					ordinal: last.ordinal + 1,
					within: 0,
				};
			}
			let (Some(begin), Some(end)) =
				(finished.resolve(item.start), finished.resolve(item.end))
			else {
				return Err(BamError::UnresolvedOffset.into());
			};
			item.entry.offsets = Chunk { begin, end };
			builder.push(&item.entry)?;
		}
		Ok((finished.inner, builder.finish(references)?))
	}
}
