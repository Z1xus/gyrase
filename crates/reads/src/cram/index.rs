use std::io::{Read, Seek, SeekFrom};

use gyrase_genome::ContigId;
use gyrase_index::{Crai, CraiEntry, CraiReference};

use crate::cram::Error;
use crate::cram::cursor::Cursor;
use crate::cram::decode::{Environment, decode_extents};
use crate::cram::reader::parse_header;
use crate::cram::slice::{self, SliceHeader};
use crate::cram::structure::{Block, BlockKind, ContainerHeader, Version};

const SLICE_HEAD_PROBE: usize = 512;
const MULTI_REFERENCE: i32 = -2;

struct Place {
	container: u64,
	landmark: u32,
	size: usize,
}

pub fn index_cram<E: From<Error>>(
	mut source: impl Read + Seek,
	mut progress: impl FnMut(u64) -> Result<(), E>,
) -> Result<Crai, E> {
	Version::read(&mut source)?;
	let (_, file_header) = ContainerHeader::read(&mut source)?.ok_or(Error::Truncated)?;
	let environment = Environment::unbacked(&parse_header(&file_header)?);
	let mut entries = Vec::new();
	let mut last = (-9, -9);
	loop {
		let container_offset = source.stream_position().map_err(Error::from)?;
		let container = ContainerHeader::read_head(&mut source)?.ok_or(Error::Truncated)?;
		if container.is_eof() {
			break;
		}
		if container.reference == last.0 && container.start < last.1 {
			return Err(Error::NotSorted.into());
		}
		last = (container.reference, container.start);
		let body = container_offset + container.size as u64;
		for (index, landmark) in container.landmarks.iter().enumerate() {
			let next = container
				.landmarks
				.get(index + 1)
				.map_or(container.length, |next| *next as usize);
			let place = Place {
				container: container_offset,
				landmark: *landmark,
				size: next
					.checked_sub(*landmark as usize)
					.ok_or(Error::Malformed("container"))?,
			};
			let header = read_slice_header(&mut source, body + u64::from(*landmark), place.size)?;
			if header.reference == MULTI_REFERENCE {
				let extents = multi_extents(&mut source, &environment, &container, body, index)?;
				split(&extents, &place, &mut entries)?;
			} else {
				entries.push(CraiEntry {
					reference: reference(header.reference)?,
					start: unsigned(i64::from(header.start))?,
					span: unsigned(i64::from(header.span))?,
					container_offset: place.container,
					slice_offset: u64::from(place.landmark),
					slice_size: place.size as u64,
				});
			}
		}
		let end = body + container.length as u64;
		source.seek(SeekFrom::Start(end)).map_err(Error::from)?;
		progress(end)?;
	}
	Ok(Crai::new(entries))
}

fn read_slice_header(
	source: &mut (impl Read + Seek),
	at: u64,
	size: usize,
) -> Result<SliceHeader, Error> {
	source.seek(SeekFrom::Start(at))?;
	let mut buffer = Vec::new();
	let mut wanted = SLICE_HEAD_PROBE.min(size);
	loop {
		let before = buffer.len();
		buffer.resize(wanted, 0);
		source.read_exact(buffer.get_mut(before..).ok_or(Error::Truncated)?)?;
		match Block::read(&mut Cursor::new(&buffer)) {
			Ok(block) if block.kind == BlockKind::SliceHeader => {
				return SliceHeader::read(&block.decode()?);
			}
			Ok(_) => return Err(Error::Malformed("container")),
			Err(Error::Truncated) if wanted < size => wanted = (wanted * 2).min(size),
			Err(error) => return Err(error),
		}
	}
}

fn multi_extents(
	source: &mut (impl Read + Seek),
	environment: &Environment,
	container: &ContainerHeader,
	body: u64,
	slice_index: usize,
) -> Result<Vec<(i32, i32, i32)>, Error> {
	source.seek(SeekFrom::Start(body))?;
	let mut bytes = Vec::new();
	source
		.take(container.length as u64)
		.read_to_end(&mut bytes)?;
	if bytes.len() != container.length {
		return Err(Error::Truncated);
	}
	let parsed = slice::parse(container, &bytes)?;
	let blocks = parsed
		.slices
		.get(slice_index)
		.ok_or(Error::Malformed("container"))?;
	decode_extents(environment, &parsed.compression, blocks)
}

fn split(
	extents: &[(i32, i32, i32)],
	place: &Place,
	entries: &mut Vec<CraiEntry>,
) -> Result<(), Error> {
	let mut current: Option<(i32, i32, i32)> = None;
	let mut last = (-9, -9);
	for &(reference, position, end) in extents {
		if reference == last.0 && position < last.1 {
			return Err(Error::NotSorted);
		}
		last = (reference, position);
		match &mut current {
			Some((current_reference, _, current_end)) if *current_reference == reference => {
				*current_end = (*current_end).max(end);
			}
			_ => {
				if let Some(group) = current.take() {
					entries.push(group_entry(group, place)?);
				}
				current = Some((reference, position, end));
			}
		}
	}
	if let Some(group) = current {
		entries.push(group_entry(group, place)?);
	}
	Ok(())
}

fn group_entry((id, start, end): (i32, i32, i32), place: &Place) -> Result<CraiEntry, Error> {
	Ok(CraiEntry {
		reference: reference(id)?,
		start: unsigned(i64::from(start))?,
		span: unsigned(i64::from(end) - i64::from(start) + 1)?,
		container_offset: place.container,
		slice_offset: u64::from(place.landmark),
		slice_size: place.size as u64,
	})
}

fn reference(id: i32) -> Result<CraiReference, Error> {
	match id {
		-1 => Ok(CraiReference::Unmapped),
		_ => u32::try_from(id)
			.map(|id| CraiReference::Contig(ContigId::new(id)))
			.map_err(|_| Error::Malformed("slice header")),
	}
}

fn unsigned(value: i64) -> Result<u64, Error> {
	u64::try_from(value).map_err(|_| Error::Malformed("slice header"))
}
