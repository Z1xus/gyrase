use crate::reference::fill_loff;
use crate::tabix::FIXED_LEN;
use crate::wire::{Cursor, MAX_REFERENCES};
use crate::{Bin, Chunk, Error, Format, Geometry, Index, RefIndex};
use gyrase_io::VirtualOffset;

const REFERENCE_MIN_BYTES: usize = 4;
const CHUNK_BYTES: usize = 16;
const LINEAR_ENTRY_BYTES: usize = 8;
const BAI_BIN_MIN_BYTES: usize = 8;
const CSI_BIN_MIN_BYTES: usize = 16;

pub(crate) fn decode(bytes: &[u8]) -> Result<Index, Error> {
	let mut cursor = Cursor::new(bytes);
	match &cursor.array::<4>()? {
		b"BAI\x01" => decode_bai(&mut cursor),
		b"TBI\x01" => decode_tbi(&mut cursor),
		b"CSI\x01" => decode_csi(&mut cursor),
		_ => Err(Error::UnknownFormat),
	}
}

fn decode_bai(cursor: &mut Cursor<'_>) -> Result<Index, Error> {
	let count = cursor.count(REFERENCE_MIN_BYTES)?;
	let references = read_references(cursor, Geometry::BAI, Format::Bai, count)?;
	Ok(Index::assemble(
		Geometry::BAI,
		references,
		cursor.optional_u64(),
		Vec::new(),
	))
}

fn decode_tbi(cursor: &mut Cursor<'_>) -> Result<Index, Error> {
	let raw_count = cursor.u32()?;
	let fixed = cursor.take(FIXED_LEN)?;
	let names_len = Cursor::new(fixed.get(FIXED_LEN - 4..).ok_or(Error::Truncated)?).u32()?;
	let names = cursor.take(usize::try_from(names_len).map_err(|_| Error::Truncated)?)?;
	let count = cursor.validate_count(raw_count, REFERENCE_MIN_BYTES)?;
	let references = read_references(cursor, Geometry::BAI, Format::Tbi, count)?;
	Ok(Index::assemble(
		Geometry::BAI,
		references,
		cursor.optional_u64(),
		[fixed, names].concat(),
	))
}

fn decode_csi(cursor: &mut Cursor<'_>) -> Result<Index, Error> {
	let geometry = Geometry::new(cursor.u32()?, cursor.u32()?)?;
	let aux_len = usize::try_from(cursor.u32()?).map_err(|_| Error::Truncated)?;
	let aux = cursor.take(aux_len)?.to_vec();
	let count = cursor.count(REFERENCE_MIN_BYTES)?;
	let references = read_references(cursor, geometry, Format::Csi, count)?;
	Ok(Index::assemble(
		geometry,
		references,
		cursor.optional_u64(),
		aux,
	))
}

fn read_references(
	cursor: &mut Cursor<'_>,
	geometry: Geometry,
	format: Format,
	count: usize,
) -> Result<Vec<RefIndex>, Error> {
	if count > MAX_REFERENCES {
		return Err(Error::InvalidCount {
			count: u64::try_from(count).unwrap_or(u64::MAX),
		});
	}
	let mut references = Vec::new();
	for _ in 0..count {
		references.push(read_reference(cursor, geometry, format)?);
	}
	Ok(references)
}

fn read_reference(
	cursor: &mut Cursor<'_>,
	geometry: Geometry,
	format: Format,
) -> Result<RefIndex, Error> {
	let bin_min_bytes = match format {
		Format::Csi => CSI_BIN_MIN_BYTES,
		Format::Bai | Format::Tbi => BAI_BIN_MIN_BYTES,
	};
	let bin_count = cursor.count(bin_min_bytes)?;
	let mut bins = Vec::new();
	for _ in 0..bin_count {
		bins.push(read_bin(cursor, format)?);
	}
	match format {
		Format::Csi => RefIndex::new(bins, Vec::new()),
		Format::Bai | Format::Tbi => {
			let linear = read_linear(cursor)?;
			fill_loff(geometry, &mut bins, &linear);
			RefIndex::new(bins, linear)
		}
	}
}

fn read_bin(cursor: &mut Cursor<'_>, format: Format) -> Result<Bin, Error> {
	let id = cursor.u32()?;
	let loff = match format {
		Format::Csi => VirtualOffset::from_raw(cursor.u64()?),
		Format::Bai | Format::Tbi => VirtualOffset::default(),
	};
	let chunk_count = cursor.count(CHUNK_BYTES)?;
	let mut chunks = Vec::new();
	for _ in 0..chunk_count {
		chunks.push(Chunk {
			begin: VirtualOffset::from_raw(cursor.u64()?),
			end: VirtualOffset::from_raw(cursor.u64()?),
		});
	}
	Ok(Bin { id, loff, chunks })
}

fn read_linear(cursor: &mut Cursor<'_>) -> Result<Vec<VirtualOffset>, Error> {
	let count = cursor.count(LINEAR_ENTRY_BYTES)?;
	let mut linear = Vec::new();
	for _ in 0..count {
		linear.push(VirtualOffset::from_raw(cursor.u64()?));
	}
	Ok(linear)
}
