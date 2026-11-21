use std::borrow::Cow;

use crate::cram::Error;
use crate::cram::codec::Method;
use crate::cram::compression::CompressionHeader;
use crate::cram::cursor::Cursor;
use crate::cram::encoding::{BitReader, Data, External};
use crate::cram::structure::{Block, BlockKind, ContainerHeader};

const MAX_SLICE_BLOCKS: usize = 10_000;

fn malformed() -> Error {
	Error::Malformed("slice header")
}

#[derive(Debug)]
pub(crate) struct SliceHeader {
	pub(crate) reference: i32,
	pub(crate) start: i32,
	pub(crate) span: i32,
	pub(crate) records: u32,
	pub(crate) record_counter: i64,
	pub(crate) blocks: usize,
	pub(crate) embedded: Option<i32>,
	pub(crate) checksum: [u8; 16],
}

impl SliceHeader {
	pub(crate) fn read(bytes: &[u8]) -> Result<Self, Error> {
		let mut cursor = Cursor::new(bytes);
		let reference = cursor.itf8()?;
		let start = cursor.itf8()?;
		let span = cursor.itf8()?;
		if start < 0 || span < 0 {
			return Err(malformed());
		}
		let records = u32::try_from(cursor.itf8()?).map_err(|_| malformed())?;
		let record_counter = cursor.ltf8()?;
		let blocks = usize::try_from(cursor.itf8()?)
			.ok()
			.filter(|count| (1..MAX_SLICE_BLOCKS).contains(count))
			.ok_or_else(malformed)?;
		let content_ids = usize::try_from(cursor.itf8()?)
			.ok()
			.filter(|count| (1..MAX_SLICE_BLOCKS).contains(count))
			.ok_or_else(malformed)?;
		for _ in 0..content_ids {
			cursor.itf8()?;
		}
		let embedded = Some(cursor.itf8()?).filter(|id| *id >= 0);
		let checksum = cursor.array()?;
		Ok(Self {
			reference,
			start,
			span,
			records,
			record_counter,
			blocks,
			embedded,
			checksum,
		})
	}
}

#[derive(Debug)]
pub(crate) struct SliceBlocks<'a> {
	pub(crate) header: SliceHeader,
	pub(crate) blocks: Vec<Block<'a>>,
}

#[derive(Debug)]
pub(crate) struct Parsed<'a> {
	pub(crate) compression: CompressionHeader,
	pub(crate) slices: Vec<SliceBlocks<'a>>,
}

pub(crate) fn parse<'a>(header: &ContainerHeader, body: &'a [u8]) -> Result<Parsed<'a>, Error> {
	let mut cursor = Cursor::new(body);
	let first = Block::read(&mut cursor)?;
	if first.kind != BlockKind::CompressionHeader {
		return Err(Error::Malformed("container"));
	}
	let compression = CompressionHeader::read(&first.decode()?)?;
	let mut slices = Vec::with_capacity(header.landmarks.len());
	for _ in 0..header.landmarks.len() {
		let head = Block::read(&mut cursor)?;
		if head.kind != BlockKind::SliceHeader {
			return Err(Error::Malformed("container"));
		}
		let header = SliceHeader::read(&head.decode()?)?;
		let mut blocks = Vec::with_capacity(header.blocks);
		for _ in 0..header.blocks {
			blocks.push(Block::read(&mut cursor)?);
		}
		slices.push(SliceBlocks { header, blocks });
	}
	Ok(Parsed {
		compression,
		slices,
	})
}

#[derive(Debug)]
pub(crate) struct Payload<'a> {
	core: Cow<'a, [u8]>,
	externals: Vec<(Cow<'a, [u8]>, Option<u8>)>,
	pub(crate) embedded: Option<Cow<'a, [u8]>>,
}

impl<'a> SliceBlocks<'a> {
	pub(crate) fn payload(&self, slots: &[i32]) -> Result<Payload<'a>, Error> {
		let core = match self
			.blocks
			.iter()
			.find(|block| block.kind == BlockKind::CoreData)
		{
			Some(block) => block.decode()?,
			None => Cow::Borrowed(&[][..]),
		};
		let find = |id: i32| {
			self.blocks
				.iter()
				.find(|block| block.kind == BlockKind::ExternalData && block.content_id == id)
		};
		let mut externals = Vec::with_capacity(slots.len());
		for id in slots {
			externals.push(match find(*id) {
				Some(block) => {
					let stop = (block.method == Method::Tok3).then_some(0);
					(block.decode()?, stop)
				}
				None => (Cow::Borrowed(&[][..]), None),
			});
		}
		let embedded = match self.header.embedded {
			Some(id) => Some(find(id).ok_or_else(malformed)?.decode()?),
			None => None,
		};
		Ok(Payload {
			core,
			externals,
			embedded,
		})
	}
}

impl Payload<'_> {
	pub(crate) fn data(&self) -> Data<'_> {
		Data::new(
			BitReader::new(&self.core),
			self.externals
				.iter()
				.map(|(bytes, stop)| External::new(bytes, *stop))
				.collect(),
		)
	}
}
