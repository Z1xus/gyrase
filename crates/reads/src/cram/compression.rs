use crate::cram::Error;
use crate::cram::cursor::Cursor;
use crate::cram::encoding::{Encoding, Slots};

const MAX_TAG_LINES: usize = 1 << 16;
const MAX_MAP_ENTRIES: usize = 1 << 16;

fn malformed() -> Error {
	Error::Malformed("compression header")
}

#[derive(Debug)]
pub(crate) struct Substitutions([[u8; 4]; 5]);

impl Substitutions {
	const BASES: [u8; 5] = *b"ACGTN";

	fn standard() -> Self {
		Self([*b"CGTN", *b"AGTN", *b"ACTN", *b"ACGN", *b"ACGT"])
	}

	fn read(bytes: [u8; 5]) -> Self {
		let mut matrix = [[b'N'; 4]; 5];
		for ((row, byte), reference) in matrix.iter_mut().zip(bytes).zip(Self::BASES) {
			let alternatives = Self::BASES
				.iter()
				.copied()
				.filter(|base| *base != reference);
			for (rank, base) in alternatives.enumerate() {
				let shift = 6 - 2 * rank;
				if let Some(slot) = row.get_mut(usize::from(byte >> shift) & 3) {
					*slot = base;
				}
			}
		}
		Self(matrix)
	}

	pub(crate) fn base(&self, reference: u8, code: u8) -> u8 {
		let row = match reference.to_ascii_uppercase() {
			b'A' => 0,
			b'C' => 1,
			b'G' => 2,
			b'T' => 3,
			_ => 4,
		};
		self.0
			.get(row)
			.and_then(|bases| bases.get(usize::from(code)))
			.copied()
			.unwrap_or(b'N')
	}
}

#[derive(Debug)]
pub(crate) struct TagBinding {
	pub(crate) name: [u8; 2],
	pub(crate) kind: u8,
	pub(crate) encoding: usize,
}

#[derive(Debug, Default)]
pub(crate) struct Series {
	pub(crate) flags: Encoding,
	pub(crate) compression_flags: Encoding,
	pub(crate) reference_id: Encoding,
	pub(crate) read_length: Encoding,
	pub(crate) position: Encoding,
	pub(crate) read_group: Encoding,
	pub(crate) name: Encoding,
	pub(crate) mate_flags: Encoding,
	pub(crate) mate_reference: Encoding,
	pub(crate) mate_position: Encoding,
	pub(crate) template_length: Encoding,
	pub(crate) next_fragment: Encoding,
	pub(crate) tag_line: Encoding,
	pub(crate) feature_count: Encoding,
	pub(crate) feature_code: Encoding,
	pub(crate) feature_position: Encoding,
	pub(crate) deletion_length: Encoding,
	pub(crate) base: Encoding,
	pub(crate) substitution: Encoding,
	pub(crate) insertion: Encoding,
	pub(crate) soft_clip: Encoding,
	pub(crate) hard_clip: Encoding,
	pub(crate) padding: Encoding,
	pub(crate) reference_skip: Encoding,
	pub(crate) mapping_quality: Encoding,
	pub(crate) bases: Encoding,
	pub(crate) quality_stretch: Encoding,
	pub(crate) quality: Encoding,
}

impl Series {
	fn slot(&mut self, key: [u8; 2]) -> Option<&mut Encoding> {
		Some(match &key {
			b"BF" => &mut self.flags,
			b"CF" => &mut self.compression_flags,
			b"RI" => &mut self.reference_id,
			b"RL" => &mut self.read_length,
			b"AP" => &mut self.position,
			b"RG" => &mut self.read_group,
			b"RN" => &mut self.name,
			b"MF" => &mut self.mate_flags,
			b"NS" => &mut self.mate_reference,
			b"NP" => &mut self.mate_position,
			b"TS" => &mut self.template_length,
			b"NF" => &mut self.next_fragment,
			b"TL" => &mut self.tag_line,
			b"FN" => &mut self.feature_count,
			b"FC" => &mut self.feature_code,
			b"FP" => &mut self.feature_position,
			b"DL" => &mut self.deletion_length,
			b"BA" => &mut self.base,
			b"BS" => &mut self.substitution,
			b"IN" => &mut self.insertion,
			b"SC" => &mut self.soft_clip,
			b"HC" => &mut self.hard_clip,
			b"PD" => &mut self.padding,
			b"RS" => &mut self.reference_skip,
			b"MQ" => &mut self.mapping_quality,
			b"BB" => &mut self.bases,
			b"QQ" => &mut self.quality_stretch,
			b"QS" => &mut self.quality,
			_ => return None,
		})
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ReadNames {
	Stored,
	Generated,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Positions {
	Delta,
	Absolute,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Bases {
	FromReference,
	Explicit,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Qualities {
	Sequenced,
	Aligned,
}

#[derive(Debug)]
pub(crate) struct CompressionHeader {
	pub(crate) read_names: ReadNames,
	pub(crate) positions: Positions,
	pub(crate) bases: Bases,
	pub(crate) qualities: Qualities,
	pub(crate) substitutions: Substitutions,
	pub(crate) tag_lines: Vec<Vec<TagBinding>>,
	pub(crate) tag_encodings: Vec<Encoding>,
	pub(crate) series: Series,
	pub(crate) slots: Vec<i32>,
}

impl CompressionHeader {
	pub(crate) fn read(bytes: &[u8]) -> Result<Self, Error> {
		let mut cursor = Cursor::new(bytes);
		let mut header = Self {
			read_names: ReadNames::Generated,
			positions: Positions::Delta,
			bases: Bases::FromReference,
			qualities: Qualities::Sequenced,
			substitutions: Substitutions::standard(),
			tag_lines: Vec::new(),
			tag_encodings: Vec::new(),
			series: Series::default(),
			slots: Vec::new(),
		};
		let mut slots = Slots::default();
		let mut dictionary = Vec::new();
		header.read_preservation(&mut cursor, &mut dictionary)?;
		header.read_series(&mut cursor, &mut slots)?;
		let tag_ids = header.read_tags(&mut cursor, &mut slots)?;
		header.tag_lines = dictionary
			.iter()
			.map(|line| bind_line(line, &tag_ids))
			.collect::<Result<_, _>>()?;
		header.slots = slots.ids().to_vec();
		Ok(header)
	}

	fn read_preservation(
		&mut self,
		cursor: &mut Cursor<'_>,
		dictionary: &mut Vec<Vec<[u8; 3]>>,
	) -> Result<(), Error> {
		let (mut map, count) = map_prologue(cursor)?;
		for _ in 0..count {
			match &map.array::<2>()? {
				b"RN" => {
					self.read_names = if map.u8()? == 0 {
						ReadNames::Generated
					} else {
						ReadNames::Stored
					};
				}
				b"AP" => {
					self.positions = if map.u8()? == 0 {
						Positions::Absolute
					} else {
						Positions::Delta
					};
				}
				b"RR" => {
					self.bases = if map.u8()? == 0 {
						Bases::Explicit
					} else {
						Bases::FromReference
					};
				}
				b"QO" => {
					self.qualities = if map.u8()? == 0 {
						Qualities::Aligned
					} else {
						Qualities::Sequenced
					};
				}
				b"SM" => self.substitutions = Substitutions::read(map.array()?),
				b"TD" => *dictionary = read_dictionary(&mut map)?,
				_ => {
					map.u8()?;
				}
			}
		}
		finish_map(&map)
	}

	fn read_series(&mut self, cursor: &mut Cursor<'_>, slots: &mut Slots) -> Result<(), Error> {
		let (mut map, count) = map_prologue(cursor)?;
		for _ in 0..count {
			let key = map.array::<2>()?;
			let encoding = Encoding::read(&mut map, slots)?;
			if let Some(target) = self.series.slot(key) {
				*target = encoding;
			}
		}
		finish_map(&map)
	}

	fn read_tags(
		&mut self,
		cursor: &mut Cursor<'_>,
		slots: &mut Slots,
	) -> Result<Vec<[u8; 3]>, Error> {
		let (mut map, count) = map_prologue(cursor)?;
		let mut ids = Vec::new();
		for _ in 0..count {
			let id = u32::try_from(map.itf8()?).map_err(|_| malformed())?;
			let [_, first, second, kind] = id.to_be_bytes();
			ids.push([first, second, kind]);
			self.tag_encodings.push(Encoding::read(&mut map, slots)?);
		}
		finish_map(&map)?;
		Ok(ids)
	}
}

fn map_prologue<'a>(cursor: &mut Cursor<'a>) -> Result<(Cursor<'a>, usize), Error> {
	let size = usize::try_from(cursor.itf8()?).map_err(|_| malformed())?;
	let mut map = Cursor::new(cursor.take(size)?);
	let count = usize::try_from(map.itf8()?)
		.ok()
		.filter(|count| *count <= MAX_MAP_ENTRIES)
		.ok_or_else(malformed)?;
	Ok((map, count))
}

fn finish_map(map: &Cursor<'_>) -> Result<(), Error> {
	if map.is_empty() {
		Ok(())
	} else {
		Err(malformed())
	}
}

fn read_dictionary(cursor: &mut Cursor<'_>) -> Result<Vec<Vec<[u8; 3]>>, Error> {
	let size = usize::try_from(cursor.itf8()?).map_err(|_| malformed())?;
	let bytes = cursor.take(size)?;
	let mut lines = Vec::new();
	for line in bytes.split(|byte| *byte == 0).take(MAX_TAG_LINES) {
		if line.len() % 3 != 0 {
			return Err(malformed());
		}
		lines.push(line.as_chunks::<3>().0.to_vec());
	}
	if bytes.last() == Some(&0) {
		lines.pop();
	}
	Ok(lines)
}

fn bind_line(line: &[[u8; 3]], ids: &[[u8; 3]]) -> Result<Vec<TagBinding>, Error> {
	line.iter()
		.map(|[first, second, kind]| {
			let encoding = ids
				.iter()
				.position(|id| id == &[*first, *second, *kind])
				.ok_or_else(malformed)?;
			Ok(TagBinding {
				name: [*first, *second],
				kind: *kind,
				encoding,
			})
		})
		.collect()
}
