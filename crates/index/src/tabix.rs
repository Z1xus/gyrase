use gyrase_genome::ContigId;

use crate::Error;
use crate::wire::Cursor;

pub(crate) const FIXED_LEN: usize = 28;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TabixHeader {
	pub preset: u32,
	pub sequence_column: u32,
	pub begin_column: u32,
	pub end_column: u32,
	pub comment: u32,
	pub skip_lines: u32,
	pub names: Vec<String>,
}

impl TabixHeader {
	pub fn parse(aux: &[u8]) -> Result<Self, Error> {
		let invalid = |_| Error::InvalidTabixHeader;
		let mut cursor = Cursor::new(aux);
		let preset = cursor.u32().map_err(invalid)?;
		let sequence_column = cursor.u32().map_err(invalid)?;
		let begin_column = cursor.u32().map_err(invalid)?;
		let end_column = cursor.u32().map_err(invalid)?;
		let comment = cursor.u32().map_err(invalid)?;
		let skip_lines = cursor.u32().map_err(invalid)?;
		let names_len = usize::try_from(cursor.u32().map_err(invalid)?)
			.map_err(|_| Error::InvalidTabixHeader)?;
		let names_block = cursor.take(names_len).map_err(invalid)?;
		Ok(Self {
			preset,
			sequence_column,
			begin_column,
			end_column,
			comment,
			skip_lines,
			names: parse_names(names_block)?,
		})
	}

	pub fn to_bytes(&self) -> Vec<u8> {
		let names_len: usize = self.names.iter().map(|name| name.len() + 1).sum();
		let mut bytes = Vec::with_capacity(FIXED_LEN + names_len);
		for value in [
			self.preset,
			self.sequence_column,
			self.begin_column,
			self.end_column,
			self.comment,
			self.skip_lines,
		] {
			bytes.extend_from_slice(&value.to_le_bytes());
		}
		bytes.extend_from_slice(&u32::try_from(names_len).unwrap_or(u32::MAX).to_le_bytes());
		for name in &self.names {
			bytes.extend_from_slice(name.as_bytes());
			bytes.push(0);
		}
		bytes
	}

	pub fn contig(&self, name: &str) -> Option<ContigId> {
		let position = self.names.iter().position(|candidate| candidate == name)?;
		u32::try_from(position).ok().map(ContigId::new)
	}
}

fn parse_names(block: &[u8]) -> Result<Vec<String>, Error> {
	let Some((&last, terminated)) = block.split_last() else {
		return Ok(Vec::new());
	};
	if last != 0 {
		return Err(Error::InvalidTabixHeader);
	}
	terminated
		.split(|&byte| byte == 0)
		.map(|name| String::from_utf8(name.to_vec()).map_err(|_| Error::InvalidTabixHeader))
		.collect()
}
