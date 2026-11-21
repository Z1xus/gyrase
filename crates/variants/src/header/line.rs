use crate::Error;
use crate::value::{Number, ValueType};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LineKind {
	FileFormat,
	Contig,
	Info,
	Format,
	Filter,
	Alt,
	Other,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldDef {
	pub number: Number,
	pub value_type: ValueType,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Line {
	pub kind: LineKind,
	pub text: String,
	pub id: Option<String>,
	pub dictionary_id: Option<u32>,
}

pub(crate) struct Parsed {
	pub line: Line,
	pub explicit_index: Option<u32>,
	pub definition: Option<FieldDef>,
	pub contig_length: Option<u64>,
}

pub(crate) fn parse_line(text: &str) -> Result<Parsed, Error> {
	let body = text
		.strip_prefix("##")
		.ok_or(Error::Header("meta line without ##"))?;
	let (key, value) = body
		.split_once('=')
		.ok_or(Error::Header("meta line without ="))?;
	let kind = match key {
		"fileformat" => LineKind::FileFormat,
		"contig" => LineKind::Contig,
		"INFO" => LineKind::Info,
		"FORMAT" => LineKind::Format,
		"FILTER" => LineKind::Filter,
		"ALT" => LineKind::Alt,
		_ => LineKind::Other,
	};
	let inner = match kind {
		LineKind::FileFormat | LineKind::Other => None,
		_ => value
			.strip_prefix('<')
			.and_then(|rest| rest.strip_suffix('>')),
	};
	let Some(inner) = inner else {
		return Ok(Parsed {
			line: Line {
				kind: if kind == LineKind::FileFormat {
					kind
				} else {
					LineKind::Other
				},
				text: text.to_owned(),
				id: None,
				dictionary_id: None,
			},
			explicit_index: None,
			definition: None,
			contig_length: None,
		});
	};
	let pairs = split_pairs(inner);
	let find = |wanted: &str| {
		pairs
			.iter()
			.find(|(name, _)| *name == wanted)
			.map(|(_, value)| *value)
	};
	let id = find("ID")
		.ok_or(Error::Header("structured line without ID"))?
		.to_owned();
	let explicit_index = find("IDX").and_then(|value| value.parse().ok());
	let contig_length = find("length").and_then(|value| value.parse().ok());
	let definition = match kind {
		LineKind::Info | LineKind::Format => Some(FieldDef {
			number: find("Number")
				.and_then(Number::parse)
				.ok_or(Error::Header("field without a valid Number"))?,
			value_type: find("Type")
				.and_then(ValueType::parse)
				.ok_or(Error::Header("field without a valid Type"))?,
		}),
		_ => None,
	};
	let rebuilt = if explicit_index.is_some() || pairs.iter().any(|(name, _)| *name == "IDX") {
		rebuild(key, &pairs)
	} else {
		text.to_owned()
	};
	Ok(Parsed {
		line: Line {
			kind,
			text: rebuilt,
			id: Some(id),
			dictionary_id: None,
		},
		explicit_index,
		definition,
		contig_length,
	})
}

fn rebuild(key: &str, pairs: &[(&str, &str)]) -> String {
	let body: Vec<String> = pairs
		.iter()
		.filter(|(name, _)| *name != "IDX")
		.map(|(name, value)| format!("{name}={value}"))
		.collect();
	format!("##{key}=<{}>", body.join(","))
}

fn split_pairs(inner: &str) -> Vec<(&str, &str)> {
	let bytes = inner.as_bytes();
	let mut pairs = Vec::new();
	let mut start = 0;
	let mut quoted = false;
	let mut escaped = false;
	for (index, byte) in bytes.iter().enumerate() {
		match byte {
			b'\\' if quoted => escaped = !escaped,
			b'"' if !escaped => quoted = !quoted,
			b',' if !quoted => {
				push_pair(&mut pairs, inner, start, index);
				start = index + 1;
			}
			_ => escaped = false,
		}
	}
	push_pair(&mut pairs, inner, start, bytes.len());
	pairs
}

fn push_pair<'a>(pairs: &mut Vec<(&'a str, &'a str)>, inner: &'a str, start: usize, end: usize) {
	if let Some((name, value)) = inner.get(start..end).and_then(|pair| pair.split_once('=')) {
		pairs.push((name, value));
	}
}
