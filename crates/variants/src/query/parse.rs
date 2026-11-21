use super::{Field, Segment};
use crate::Error;
use crate::header::Header;

pub(super) fn parse(format: &str, header: &Header) -> Result<Vec<Segment>, Error> {
	let bytes = format.as_bytes();
	let mut at = 0;
	let segments = segments(bytes, &mut at, header, false)?;
	if at != bytes.len() {
		return Err(error(at, "unbalanced bracket"));
	}
	Ok(segments)
}

fn error(position: usize, reason: &'static str) -> Error {
	Error::FormatString { position, reason }
}

fn segments(
	bytes: &[u8],
	at: &mut usize,
	header: &Header,
	per_sample: bool,
) -> Result<Vec<Segment>, Error> {
	let mut out = Vec::new();
	let mut literal = Vec::new();
	while let Some(&byte) = bytes.get(*at) {
		match byte {
			b']' if per_sample => break,
			b']' => return Err(error(*at, "unbalanced bracket")),
			b'[' if per_sample => return Err(error(*at, "nested brackets")),
			b'[' => {
				flush(&mut literal, &mut out);
				*at += 1;
				let inner = segments(bytes, at, header, true)?;
				if bytes.get(*at) != Some(&b']') {
					return Err(error(*at, "unbalanced bracket"));
				}
				*at += 1;
				out.push(Segment::PerSample(inner));
			}
			b'%' => {
				flush(&mut literal, &mut out);
				*at += 1;
				out.push(Segment::Field(field(bytes, at, header, per_sample)?));
			}
			b'\\' => {
				*at += 1;
				let escaped = bytes
					.get(*at)
					.ok_or_else(|| error(*at, "dangling escape"))?;
				literal.push(match escaped {
					b'n' => b'\n',
					b't' => b'\t',
					b'r' => b'\r',
					other => *other,
				});
				*at += 1;
			}
			other => {
				literal.push(other);
				*at += 1;
			}
		}
	}
	flush(&mut literal, &mut out);
	Ok(out)
}

fn flush(literal: &mut Vec<u8>, out: &mut Vec<Segment>) {
	if !literal.is_empty() {
		out.push(Segment::Literal(std::mem::take(literal)));
	}
}

fn field(bytes: &[u8], at: &mut usize, header: &Header, per_sample: bool) -> Result<Field, Error> {
	let start = *at;
	let length = bytes
		.get(start..)
		.unwrap_or_default()
		.iter()
		.take_while(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'/' | b'.'))
		.count();
	let name = std::str::from_utf8(bytes.get(start..start + length).unwrap_or_default())
		.map_err(|_| error(start, "tag is not utf-8"))?;
	*at = start + length;
	let subscript = subscript(bytes, at)?;
	let known = match name {
		"CHROM" => Some(Field::Chrom),
		"POS" => Some(Field::Pos),
		"POS0" => Some(Field::Pos0),
		"END" => Some(Field::End),
		"END0" => Some(Field::End0),
		"ID" => Some(Field::Id),
		"REF" => Some(Field::Ref),
		"ALT" => Some(Field::Alt(subscript)),
		"FIRST_ALT" => Some(Field::FirstAlt),
		"QUAL" => Some(Field::Qual),
		"FILTER" => Some(Field::Filter),
		"INFO" => Some(Field::InfoAll),
		"FORMAT" => Some(Field::FormatKeys),
		"GT" if per_sample => Some(Field::Genotype),
		"TGT" if per_sample => Some(Field::TargetGenotype),
		"SAMPLE" if per_sample => Some(Field::Sample),
		"TYPE" => Some(Field::Type),
		"LINE" => Some(Field::Line),
		_ => None,
	};
	if let Some(field) = known {
		return Ok(field);
	}
	let (scope, tag) = match name.split_once('/') {
		Some(("INFO", tag)) => (Some(false), tag),
		Some(("FORMAT" | "FMT", tag)) => (Some(true), tag),
		_ => (None, name),
	};
	let format = scope.unwrap_or(per_sample);
	let key = header
		.key(tag)
		.ok_or_else(|| Error::UndefinedTag(tag.to_owned()))?;
	if format {
		header
			.type_of_format(key)
			.ok_or_else(|| Error::UndefinedTag(tag.to_owned()))?;
		Ok(Field::Format { key, subscript })
	} else {
		header
			.type_of_info(key)
			.ok_or_else(|| Error::UndefinedTag(tag.to_owned()))?;
		Ok(Field::Info { key, subscript })
	}
}

fn subscript(bytes: &[u8], at: &mut usize) -> Result<Option<usize>, Error> {
	if bytes.get(*at) != Some(&b'{') {
		return Ok(None);
	}
	let start = *at + 1;
	let length = bytes
		.get(start..)
		.unwrap_or_default()
		.iter()
		.take_while(|byte| byte.is_ascii_digit())
		.count();
	let end = start + length;
	if bytes.get(end) != Some(&b'}') || length == 0 {
		return Err(error(*at, "invalid subscript"));
	}
	*at = end + 1;
	std::str::from_utf8(bytes.get(start..end).unwrap_or_default())
		.ok()
		.and_then(|digits| digits.parse().ok())
		.map(Some)
		.ok_or_else(|| error(start, "invalid subscript"))
}
