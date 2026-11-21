use super::fields::{Field, Fields};
use crate::error::HeaderError;

#[derive(Debug)]
pub(crate) enum ParsedLine<'a> {
	Header(Fields),
	Reference(Fields),
	ReadGroup(Fields),
	Program(Fields),
	Comment(&'a str),
	Other(&'a str),
}

pub(crate) fn parse_line(text: &str, line: usize) -> Result<ParsedLine<'_>, HeaderError> {
	let malformed = || HeaderError::MalformedLine { line };
	let body = text.strip_prefix('@').ok_or_else(malformed)?;
	let (code, rest) = match body.split_once('\t') {
		Some((code, rest)) => (code, Some(rest)),
		None => (body, None),
	};
	if code.len() != 2 || !code.bytes().all(|byte| byte.is_ascii_alphabetic()) {
		return Err(malformed());
	}
	Ok(match code {
		"HD" => ParsedLine::Header(parse_fields(rest, line)?),
		"SQ" => ParsedLine::Reference(parse_fields(rest, line)?),
		"RG" => ParsedLine::ReadGroup(parse_fields(rest, line)?),
		"PG" => ParsedLine::Program(parse_fields(rest, line)?),
		"CO" => ParsedLine::Comment(rest.unwrap_or_default()),
		_ => ParsedLine::Other(text),
	})
}

fn parse_fields(rest: Option<&str>, line: usize) -> Result<Fields, HeaderError> {
	let mut fields = Fields::new();
	for token in rest.into_iter().flat_map(|rest| rest.split('\t')) {
		let bytes = token.as_bytes();
		let [first, second, b':', ..] = *bytes else {
			return Err(HeaderError::MalformedLine { line });
		};
		fields.push(Field {
			tag: [first, second],
			value: token.get(3..).unwrap_or_default().to_owned(),
		});
	}
	Ok(fields)
}

pub(crate) fn write_fields(out: &mut String, code: &str, fields: &Fields) {
	out.push('@');
	out.push_str(code);
	for field in fields.iter() {
		out.push('\t');
		out.extend(field.tag.iter().map(|&byte| char::from(byte)));
		out.push(':');
		out.push_str(&field.value);
	}
	out.push('\n');
}
