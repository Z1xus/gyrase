use gyrase_genome::ContigId;
use memchr::memchr;

use crate::error::{Error, SamError, SamField};
use crate::header::Header;
use crate::record::{CigarKind, CigarOp, Element, Flags, QualityInput, Record, Value};

const MISSING: &[u8] = b"*";
const SAME_AS_REFERENCE: &[u8] = b"=";
const MANDATORY_FIELDS: usize = 11;
const MAX_OP_LEN: u64 = (1 << 28) - 1;
const MAX_POSITION: u64 = 1 << 31;

#[derive(Debug, Default)]
pub(super) struct Parser {
	ops: Vec<CigarOp>,
}

struct Line<'a> {
	number: u64,
	rest: &'a [u8],
}

impl<'a> Line<'a> {
	fn next_field(&mut self) -> Option<&'a [u8]> {
		if self.rest.is_empty() {
			return None;
		}
		let end = memchr(b'\t', self.rest).unwrap_or(self.rest.len());
		let (field, rest) = self.rest.split_at(end);
		self.rest = rest.get(1..).unwrap_or_default();
		Some(field)
	}

	fn mandatory(&mut self) -> Result<[&'a [u8]; MANDATORY_FIELDS], SamError> {
		let mut fields = [&[][..]; MANDATORY_FIELDS];
		for slot in &mut fields {
			*slot = self
				.next_field()
				.ok_or(SamError::TooFewFields { line: self.number })?;
		}
		Ok(fields)
	}

	fn invalid(&self, field: SamField) -> SamError {
		SamError::Invalid {
			line: self.number,
			field,
		}
	}
}

impl Parser {
	pub(super) fn parse(
		&mut self,
		text: &[u8],
		number: u64,
		header: &Header,
		record: &mut Record,
	) -> Result<(), Error> {
		let mut line = Line { number, rest: text };
		let [
			name,
			flags,
			reference,
			position,
			mapq,
			cigar,
			mate,
			mate_position,
			template,
			bases,
			qualities,
		] = line.mandatory()?;
		record.clear();
		if name.first() == Some(&b'@') {
			return Err(SamError::LateHeader { line: number }.into());
		}
		record
			.set_name(name)
			.map_err(|_| line.invalid(SamField::Name))?;
		let mut flags = parse_flags(&line, flags)?;
		let mut contig = parse_reference(&line, header, reference, SamField::Reference)?;
		let position = parse_position(&line, position, SamField::Position)?;
		if position.is_none() {
			contig = None;
		}
		if contig.is_none() {
			flags |= Flags::UNMAPPED;
		}
		let mapq = parse_unsigned(mapq)
			.and_then(|value| u8::try_from(value).ok())
			.ok_or_else(|| line.invalid(SamField::MapQ))?;
		self.parse_cigar(&line, cigar)?;
		if self.ops.is_empty() {
			flags |= Flags::UNMAPPED;
		}
		let mate_contig = if mate == SAME_AS_REFERENCE {
			contig
		} else {
			parse_reference(&line, header, mate, SamField::MateReference)?
		};
		let mate_position = parse_position(&line, mate_position, SamField::MatePosition)?;
		let mate_contig = mate_contig.filter(|_| mate_position.is_some());
		let template = parse_signed(template)
			.and_then(|value| i32::try_from(value).ok())
			.ok_or_else(|| line.invalid(SamField::TemplateLength))?;
		record.set_mapq(mapq);
		record.set_contig(contig)?;
		record.set_position(position)?;
		record.set_flags(flags);
		record.set_mate_contig(mate_contig)?;
		record.set_mate_position(mate_position)?;
		record.set_template_length(template);
		record.set_cigar(&self.ops)?;
		parse_read(&line, bases, qualities, record)?;
		check_cigar_length(&line, record)?;
		parse_tags(&mut line, record)
	}

	fn parse_cigar(&mut self, line: &Line<'_>, text: &[u8]) -> Result<(), SamError> {
		self.ops.clear();
		if text == MISSING {
			return Ok(());
		}
		let invalid = || line.invalid(SamField::Cigar);
		let mut rest = text;
		while !rest.is_empty() {
			let digits = rest.iter().take_while(|byte| byte.is_ascii_digit()).count();
			let (number, tail) = rest.split_at(digits);
			let (&symbol, tail) = tail.split_first().ok_or_else(invalid)?;
			let len = parse_unsigned(number)
				.filter(|&len| len <= MAX_OP_LEN)
				.and_then(|len| u32::try_from(len).ok())
				.ok_or_else(invalid)?;
			let kind = CigarKind::from_symbol(symbol).ok_or_else(invalid)?;
			self.ops.push(CigarOp::new(len, kind).ok_or_else(invalid)?);
			rest = tail;
		}
		Ok(())
	}
}

fn parse_unsigned(text: &[u8]) -> Option<u64> {
	if text.is_empty() {
		return None;
	}
	text.iter().try_fold(0u64, |value, &byte| {
		let digit = char::from(byte).to_digit(10)?;
		value.checked_mul(10)?.checked_add(u64::from(digit))
	})
}

fn parse_signed(text: &[u8]) -> Option<i64> {
	match text.split_first() {
		Some((b'-', digits)) => i64::try_from(parse_unsigned(digits)?)
			.ok()
			.map(|value| -value),
		Some((b'+', digits)) => i64::try_from(parse_unsigned(digits)?).ok(),
		_ => i64::try_from(parse_unsigned(text)?).ok(),
	}
}

fn parse_float(text: &[u8]) -> Option<f64> {
	std::str::from_utf8(text).ok()?.parse().ok()
}

fn parse_flags(line: &Line<'_>, text: &[u8]) -> Result<Flags, SamError> {
	parse_unsigned(text)
		.and_then(|value| u16::try_from(value).ok())
		.map(Flags::new)
		.ok_or_else(|| line.invalid(SamField::Flag))
}

fn parse_reference(
	line: &Line<'_>,
	header: &Header,
	text: &[u8],
	field: SamField,
) -> Result<Option<ContigId>, SamError> {
	if text == MISSING {
		return Ok(None);
	}
	let name = std::str::from_utf8(text).map_err(|_| line.invalid(field))?;
	header
		.contig_id(name)
		.map(Some)
		.ok_or_else(|| SamError::UnknownReference {
			line: line.number,
			name: name.to_owned(),
		})
}

fn parse_position(line: &Line<'_>, text: &[u8], field: SamField) -> Result<Option<u32>, SamError> {
	let value = parse_unsigned(text)
		.filter(|&value| value <= MAX_POSITION)
		.ok_or_else(|| line.invalid(field))?;
	match value.checked_sub(1) {
		Some(zero_based) => u32::try_from(zero_based)
			.map(Some)
			.map_err(|_| line.invalid(field)),
		None => Ok(None),
	}
}

fn parse_read(
	line: &Line<'_>,
	bases: &[u8],
	qualities: &[u8],
	record: &mut Record,
) -> Result<(), SamError> {
	let bases = if bases == MISSING { &[][..] } else { bases };
	let qualities = if qualities == MISSING {
		QualityInput::Missing
	} else {
		QualityInput::Ascii(qualities)
	};
	record
		.set_sequence(bases, qualities)
		.map_err(|_| line.invalid(SamField::Qualities))
}

fn check_cigar_length(line: &Line<'_>, record: &Record) -> Result<(), SamError> {
	let view = record.view();
	let cigar = view.stored_cigar();
	let consistent = cigar.is_empty()
		|| view.sequence().is_empty()
		|| cigar.query_len() == view.sequence().len() as u64;
	consistent
		.then_some(())
		.ok_or_else(|| line.invalid(SamField::Cigar))
}

fn parse_tags(line: &mut Line<'_>, record: &mut Record) -> Result<(), Error> {
	while let Some(token) = line.next_field() {
		parse_tag(token, record).ok_or_else(|| line.invalid(SamField::Tag))?;
	}
	Ok(())
}

fn parse_tag(token: &[u8], record: &mut Record) -> Option<()> {
	let [first, second, b':', kind, b':', value @ ..] = token else {
		return None;
	};
	let tag = [*first, *second];
	match kind {
		b'A' => match value {
			[character] => record.push_tag(tag, Value::Char(*character)).ok(),
			_ => None,
		},
		b'i' | b'I' => record.push_integer(tag, parse_signed(value)?).ok(),
		b'f' => record
			.push_tag(tag, Value::Float(to_single(parse_float(value)?)))
			.ok(),
		b'd' => record
			.push_tag(tag, Value::Double(parse_float(value)?))
			.ok(),
		b'Z' => record.push_tag(tag, Value::String(value)).ok(),
		b'H' if value.len() % 2 == 0 => record.push_tag(tag, Value::Hex(value)).ok(),
		b'B' => parse_array(tag, value, record),
		_ => None,
	}
}

#[allow(clippy::cast_possible_truncation)]
fn to_single(value: f64) -> f32 {
	value as f32
}

fn parse_array(tag: [u8; 2], value: &[u8], record: &mut Record) -> Option<()> {
	let (code, items) = value.split_first()?;
	match code {
		b'c' => push_array::<i8>(record, tag, items, |item| {
			i8::try_from(parse_signed(item)?).ok()
		}),
		b'C' => push_array::<u8>(record, tag, items, |item| {
			u8::try_from(parse_unsigned(item)?).ok()
		}),
		b's' => push_array::<i16>(record, tag, items, |item| {
			i16::try_from(parse_signed(item)?).ok()
		}),
		b'S' => push_array::<u16>(record, tag, items, |item| {
			u16::try_from(parse_unsigned(item)?).ok()
		}),
		b'i' => push_array::<i32>(record, tag, items, |item| {
			i32::try_from(parse_signed(item)?).ok()
		}),
		b'I' => push_array::<u32>(record, tag, items, |item| {
			u32::try_from(parse_unsigned(item)?).ok()
		}),
		b'f' => push_array::<f32>(record, tag, items, |item| parse_float(item).map(to_single)),
		_ => None,
	}
}

fn push_array<T: Element>(
	record: &mut Record,
	tag: [u8; 2],
	items: &[u8],
	parse: impl Fn(&[u8]) -> Option<T>,
) -> Option<()> {
	let mut writer = record.begin_array::<T>(tag);
	if items.is_empty() {
		return Some(());
	}
	for item in items.strip_prefix(b",")?.split(|&byte| byte == b',') {
		writer.push(parse(item)?).ok()?;
	}
	Some(())
}
