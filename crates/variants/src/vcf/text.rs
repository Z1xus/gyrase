use crate::Error;
use crate::number::narrow_float;
use crate::value::{Cells, MISSING_INT, Store, ValueType, missing_float};

const MIN_INT: i64 = -2_147_483_640;
const MAX_INT: i64 = 2_147_483_647;

pub(crate) fn parse_position(text: &[u8]) -> Result<u32, Error> {
	if text.is_empty() {
		return Err(Error::Vcf("empty position"));
	}
	let mut value: u32 = 0;
	for &byte in text {
		if !byte.is_ascii_digit() {
			return Err(Error::Vcf("position is not a number"));
		}
		value = value
			.checked_mul(10)
			.and_then(|shifted| shifted.checked_add(u32::from(byte - b'0')))
			.filter(|&value| i32::try_from(value).is_ok())
			.ok_or(Error::Vcf("position is too large"))?;
	}
	Ok(value)
}

pub(crate) fn parse_leading_uint(text: &[u8]) -> Option<u32> {
	let digits = text.iter().take_while(|byte| byte.is_ascii_digit()).count();
	std::str::from_utf8(text.get(..digits)?).ok()?.parse().ok()
}

pub(crate) fn parse_qual(text: &[u8]) -> Result<Option<f32>, Error> {
	if text == b"." {
		return Ok(None);
	}
	parse_f64(text)
		.map(|value| Some(narrow_float(value)))
		.ok_or(Error::Vcf("QUAL is not a number"))
}

fn parse_f64(text: &[u8]) -> Option<f64> {
	std::str::from_utf8(text).ok()?.parse().ok()
}

pub(crate) fn parse_int_item(item: &[u8]) -> i32 {
	let (negative, digits) = match item {
		[b'-', rest @ ..] => (true, rest),
		[b'+', rest @ ..] => (false, rest),
		other => (false, other),
	};
	let count = digits
		.iter()
		.take_while(|byte| byte.is_ascii_digit())
		.count();
	if count == 0 {
		return MISSING_INT;
	}
	let Some(text) = digits
		.get(..count)
		.and_then(|slice| std::str::from_utf8(slice).ok())
	else {
		return MISSING_INT;
	};
	let Ok(magnitude) = text.parse::<i64>() else {
		return MISSING_INT;
	};
	let value = if negative { -magnitude } else { magnitude };
	if (MIN_INT..=MAX_INT).contains(&value) {
		i32::try_from(value).unwrap_or(MISSING_INT)
	} else {
		MISSING_INT
	}
}

pub(crate) fn parse_float_item(item: &[u8]) -> f32 {
	if item == b"." {
		return missing_float();
	}
	parse_f64(item).map_or_else(missing_float, narrow_float)
}

pub(crate) fn parse_info_value(
	value_type: ValueType,
	text: &[u8],
	store: &mut Store,
) -> (Cells, usize) {
	match value_type {
		ValueType::Integer => {
			let cells = store.push_ints(text.split(|byte| *byte == b',').map(parse_int_item));
			(cells.clone(), cell_len(&cells))
		}
		ValueType::Float => {
			let cells = store.push_floats(text.split(|byte| *byte == b',').map(parse_float_item));
			(cells.clone(), cell_len(&cells))
		}
		ValueType::Flag | ValueType::Character | ValueType::String => {
			(store.push_text(text), text.len())
		}
	}
}

fn cell_len(cells: &Cells) -> usize {
	match cells {
		Cells::Flag => 0,
		Cells::Ints(range) | Cells::Floats(range) | Cells::Text(range) => range.len(),
	}
}
