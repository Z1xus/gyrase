use std::fmt;
use std::io::{Cursor, Write};

const FIXED_DECIMALS: [(f64, f64, usize); 9] = [
	(0.001, 1e9, 9),
	(0.01, 1e8, 8),
	(0.1, 1e7, 7),
	(1.0, 1e6, 6),
	(10.0, 1e5, 5),
	(100.0, 1e4, 4),
	(1000.0, 1e3, 3),
	(10000.0, 1e2, 2),
	(100_000.0, 1e1, 1),
];
const SHORT_FORM_MIN: f64 = 0.0001;
const SHORT_FORM_MAX: f64 = 999_999.0;
const SIGNIFICANT_DIGITS: i32 = 6;
const SCRATCH_LEN: usize = 64;

pub(crate) fn push_unsigned(out: &mut Vec<u8>, mut value: u64) {
	let mut digits = [0u8; 20];
	let mut start = digits.len();
	loop {
		start -= 1;
		if let Some(slot) = digits.get_mut(start) {
			*slot = b'0' + u8::try_from(value % 10).unwrap_or(0);
		}
		value /= 10;
		if value == 0 {
			break;
		}
	}
	out.extend_from_slice(digits.get(start..).unwrap_or_default());
}

pub(crate) fn push_signed(out: &mut Vec<u8>, value: i64) {
	if value < 0 {
		out.push(b'-');
	}
	push_unsigned(out, value.unsigned_abs());
}

pub(crate) fn push_general(out: &mut Vec<u8>, value: f64) {
	if value.is_sign_negative() {
		out.push(b'-');
	}
	push_general_magnitude(out, value.abs());
}

fn push_general_magnitude(out: &mut Vec<u8>, magnitude: f64) {
	if magnitude.is_nan() {
		out.extend_from_slice(b"nan");
		return;
	}
	if magnitude.is_infinite() {
		out.extend_from_slice(b"inf");
		return;
	}
	if magnitude == 0.0 {
		out.push(b'0');
		return;
	}
	let mut scratch = [0u8; SCRATCH_LEN];
	let scientific = render(&mut scratch, format_args!("{magnitude:.5e}"));
	let (mantissa, exponent) = split_exponent(scientific);
	if (-4..SIGNIFICANT_DIGITS).contains(&exponent) {
		let decimals = usize::try_from(SIGNIFICANT_DIGITS - 1 - exponent).unwrap_or(0);
		let mut fixed_scratch = [0u8; SCRATCH_LEN];
		let fixed = render(&mut fixed_scratch, format_args!("{magnitude:.decimals$}"));
		push_trimmed(out, fixed);
	} else {
		push_trimmed(out, mantissa);
		out.push(b'e');
		out.push(if exponent < 0 { b'-' } else { b'+' });
		if exponent.unsigned_abs() < 10 {
			out.push(b'0');
		}
		push_unsigned(out, u64::from(exponent.unsigned_abs()));
	}
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub(crate) fn push_short(out: &mut Vec<u8>, value: f64) {
	if value == 0.0 {
		out.extend_from_slice(if value.is_sign_negative() {
			b"-0"
		} else {
			b"0"
		});
		return;
	}
	let magnitude = value.abs();
	if value < 0.0 {
		out.push(b'-');
	}
	if !(SHORT_FORM_MIN..=SHORT_FORM_MAX).contains(&magnitude) {
		push_general_magnitude(out, magnitude);
		return;
	}
	let (scale, decimals) = FIXED_DECIMALS
		.iter()
		.find(|(limit, _, _)| magnitude < *limit)
		.map_or((1.0, 0), |&(_, scale, decimals)| (scale, decimals));
	let scaled = (magnitude * scale).round_ties_even() as u64;
	let divisor = 10u64.pow(u32::try_from(decimals).unwrap_or(0));
	push_unsigned(out, scaled / divisor);
	let mut fraction = scaled % divisor;
	if fraction == 0 {
		return;
	}
	let mut digits = [b'0'; 9];
	for slot in digits.iter_mut().take(decimals).rev() {
		*slot = b'0' + u8::try_from(fraction % 10).unwrap_or(0);
		fraction /= 10;
	}
	out.push(b'.');
	out.extend_from_slice(trim_zeros(digits.get(..decimals).unwrap_or_default()));
}

fn render<'a>(scratch: &'a mut [u8; SCRATCH_LEN], arguments: fmt::Arguments<'_>) -> &'a [u8] {
	let mut cursor = Cursor::new(&mut scratch[..]);
	let _ = cursor.write_fmt(arguments);
	let written = usize::try_from(cursor.position()).unwrap_or(0);
	scratch.get(..written).unwrap_or_default()
}

fn split_exponent(scientific: &[u8]) -> (&[u8], i32) {
	let position = scientific.iter().position(|&byte| byte == b'e');
	let (mantissa, exponent) = match position {
		Some(at) => scientific.split_at(at),
		None => (scientific, &[][..]),
	};
	let exponent = exponent
		.get(1..)
		.and_then(|digits| std::str::from_utf8(digits).ok())
		.and_then(|digits| digits.parse().ok())
		.unwrap_or(0);
	(mantissa, exponent)
}

fn trim_zeros(fraction: &[u8]) -> &[u8] {
	let end = fraction
		.iter()
		.rposition(|&byte| byte != b'0')
		.map_or(0, |last| last + 1);
	fraction.get(..end).unwrap_or_default()
}

fn push_trimmed(out: &mut Vec<u8>, number: &[u8]) {
	match number.iter().position(|&byte| byte == b'.') {
		Some(dot) => {
			let (whole, fraction) = number.split_at(dot);
			out.extend_from_slice(whole);
			let fraction = trim_zeros(fraction.get(1..).unwrap_or_default());
			if !fraction.is_empty() {
				out.push(b'.');
				out.extend_from_slice(fraction);
			}
		}
		None => out.extend_from_slice(number),
	}
}
