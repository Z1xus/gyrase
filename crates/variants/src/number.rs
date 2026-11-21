const SMALL_MIN: f64 = 0.0001;
const SMALL_MAX: f64 = 999_999.0;
const SIGNIFICANT_DIGITS: i32 = 6;
const PLACES_BY_MAGNITUDE: [(f64, i32); 9] = [
	(0.001, 9),
	(0.01, 8),
	(0.1, 7),
	(1.0, 6),
	(10.0, 5),
	(100.0, 4),
	(1000.0, 3),
	(10000.0, 2),
	(100_000.0, 1),
];

pub(crate) fn push_uint(out: &mut Vec<u8>, mut value: u64) {
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

pub(crate) fn push_int(out: &mut Vec<u8>, value: i64) {
	if value < 0 {
		out.push(b'-');
	}
	push_uint(out, value.unsigned_abs());
}

pub(crate) fn push_float(out: &mut Vec<u8>, value: f32) {
	let mut number = f64::from(value);
	if number == 0.0 {
		out.extend_from_slice(if number.is_sign_negative() {
			b"-0"
		} else {
			b"0"
		});
		return;
	}
	if number.is_nan() {
		out.extend_from_slice(if number.is_sign_negative() {
			b"-nan"
		} else {
			b"nan"
		});
		return;
	}
	if number < 0.0 {
		out.push(b'-');
		number = -number;
	}
	if number.is_infinite() {
		out.extend_from_slice(b"inf");
	} else if (SMALL_MIN..=SMALL_MAX).contains(&number) {
		push_fixed(out, number);
	} else {
		push_general(out, number);
	}
}

fn push_fixed(out: &mut Vec<u8>, number: f64) {
	let places = PLACES_BY_MAGNITUDE
		.iter()
		.find(|(limit, _)| number < *limit)
		.map_or(0, |(_, places)| *places);
	let scaled = (number * 10f64.powi(places)).round_ties_even();
	let mut text = Vec::with_capacity(16);
	push_uint(&mut text, whole_number(scaled));
	let places = usize::try_from(places).unwrap_or(0);
	while text.len() <= places {
		text.insert(0, b'0');
	}
	let point = text.len() - places;
	let (whole, fraction) = text.split_at(point);
	out.extend_from_slice(whole);
	let kept = fraction
		.iter()
		.rposition(|digit| *digit != b'0')
		.map_or(0, |last| last + 1);
	if kept > 0 {
		out.push(b'.');
		out.extend_from_slice(fraction.get(..kept).unwrap_or_default());
	}
}

fn push_general(out: &mut Vec<u8>, number: f64) {
	let scientific = format!("{:.*e}", (SIGNIFICANT_DIGITS - 1) as usize, number);
	let Some((mantissa, exponent)) = scientific.split_once('e') else {
		return;
	};
	let exponent: i32 = exponent.parse().unwrap_or(0);
	if (-4..SIGNIFICANT_DIGITS).contains(&exponent) {
		let places = usize::try_from(SIGNIFICANT_DIGITS - 1 - exponent).unwrap_or(0);
		let fixed = format!("{number:.places$}");
		out.extend_from_slice(trim_fraction(&fixed).as_bytes());
		return;
	}
	out.extend_from_slice(trim_fraction(mantissa).as_bytes());
	out.push(b'e');
	out.push(if exponent < 0 { b'-' } else { b'+' });
	let magnitude = exponent.unsigned_abs();
	if magnitude < 10 {
		out.push(b'0');
	}
	push_uint(out, u64::from(magnitude));
}

fn trim_fraction(text: &str) -> &str {
	if text.contains('.') {
		text.trim_end_matches('0').trim_end_matches('.')
	} else {
		text
	}
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub(crate) const fn whole_number(value: f64) -> u64 {
	value as u64
}

#[allow(clippy::cast_possible_truncation)]
pub(crate) const fn narrow_float(value: f64) -> f32 {
	value as f32
}
