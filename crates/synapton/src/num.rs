pub(crate) fn u32_of(value: usize) -> u32 {
	u32::try_from(value).unwrap_or(u32::MAX)
}

pub(crate) fn u16_of(value: usize) -> u16 {
	u16::try_from(value).unwrap_or(u16::MAX)
}

pub(crate) fn u32_from(value: u64) -> u32 {
	u32::try_from(value).unwrap_or(u32::MAX)
}

pub(crate) fn len_u64(value: usize) -> u64 {
	u64::try_from(value).unwrap_or(u64::MAX)
}

pub(crate) fn usize_of(value: u64) -> usize {
	usize::try_from(value).unwrap_or(usize::MAX)
}

pub(crate) fn i64_of(value: u64) -> i64 {
	i64::try_from(value).unwrap_or(i64::MAX)
}

pub(crate) fn i64_len(value: usize) -> i64 {
	i64::try_from(value).unwrap_or(i64::MAX)
}

pub(crate) fn u64_of(value: i64) -> u64 {
	u64::try_from(value).unwrap_or(0)
}

pub(crate) fn i32_of(value: i64) -> i32 {
	i32::try_from(value).unwrap_or(if value < 0 { i32::MIN } else { i32::MAX })
}

pub(crate) fn float_of(value: i64) -> f64 {
	f64::from(i32_of(value))
}

pub(crate) fn float_len(value: usize) -> f64 {
	f64::from(u32_of(value))
}

#[expect(
	clippy::cast_possible_truncation,
	reason = "the value is clamped to the i32 range first"
)]
pub(crate) fn round_i32(value: f64) -> i32 {
	value.clamp(f64::from(i32::MIN), f64::from(i32::MAX)) as i32
}

pub(crate) fn round_i64(value: f64) -> i64 {
	i64::from(round_i32(value))
}
