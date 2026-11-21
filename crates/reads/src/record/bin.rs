const LEVELS: [(u32, i64); 5] = [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)];

pub(crate) fn region_bin(start: i64, end: i64) -> u16 {
	let last = end - 1;
	LEVELS
		.iter()
		.find(|(shift, _)| start >> shift == last >> shift)
		.and_then(|(shift, offset)| u16::try_from(offset + (start >> shift)).ok())
		.unwrap_or(0)
}
