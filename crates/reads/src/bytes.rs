pub(crate) fn u16_le(bytes: &[u8], at: usize) -> Option<u16> {
	Some(u16::from_le_bytes(*bytes.get(at..)?.first_chunk()?))
}

pub(crate) fn u32_le(bytes: &[u8], at: usize) -> Option<u32> {
	Some(u32::from_le_bytes(*bytes.get(at..)?.first_chunk()?))
}

pub(crate) fn i32_le(bytes: &[u8], at: usize) -> Option<i32> {
	Some(i32::from_le_bytes(*bytes.get(at..)?.first_chunk()?))
}
