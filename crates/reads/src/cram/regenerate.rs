use std::io::Write;

use crate::cram::Error;
use crate::cram::decode::advance;
use crate::cram::reference::Window;

fn malformed() -> Error {
	Error::Malformed("record")
}

pub(crate) fn describe(
	window: &Window,
	start: u64,
	operations: impl IntoIterator<Item = u32>,
	bases: &[u8],
	md: &mut Vec<u8>,
	reference_bases: &mut Vec<u8>,
) -> Result<u32, Error> {
	md.clear();
	let mut run = 0u32;
	let mut distance = 0u32;
	let mut read = 0usize;
	let mut reference = start;
	for op in operations {
		let length = (op >> 4) as usize;
		match op & 0xf {
			0 | 7 | 8 => {
				reference_bases.resize(length, 0);
				window.copy(reference, reference_bases)?;
				let aligned = bases.get(read..read + length).ok_or_else(malformed)?;
				if aligned == reference_bases.as_slice() && !reference_bases.contains(&b'N') {
					run = run.saturating_add(u32::try_from(length).unwrap_or(u32::MAX));
				} else {
					for (base, expected) in aligned.iter().zip(reference_bases.iter()) {
						if base.eq_ignore_ascii_case(expected) && *expected != b'N' {
							run += 1;
						} else {
							write!(md, "{run}")?;
							md.push(*expected);
							run = 0;
							distance += 1;
						}
					}
				}
				read += length;
				reference = advance(reference, length)?;
			}
			1 => {
				read += length;
				distance = distance.saturating_add(op >> 4);
			}
			2 => {
				reference_bases.resize(length, 0);
				window.copy(reference, reference_bases)?;
				write!(md, "{run}^")?;
				md.extend_from_slice(reference_bases);
				run = 0;
				distance = distance.saturating_add(op >> 4);
				reference = advance(reference, length)?;
			}
			3 => reference = advance(reference, length)?,
			4 => read += length,
			_ => {}
		}
	}
	write!(md, "{run}")?;
	Ok(distance)
}
