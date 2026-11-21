use std::io::Write;

use crate::cram::Error;
use crate::cram::decode::advance;

const CHUNK: usize = 1 << 16;

fn malformed() -> Error {
	Error::Malformed("record")
}

pub(crate) fn describe(
	source: &mut impl FnMut(u64, &mut [u8]) -> Result<(), Error>,
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
				let mut done = 0;
				while done < length {
					let step = (length - done).min(CHUNK);
					reference_bases.resize(step, 0);
					source(reference, reference_bases)?;
					let aligned = bases.get(read..read + step).ok_or_else(malformed)?;
					for (base, expected) in aligned.iter().zip(reference_bases.iter()) {
						if *base == b'='
							|| (base.eq_ignore_ascii_case(expected) && *expected != b'N')
						{
							run += 1;
						} else {
							write!(md, "{run}")?;
							md.push(*expected);
							run = 0;
							distance += 1;
						}
					}
					read += step;
					reference = advance(reference, step)?;
					done += step;
				}
			}
			1 => {
				read += length;
				distance = distance.saturating_add(op >> 4);
			}
			2 => {
				write!(md, "{run}^")?;
				let mut done = 0;
				while done < length {
					let step = (length - done).min(CHUNK);
					reference_bases.resize(step, 0);
					source(reference, reference_bases)?;
					md.extend_from_slice(reference_bases);
					reference = advance(reference, step)?;
					done += step;
				}
				run = 0;
				distance = distance.saturating_add(op >> 4);
			}
			3 => reference = advance(reference, length)?,
			4 => read += length,
			_ => {}
		}
	}
	write!(md, "{run}")?;
	Ok(distance)
}
