use memchr::memchr_iter;

use crate::Error;
use crate::header::{Header, Version};
use crate::value::{Cells, END_INT, MISSING_INT, Store, ValueType, end_float, missing_float};

use super::text::{parse_float_item, parse_int_item};

const PREFIXED_PHASING: Version = Version { major: 4, minor: 4 };
const MAX_ALLELE_INDEX: u32 = (i32::MAX as u32 >> 1) - 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ColumnKind {
	Genotype,
	Ints,
	Floats,
	Text,
}

pub(crate) fn kind_of(header: &Header, key: u32) -> Option<ColumnKind> {
	let name = header.name(key)?;
	Some(match header.type_of_format(key)? {
		ValueType::Integer => ColumnKind::Ints,
		ValueType::Float => ColumnKind::Floats,
		ValueType::String if name == "GT" => ColumnKind::Genotype,
		ValueType::String | ValueType::Character => ColumnKind::Text,
		ValueType::Flag => return None,
	})
}

fn subfields(samples: &[u8], index: usize) -> impl Iterator<Item = Option<&[u8]>> {
	samples
		.split(|byte| *byte == b'\t')
		.map(move |sample| sample.split(|byte| *byte == b':').nth(index))
}

pub(crate) fn load_column(
	samples: &[u8],
	index: usize,
	kind: ColumnKind,
	header: &Header,
	store: &mut Store,
) -> Result<(Cells, usize), Error> {
	let width = column_width(samples, index, kind);
	let version = header.version();
	match kind {
		ColumnKind::Ints => {
			let start = store.ints.len();
			for item in subfields(samples, index) {
				fill_numbers(
					&mut store.ints,
					item,
					width,
					parse_int_item,
					MISSING_INT,
					END_INT,
				);
			}
			Ok((Cells::Ints(start..store.ints.len()), width))
		}
		ColumnKind::Floats => {
			let start = store.floats.len();
			for item in subfields(samples, index) {
				fill_numbers(
					&mut store.floats,
					item,
					width,
					parse_float_item,
					missing_float(),
					end_float(),
				);
			}
			Ok((Cells::Floats(start..store.floats.len()), width))
		}
		ColumnKind::Genotype => {
			let start = store.ints.len();
			for item in subfields(samples, index) {
				fill_genotype(&mut store.ints, item, width, version)?;
			}
			Ok((Cells::Ints(start..store.ints.len()), width))
		}
		ColumnKind::Text => {
			let start = store.text.len();
			for item in subfields(samples, index) {
				fill_text(&mut store.text, item, width);
			}
			Ok((Cells::Text(start..store.text.len()), width))
		}
	}
}

fn column_width(samples: &[u8], index: usize, kind: ColumnKind) -> usize {
	let widest = subfields(samples, index)
		.map(|item| {
			let item = item.unwrap_or_default();
			match kind {
				ColumnKind::Ints | ColumnKind::Floats => memchr_iter(b',', item).count() + 1,
				ColumnKind::Genotype => {
					memchr_iter(b'/', item).count() + memchr_iter(b'|', item).count() + 1
				}
				ColumnKind::Text => item.len(),
			}
		})
		.max()
		.unwrap_or(0);
	widest.max(1)
}

fn fill_numbers<T: Copy>(
	out: &mut Vec<T>,
	item: Option<&[u8]>,
	width: usize,
	parse: impl Fn(&[u8]) -> T,
	missing: T,
	end: T,
) {
	let start = out.len();
	match item {
		Some(item) => out.extend(item.split(|byte| *byte == b',').map(parse)),
		None => out.push(missing),
	}
	let written = out.len() - start;
	for _ in written..width {
		out.push(end);
	}
}

fn fill_text(out: &mut Vec<u8>, item: Option<&[u8]>, width: usize) {
	let start = out.len();
	match item {
		Some(item) => out.extend_from_slice(item),
		None => out.push(b'.'),
	}
	let written = out.len() - start;
	out.resize(start + written.max(width), 0);
}

fn fill_genotype(
	out: &mut Vec<i32>,
	item: Option<&[u8]>,
	width: usize,
	version: Version,
) -> Result<(), Error> {
	let start = out.len();
	match item {
		Some(item) => parse_genotype(item, out, version)?,
		None => out.push(0),
	}
	let written = out.len() - start;
	for _ in written..width {
		out.push(END_INT);
	}
	Ok(())
}

pub(crate) fn parse_genotype(
	text: &[u8],
	out: &mut Vec<i32>,
	version: Version,
) -> Result<(), Error> {
	let first = out.len();
	let mut rest = text;
	let mut phased = false;
	let mut prefixed = false;
	if version >= PREFIXED_PHASING
		&& let [separator @ (b'|' | b'/'), tail @ ..] = rest
	{
		phased = *separator == b'|';
		prefixed = true;
		rest = tail;
	}
	let mut ploidy = 0u32;
	let mut unknown_first = false;
	let mut any_unphased = false;
	loop {
		ploidy += 1;
		let code = if let [b'.', tail @ ..] = rest {
			rest = tail;
			unknown_first |= ploidy == 1;
			i32::from(phased)
		} else {
			let digits = rest.iter().take_while(|byte| byte.is_ascii_digit()).count();
			let number = rest
				.get(..digits)
				.filter(|slice| !slice.is_empty())
				.and_then(|slice| std::str::from_utf8(slice).ok())
				.and_then(|slice| slice.parse::<u32>().ok())
				.filter(|&value| value <= MAX_ALLELE_INDEX)
				.ok_or(Error::Vcf("invalid genotype"))?;
			rest = rest.get(digits..).unwrap_or_default();
			i32::try_from((number + 1) << 1).map_err(|_| Error::Vcf("genotype is too large"))?
				| i32::from(phased)
		};
		out.push(code);
		any_unphased |= ploidy != 1 && !phased;
		match rest {
			[separator @ (b'|' | b'/'), tail @ ..] => {
				phased = *separator == b'|';
				rest = tail;
			}
			[] => break,
			_ => return Err(Error::Vcf("invalid genotype")),
		}
	}
	if !prefixed {
		let mark = if ploidy == 1 {
			!unknown_first
		} else {
			!any_unphased
		};
		if mark && let Some(slot) = out.get_mut(first) {
			*slot |= 1;
		}
	}
	Ok(())
}
