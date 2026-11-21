use crate::config::Scoring;
use crate::num::usize_of;
use crate::reference::{ContigSpan, PackedReference};

const LOW_BITS: u64 = 0x5555_5555_5555_5555;
const MAX_MISMATCHES: u32 = 2;
const END_MARGIN: usize = 10;
const WORD: usize = 32;

pub(crate) fn pack(codes: &[u8], out: &mut Vec<u64>) -> bool {
	let mut clean = true;
	for chunk in codes.chunks(WORD) {
		let mut word = 0u64;
		for (index, &code) in chunk.iter().enumerate() {
			clean &= code < 4;
			word |= u64::from(code & 3) << (2 * index);
		}
		out.push(word);
	}
	clean
}

fn mismatch_mask(
	reference: &PackedReference,
	position: u64,
	read: &[u64],
	length: usize,
	index: usize,
) -> u64 {
	let bases = (length - index * WORD).min(WORD);
	let word = read.get(index).copied().unwrap_or(0);
	let diff = reference.word_at(position + (index * WORD) as u64) ^ word;
	let mask = (diff | diff >> 1) & LOW_BITS;
	if bases < WORD {
		mask & ((1u64 << (2 * bases)) - 1)
	} else {
		mask
	}
}

fn count_mismatches(
	reference: &PackedReference,
	position: u64,
	read: &[u64],
	length: usize,
	limit: u32,
) -> Option<u32> {
	let mut total = 0;
	for index in 0..read.len() {
		total += mismatch_mask(reference, position, read, length, index).count_ones();
		if total > limit {
			return None;
		}
	}
	Some(total)
}

fn interior(reference: &PackedReference, position: u64, read: &[u64], length: usize) -> bool {
	(0..read.len()).all(|index| {
		let mut mask = mismatch_mask(reference, position, read, length, index);
		while mask != 0 {
			let offset = index * WORD + usize_of(u64::from(mask.trailing_zeros() / 2));
			if offset < END_MARGIN || offset + END_MARGIN >= length {
				return false;
			}
			mask &= mask - 1;
		}
		true
	})
}

pub(crate) fn ungapped(
	reference: &PackedReference,
	scoring: &Scoring,
	span: &ContigSpan,
	flat: i64,
	read: &[u64],
	length: usize,
) -> Option<i32> {
	let position = u64::try_from(flat).ok()?;
	let end = position + length as u64;
	if position < span.start || end > span.end() || reference.has_ambiguous(position, end) {
		return None;
	}
	let mismatches = count_mismatches(reference, position, read, length, MAX_MISMATCHES)?;
	if mismatches == MAX_MISMATCHES && !interior(reference, position, read, length) {
		return None;
	}
	let matched = i32::from(scoring.matched);
	let penalty = matched + i32::from(scoring.mismatch);
	let length = i32::try_from(length).ok()?;
	Some(length * matched - i32::try_from(mismatches).ok()? * penalty)
}

pub(crate) fn scan(
	reference: &PackedReference,
	scoring: &Scoring,
	span: &ContigSpan,
	range: (i64, i64),
	read: &[u64],
	length: usize,
) -> Option<(i64, i32)> {
	let first = read.first().copied()?;
	let mut found = None;
	for flat in range.0..=range.1 {
		let Ok(position) = u64::try_from(flat) else {
			continue;
		};
		let diff = reference.word_at(position) ^ first;
		if ((diff | diff >> 1) & LOW_BITS).count_ones() > MAX_MISMATCHES {
			continue;
		}
		if let Some(score) = ungapped(reference, scoring, span, flat, read, length) {
			if found.is_some() {
				return None;
			}
			found = Some((flat, score));
		}
	}
	found
}
