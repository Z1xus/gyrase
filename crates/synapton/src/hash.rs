use rayon::prelude::*;

const MIX_FIRST: u64 = 0xbf58_476d_1ce4_e5b9;
const MIX_SECOND: u64 = 0x94d0_49bb_1331_11eb;
const FOLD: u64 = 0x9e37_79b1_85eb_ca87;
const CHECKSUM_CHUNK_WORDS: usize = 1 << 22;

pub(crate) const fn mix(value: u64) -> u64 {
	let mut x = value ^ (value >> 30);
	x = x.wrapping_mul(MIX_FIRST);
	x ^= x >> 27;
	x = x.wrapping_mul(MIX_SECOND);
	x ^ (x >> 31)
}

pub(crate) fn checksum_words(words: &[u64]) -> u64 {
	let chunks: Vec<u64> = words
		.par_chunks(CHECKSUM_CHUNK_WORDS)
		.map(checksum_chunk)
		.collect();
	checksum_chunk(&chunks)
}

pub(crate) fn checksum_chunk(words: &[u64]) -> u64 {
	let mut lanes = [FOLD; 4];
	let (groups, remainder) = words.as_chunks::<4>();
	for group in groups {
		for (lane, &word) in lanes.iter_mut().zip(group) {
			*lane = (lane.rotate_left(23) ^ word).wrapping_mul(MIX_FIRST);
		}
	}
	let mut folded = words.len() as u64;
	for &word in remainder {
		folded = (folded.rotate_left(23) ^ word).wrapping_mul(MIX_FIRST);
	}
	for lane in lanes {
		folded = mix(folded ^ lane);
	}
	folded
}
