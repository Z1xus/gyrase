use gyrase_genome::Base;

use crate::Observation;

const DISCORDANT_MATE_KEEP_TENTHS: u16 = 8;

#[derive(Debug, Clone, Copy)]
pub struct Column<'a> {
	pub position: u64,
	pub reference: Base,
	pub observations: &'a [Observation],
}

pub(crate) fn merge_mates(observations: &mut [Observation]) -> usize {
	let mut unpaired = 0;
	for index in 0..observations.len() {
		if observations[index].fragment.is_none() {
			observations.swap(unpaired, index);
			unpaired += 1;
		}
	}
	let tail = &mut observations[unpaired..];
	if tail.len() < 2 {
		return observations.len();
	}
	tail.sort_unstable_by_key(|observation| observation.fragment);
	let mut kept = unpaired;
	for index in unpaired..observations.len() {
		let current = observations[index];
		match observations[unpaired..kept].last_mut() {
			Some(previous) if previous.fragment == current.fragment => {
				*previous = combine(*previous, current);
			}
			_ => {
				observations[kept] = current;
				kept += 1;
			}
		}
	}
	kept
}

fn combine(first: Observation, second: Observation) -> Observation {
	let (stronger, weaker) = if first.quality >= second.quality {
		(first, second)
	} else {
		(second, first)
	};
	if stronger.base == weaker.base {
		stronger
	} else {
		let quality = u16::from(stronger.quality) * DISCORDANT_MATE_KEEP_TENTHS / 10;
		Observation {
			quality: u8::try_from(quality).unwrap_or(u8::MAX),
			..stronger
		}
	}
}
