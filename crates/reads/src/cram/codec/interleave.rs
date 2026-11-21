pub(super) fn deinterleave<const N: usize>(interleaved: &[u8], output: &mut [u8]) {
	let lane_len = interleaved.len() / N;
	if lane_len == 0 {
		return;
	}
	for (lane, out_lane) in output.chunks_exact_mut(lane_len).take(N).enumerate() {
		for (slot, group) in out_lane.iter_mut().zip(interleaved.as_chunks::<N>().0) {
			*slot = group.get(lane).copied().unwrap_or(0);
		}
	}
}
