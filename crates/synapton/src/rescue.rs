const KMER_LENGTH: usize = 11;
const KMER_MASK: u32 = (1 << (2 * KMER_LENGTH)) - 1;
const MIN_VOTES: usize = 3;
const VOTE_WIDTH: i32 = 16;

#[derive(Debug, Default)]
pub(crate) struct Voter {
	table: Vec<(u32, u32)>,
	diagonals: Vec<i32>,
}

impl Voter {
	pub(crate) fn best_diagonal(&mut self, mate: &[u8], window: &[u8]) -> Option<i32> {
		self.table.clear();
		self.diagonals.clear();
		let mut kmer = 0u32;
		let mut valid = 0usize;
		for (position, &code) in mate.iter().enumerate() {
			if code > 3 {
				valid = 0;
				continue;
			}
			kmer = ((kmer << 2) | u32::from(code)) & KMER_MASK;
			valid += 1;
			if valid >= KMER_LENGTH {
				self.table
					.push((kmer, u32::try_from(position + 1 - KMER_LENGTH).ok()?));
			}
		}
		self.table.sort_unstable();
		kmer = 0;
		valid = 0;
		for (position, &code) in window.iter().enumerate() {
			if code > 3 {
				valid = 0;
				continue;
			}
			kmer = ((kmer << 2) | u32::from(code)) & KMER_MASK;
			valid += 1;
			if valid < KMER_LENGTH {
				continue;
			}
			let first = self.table.partition_point(|&(found, _)| found < kmer);
			let window_position = i32::try_from(position + 1 - KMER_LENGTH).ok()?;
			for &(found, mate_position) in self.table.get(first..).unwrap_or_default() {
				if found != kmer {
					break;
				}
				self.diagonals
					.push(window_position - i32::try_from(mate_position).ok()?);
			}
		}
		self.diagonals.sort_unstable();
		let mut best: Option<(usize, i32)> = None;
		let mut tail = 0usize;
		for (head, &diagonal) in self.diagonals.iter().enumerate() {
			while self
				.diagonals
				.get(tail)
				.is_some_and(|&earlier| diagonal - earlier > VOTE_WIDTH)
			{
				tail += 1;
			}
			let votes = head + 1 - tail;
			if best.is_none_or(|(top, _)| votes > top) {
				let middle = self
					.diagonals
					.get(tail + votes / 2)
					.copied()
					.unwrap_or(diagonal);
				best = Some((votes, middle));
			}
		}
		best.filter(|&(votes, _)| votes >= MIN_VOTES)
			.map(|(_, diagonal)| diagonal)
	}
}
