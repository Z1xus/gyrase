const ROW: usize = 8;

#[derive(Debug)]
pub(crate) struct FrontBatch {
	table: Vec<u32>,
	pool: Vec<u8>,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct ReadFront<'a> {
	pub(crate) anchors: [&'a [u8]; 2],
	pub(crate) skipped: [u32; 2],
	pub(crate) total: [u32; 2],
	pub(crate) quick: &'a [u8],
}

impl FrontBatch {
	#[cfg(feature = "gpu")]
	pub(crate) fn new(table: Vec<u32>, pool: Vec<u8>) -> Self {
		Self { table, pool }
	}

	pub(crate) fn read(&self, index: usize) -> Option<ReadFront<'_>> {
		let row: [u32; ROW] = self
			.table
			.get(index * ROW..(index + 1) * ROW)?
			.try_into()
			.ok()?;
		let [
			first_start,
			first_count,
			second_start,
			second_count,
			quick_start,
			quick_count,
			stats,
			flags,
		] = row;
		if flags != 0 {
			return None;
		}
		let words = |start: u32, count: u32| {
			let start = start as usize * 8;
			self.pool.get(start..start + count as usize * 8)
		};
		let counts = |lane: u32| {
			[
				(stats >> (16 * lane)) & 0xff,
				(stats >> (16 * lane + 8)) & 0xff,
			]
		};
		let [first_skipped, first_total] = counts(0);
		let [second_skipped, second_total] = counts(1);
		Some(ReadFront {
			anchors: [
				words(first_start, first_count)?,
				words(second_start, second_count)?,
			],
			skipped: [first_skipped, second_skipped],
			total: [first_total, second_total],
			quick: words(quick_start, quick_count)?,
		})
	}
}

pub(crate) fn words(bytes: &[u8]) -> impl Iterator<Item = u64> + '_ {
	bytes
		.as_chunks::<8>()
		.0
		.iter()
		.map(|&word| u64::from_le_bytes(word))
}

pub(crate) fn skipped_fraction(skipped: u32, total: u32) -> f32 {
	if total == 0 {
		0.0
	} else {
		let total = u16::try_from(total).unwrap_or(u16::MAX);
		f32::from(u16::try_from(skipped).unwrap_or(u16::MAX)) / f32::from(total)
	}
}
