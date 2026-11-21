const LOAD_PERCENT: usize = 77;
const MIN_BUCKETS: usize = 4;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Slot {
	Empty,
	Occupied,
	Moved,
}

#[derive(Debug, Default)]
struct Table {
	keys: Vec<u32>,
	occupied: Vec<bool>,
	size: usize,
}

fn upper_bound(buckets: usize) -> usize {
	(buckets * LOAD_PERCENT + 50) / 100
}

fn find_bucket(keys: &[u32], occupied: &[bool], key: u32) -> (usize, bool) {
	let mask = occupied.len() - 1;
	let mut bucket = usize::try_from(key).unwrap_or(usize::MAX) & mask;
	let mut step = 0;
	while occupied.get(bucket).copied().unwrap_or(false) {
		if keys.get(bucket) == Some(&key) {
			return (bucket, true);
		}
		step += 1;
		bucket = (bucket + step) & mask;
	}
	(bucket, false)
}

fn free_bucket(occupied: &[bool], key: u32) -> usize {
	find_bucket(&[], occupied, key).0
}

impl Table {
	fn put(&mut self, key: u32) {
		if self.size >= upper_bound(self.keys.len()) {
			self.grow();
		}
		let (bucket, present) = find_bucket(&self.keys, &self.occupied, key);
		if present {
			return;
		}
		if let (Some(taken), Some(stored)) =
			(self.occupied.get_mut(bucket), self.keys.get_mut(bucket))
		{
			*taken = true;
			*stored = key;
			self.size += 1;
		}
	}

	fn grow(&mut self) {
		let old_len = self.keys.len();
		let new_len = (old_len + 1).next_power_of_two().max(MIN_BUCKETS);
		let mut old: Vec<Slot> = self
			.occupied
			.iter()
			.map(|&taken| if taken { Slot::Occupied } else { Slot::Empty })
			.collect();
		let mut placed = vec![false; new_len];
		self.keys.resize(new_len, 0);
		for index in 0..old_len {
			if old.get(index) != Some(&Slot::Occupied) {
				continue;
			}
			let mut key = self.keys.get(index).copied().unwrap_or_default();
			if let Some(slot) = old.get_mut(index) {
				*slot = Slot::Moved;
			}
			loop {
				let target = free_bucket(&placed, key);
				if let Some(taken) = placed.get_mut(target) {
					*taken = true;
				}
				let displaced = old.get(target) == Some(&Slot::Occupied);
				let Some(stored) = self.keys.get_mut(target) else {
					break;
				};
				if !displaced {
					*stored = key;
					break;
				}
				std::mem::swap(stored, &mut key);
				if let Some(slot) = old.get_mut(target) {
					*slot = Slot::Moved;
				}
			}
		}
		self.occupied = placed;
	}

	fn in_bucket_order(&self) -> Vec<u32> {
		self.keys
			.iter()
			.zip(&self.occupied)
			.filter(|(_, taken)| **taken)
			.map(|(key, _)| *key)
			.collect()
	}
}

pub(crate) fn bucket_order(puts: impl IntoIterator<Item = u32>) -> Vec<u32> {
	let mut table = Table::default();
	for key in puts {
		table.put(key);
	}
	table.in_bucket_order()
}
