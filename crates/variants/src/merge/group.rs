use std::cmp::Reverse;

use super::input::Taken;
use crate::Error;
use crate::record::Record;

#[derive(Debug, Clone, Copy)]
pub(super) struct Member {
	pub(super) input: usize,
	pub(super) slot: usize,
}

#[derive(Debug, Default)]
pub(super) struct Groups {
	members: Vec<Vec<Member>>,
	used: usize,
	order: Vec<usize>,
}

fn same_allele(left: &[u8], right: &[u8]) -> bool {
	if left.first() == Some(&b'<') || left.contains(&b'[') || left.contains(&b']') {
		left == right
	} else {
		left.eq_ignore_ascii_case(right)
	}
}

fn same_alleles(left: &Record, right: &Record) -> bool {
	left.allele_count() == right.allele_count()
		&& left
			.alleles()
			.zip(right.alleles())
			.all(|(left, right)| same_allele(left, right))
}

impl Groups {
	pub(super) fn build(&mut self, taken: &[Taken], ready: &[usize]) -> Result<(), Error> {
		self.used = 0;
		for &input in ready {
			let source = taken.get(input).ok_or(Error::Limits("inputs"))?;
			for slot in 0..source.count() {
				let record = source.record(slot).ok_or(Error::Limits("records"))?;
				let found = self.members.get(..self.used).and_then(|groups| {
					groups.iter().position(|group| {
						group.iter().all(|member| member.input != input)
							&& group
								.first()
								.and_then(|first| taken.get(first.input)?.record(first.slot))
								.is_some_and(|first| same_alleles(first, record))
					})
				});
				let index = match found {
					Some(index) => index,
					None => self.open()?,
				};
				self.members
					.get_mut(index)
					.ok_or(Error::Limits("groups"))?
					.push(Member { input, slot });
			}
		}
		self.order.clear();
		self.order.extend(0..self.used);
		let members = &self.members;
		self.order
			.sort_by_key(|&index| Reverse(members.get(index).map_or(0, Vec::len)));
		Ok(())
	}

	fn open(&mut self) -> Result<usize, Error> {
		if self.used == self.members.len() {
			self.members.push(Vec::new());
		}
		let index = self.used;
		self.used += 1;
		self.members
			.get_mut(index)
			.ok_or(Error::Limits("groups"))?
			.clear();
		Ok(index)
	}

	pub(super) fn ordered(&self) -> impl Iterator<Item = &[Member]> {
		self.order
			.iter()
			.filter_map(|&index| self.members.get(index).map(Vec::as_slice))
	}
}
