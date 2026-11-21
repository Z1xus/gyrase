use std::iter::repeat_n;

use super::Merge;
use super::group::Group;
use super::input::Taken;
use crate::Error;
use crate::record::{Field, Record};
use crate::value::{Cells, END_INT, MISSING_INT, Number, Numeric, ValueType};

#[derive(Debug, Default)]
struct Sum {
	ints: Vec<i32>,
	floats: Vec<f32>,
	blocks: usize,
	width: usize,
}

#[derive(Debug)]
enum Shape {
	Ints,
	Floats,
	Text,
}

#[derive(Debug)]
struct Spread {
	key: u32,
	shape: Shape,
	ints: Vec<i32>,
	floats: Vec<f32>,
	text: Vec<u8>,
}

#[derive(Debug)]
pub(super) struct Combiner {
	samples: Vec<usize>,
	present: Vec<Option<usize>>,
	sums: Vec<Sum>,
	spreads: Vec<Spread>,
	spread_count: usize,
	seen: Vec<u32>,
	keys: Vec<u32>,
	filters: Vec<u32>,
	ids: Vec<u8>,
	ints: Vec<i32>,
	floats: Vec<f32>,
	text: Vec<u8>,
	counts: Vec<i32>,
	slots: Vec<usize>,
	ploidy: Vec<usize>,
}

fn nonzero(count: usize) -> usize {
	count.max(1)
}

fn copy_info(out: &mut Record, source: &Record, field: &Field) {
	match &field.cells {
		Cells::Flag => out.push_info_flag(field.key),
		Cells::Ints(_) => out.push_info_ints(field.key, source.store.ints(&field.cells)),
		Cells::Floats(_) => out.push_info_floats(field.key, source.store.floats(&field.cells)),
		Cells::Text(_) => out.push_info_text(field.key, source.store.text(&field.cells)),
	}
}

fn live<T: Numeric>(values: &[T]) -> &[T] {
	let count = values.iter().take_while(|value| !value.is_end()).count();
	values.get(..count).unwrap_or(values)
}

fn spread_size(number: Number, alleles: usize) -> usize {
	match number {
		Number::PerAlt => alleles.saturating_sub(1),
		Number::PerAllele => alleles,
		Number::PerGenotype => alleles * (alleles + 1) / 2,
		Number::Fixed(_) | Number::Variable => 0,
	}
}

fn scatter<T: Numeric>(target: &mut [T], values: &[T], slots: &[usize]) {
	for (&slot, &value) in slots.iter().zip(live(values)) {
		if !value.is_missing()
			&& let Some(target) = target.get_mut(slot)
		{
			*target = value;
		}
	}
}

fn size_alleles(number: Number, size: usize) -> usize {
	match number {
		Number::PerAlt => size + 1,
		Number::PerGenotype => (1..=size)
			.find(|&n| n * (n + 1) / 2 >= size)
			.unwrap_or(size),
		Number::PerAllele | Number::Fixed(_) | Number::Variable => size,
	}
}

fn remapped(map: &[usize], alleles: usize) -> bool {
	map.len() != alleles || map.iter().enumerate().any(|(at, &target)| at != target)
}

fn slots(number: Number, map: &[usize], out: &mut Vec<usize>) {
	out.clear();
	match number {
		Number::PerAlt => out.extend(map.iter().skip(1).map(|&target| target.saturating_sub(1))),
		Number::PerAllele => out.extend_from_slice(map),
		Number::PerGenotype => {
			for (row, &first) in map.iter().enumerate() {
				for &second in map.iter().take(row + 1) {
					let (high, low) = if first > second {
						(first, second)
					} else {
						(second, first)
					};
					out.push(high * (high + 1) / 2 + low);
				}
			}
		}
		Number::Fixed(_) | Number::Variable => {}
	}
}

fn genotype_allele(value: i32, map: &[usize]) -> i32 {
	let allele = (value >> 1) - 1;
	let Ok(index) = usize::try_from(allele) else {
		return value;
	};
	match map.get(index) {
		Some(&target) if index > 0 => {
			i32::try_from(target + 1).map_or(value, |target| (target << 1) | (value & 1))
		}
		_ => value,
	}
}

#[derive(Debug)]
struct Plan<'a> {
	number: Number,
	map: &'a [usize],
	slots: &'a [usize],
	alleles: usize,
}

fn gather<T: Numeric>(
	out: &mut Vec<T>,
	samples: usize,
	found: Option<(&[T], usize)>,
	plan: Option<&Plan<'_>>,
	ploidy: &[usize],
	size: usize,
) -> Result<(), Error> {
	let Some((values, width)) = found else {
		for _ in 0..samples {
			out.push(T::MISSING);
			out.extend(repeat_n(T::END, size.saturating_sub(1)));
		}
		return Ok(());
	};
	let chunks = values.chunks_exact(width.max(1));
	if chunks.len() < samples {
		return Err(Error::Bcf("FORMAT field does not cover every sample"));
	}
	for (index, chunk) in chunks.take(samples).enumerate() {
		let Some(plan) = plan else {
			let kept = live(chunk);
			let kept = kept.get(..kept.len().min(size)).unwrap_or(kept);
			out.extend_from_slice(kept);
			out.extend(repeat_n(T::END, size - kept.len()));
			continue;
		};
		let lone = chunk.first().is_some_and(|first| first.is_missing())
			&& chunk.get(1).is_none_or(|next| next.is_end());
		if lone {
			out.push(T::MISSING);
			out.extend(repeat_n(T::END, size.saturating_sub(1)));
			continue;
		}
		let haploid = plan.number == Number::PerGenotype && ploidy.get(index) == Some(&1);
		let used = if haploid { plan.alleles } else { size };
		let start = out.len();
		out.extend(repeat_n(T::MISSING, used));
		out.extend(repeat_n(T::END, size.saturating_sub(used)));
		let positions = if haploid { plan.map } else { plan.slots };
		for (&slot, &value) in positions.iter().zip(chunk) {
			if value.is_end() {
				break;
			}
			if !value.is_missing()
				&& let Some(target) = out.get_mut(start + slot)
			{
				*target = value;
			}
		}
	}
	Ok(())
}

fn sum_floats(sum: &mut Sum) -> Vec<f32> {
	for value in &mut sum.floats {
		if value.is_missing() {
			*value = 0.0;
		}
	}
	let mut blocks = sum.floats.chunks_exact(sum.width);
	let mut total = blocks.next().unwrap_or_default().to_vec();
	for block in blocks {
		for (total, value) in total.iter_mut().zip(block) {
			*total += value;
		}
	}
	total
}

fn sum_ints(sum: &mut Sum) -> Vec<i32> {
	for value in &mut sum.ints {
		if value.is_missing() {
			*value = 0;
		}
	}
	let mut blocks = sum.ints.chunks_exact(sum.width);
	let mut total = blocks.next().unwrap_or_default().to_vec();
	for block in blocks {
		for (total, value) in total.iter_mut().zip(block) {
			*total = total.wrapping_add(*value);
		}
	}
	total
}

fn records<'a>(taken: &'a [Taken], group: &Group) -> Result<Vec<&'a Record>, Error> {
	group
		.members
		.iter()
		.map(|member| {
			taken
				.get(member.input)
				.and_then(|source| source.record(member.slot))
				.ok_or(Error::Limits("records"))
		})
		.collect()
}

#[derive(Debug)]
struct Layout<'a> {
	samples: &'a [usize],
	ploidy: &'a [usize],
	present: &'a [Option<usize>],
	taken: &'a [Taken],
	group: &'a Group,
	key: u32,
	number: Number,
	alleles: usize,
	size: usize,
}

fn gather_all<'a, T: Numeric + 'a>(
	out: &mut Vec<T>,
	buffer: &mut Vec<usize>,
	layout: &Layout<'a>,
	values: impl Fn(&'a Record, &'a Field) -> &'a [T],
) -> Result<(), Error> {
	let mut offset = 0;
	for input in 0..layout.samples.len() {
		let samples = layout.samples.get(input).copied().unwrap_or(0);
		let found = source(layout.present, layout.taken, input, layout.key)
			.map(|(record, field)| (values(record, field), field.width));
		let map = layout.group.map(input).unwrap_or_default();
		let remap = matches!(
			layout.number,
			Number::PerAlt | Number::PerAllele | Number::PerGenotype
		) && remapped(map, layout.alleles);
		if remap {
			slots(layout.number, map, buffer);
		}
		let plan = remap.then_some(Plan {
			number: layout.number,
			map,
			slots: buffer.as_slice(),
			alleles: layout.alleles,
		});
		let ploidy = layout.ploidy.get(offset..).unwrap_or_default();
		gather(out, samples, found, plan.as_ref(), ploidy, layout.size)?;
		offset += samples;
	}
	Ok(())
}

fn source<'a>(
	present: &[Option<usize>],
	taken: &'a [Taken],
	input: usize,
	key: u32,
) -> Option<(&'a Record, &'a Field)> {
	let slot = present.get(input).copied().flatten()?;
	let record = taken.get(input)?.record(slot)?;
	Some((record, record.format_field(key)?))
}

impl Combiner {
	pub(super) fn new(merge: &Merge, samples: &[usize]) -> Self {
		Self {
			samples: samples.to_vec(),
			present: vec![None; samples.len()],
			sums: merge.rules.iter().map(|_| Sum::default()).collect(),
			spreads: Vec::new(),
			spread_count: 0,
			seen: Vec::new(),
			keys: Vec::new(),
			filters: Vec::new(),
			ids: Vec::new(),
			ints: Vec::new(),
			floats: Vec::new(),
			text: Vec::new(),
			counts: Vec::new(),
			slots: Vec::new(),
			ploidy: Vec::new(),
		}
	}

	pub(super) fn combine(
		&mut self,
		merge: &Merge,
		taken: &[Taken],
		group: &Group,
		out: &mut Record,
	) -> Result<(), Error> {
		self.present.fill(None);
		for member in &group.members {
			let slot = self
				.present
				.get_mut(member.input)
				.ok_or(Error::Limits("inputs"))?;
			*slot = Some(member.slot);
		}
		let records = records(taken, group)?;
		let first = *records.first().ok_or(Error::Limits("empty group"))?;
		out.clear();
		out.set_site(first.contig(), first.start());
		self.identity(out, &records);
		out.set_alleles(group.alleles.iter().map(Vec::as_slice));
		self.filters(merge, out, &records);
		let alleles = group.alleles.len();
		self.info(merge, group, out, &records, alleles)?;
		self.format(merge, taken, group, out, &records, alleles)
	}

	fn identity(&mut self, out: &mut Record, records: &[&Record]) {
		self.ids.clear();
		let mut quality: Option<f32> = None;
		for record in records {
			let id = record.id();
			let known = self.ids.split(|byte| *byte == b';').any(|part| part == id);
			if !id.is_empty() && id != b"." && !known {
				if !self.ids.is_empty() {
					self.ids.push(b';');
				}
				self.ids.extend_from_slice(id);
			}
			if let Some(value) = record.qual() {
				quality =
					Some(quality.map_or(value, |best| if best < value { value } else { best }));
			}
		}
		out.set_id(&self.ids);
		out.set_qual(quality);
	}

	fn filters(&mut self, merge: &Merge, out: &mut Record, records: &[&Record]) {
		self.filters.clear();
		for record in records {
			for &filter in record.filters() {
				if !self.filters.contains(&filter) {
					self.filters.push(filter);
				}
			}
		}
		if self.filters.len() > 1 {
			self.filters.retain(|&filter| Some(filter) != merge.pass);
		}
		for &filter in &self.filters {
			out.push_filter(filter);
		}
	}

	fn info(
		&mut self,
		merge: &Merge,
		group: &Group,
		out: &mut Record,
		records: &[&Record],
		alleles: usize,
	) -> Result<(), Error> {
		for sum in &mut self.sums {
			sum.ints.clear();
			sum.floats.clear();
			sum.blocks = 0;
			sum.width = 0;
		}
		self.spread_count = 0;
		self.seen.clear();
		for (member, record) in group.members.iter().zip(records) {
			for field in record.info() {
				if merge.skipped.contains(&field.key) {
					continue;
				}
				if let Some(index) = merge.rules.iter().position(|rule| rule.key == field.key) {
					self.add_sum(index, record, field)?;
					continue;
				}
				let number = merge
					.header
					.info_def(field.key)
					.map_or(Number::Variable, |definition| definition.number);
				let size = spread_size(number, alleles);
				if size > 0 {
					self.add_spread(record, field, size, number, &member.map)?;
				} else if !self.seen.contains(&field.key) {
					self.seen.push(field.key);
					copy_info(out, record, field);
				}
			}
		}
		self.finish_sums(merge, out);
		self.finish_spreads(out);
		Ok(())
	}

	fn add_sum(&mut self, index: usize, record: &Record, field: &Field) -> Result<(), Error> {
		let sum = self.sums.get_mut(index).ok_or(Error::Limits("rules"))?;
		let width = match &field.cells {
			Cells::Ints(_) => {
				let values = live(record.store.ints(&field.cells));
				sum.ints.extend_from_slice(values);
				values.len()
			}
			Cells::Floats(_) => {
				let values = live(record.store.floats(&field.cells));
				sum.floats.extend_from_slice(values);
				values.len()
			}
			Cells::Flag | Cells::Text(_) => return Err(Error::Incompatible("INFO type")),
		};
		if sum.blocks > 0 && sum.width != width {
			return Err(Error::Incompatible("INFO value counts differ"));
		}
		sum.width = width;
		sum.blocks += 1;
		Ok(())
	}

	fn finish_sums(&mut self, merge: &Merge, out: &mut Record) {
		for (rule, sum) in merge.rules.iter().zip(&mut self.sums) {
			if sum.blocks == 0 || sum.width == 0 {
				continue;
			}
			if rule.value_type == ValueType::Float {
				out.push_info_floats(rule.key, &sum_floats(sum));
			} else {
				out.push_info_ints(rule.key, &sum_ints(sum));
			}
		}
	}

	fn spread_index(&mut self, field: &Field, size: usize) -> usize {
		let known = self
			.spreads
			.get(..self.spread_count)
			.and_then(|spreads| spreads.iter().position(|spread| spread.key == field.key));
		if let Some(index) = known {
			return index;
		}
		if self.spread_count == self.spreads.len() {
			self.spreads.push(Spread {
				key: field.key,
				shape: Shape::Text,
				ints: Vec::new(),
				floats: Vec::new(),
				text: Vec::new(),
			});
		}
		let index = self.spread_count;
		self.spread_count += 1;
		if let Some(spread) = self.spreads.get_mut(index) {
			spread.key = field.key;
			spread.ints.clear();
			spread.floats.clear();
			spread.text.clear();
			spread.shape = match &field.cells {
				Cells::Ints(_) | Cells::Flag => {
					spread.ints.resize(size, MISSING_INT);
					Shape::Ints
				}
				Cells::Floats(_) => {
					spread.floats.resize(size, f32::MISSING);
					Shape::Floats
				}
				Cells::Text(_) => Shape::Text,
			};
		}
		index
	}

	fn add_spread(
		&mut self,
		record: &Record,
		field: &Field,
		size: usize,
		number: Number,
		map: &[usize],
	) -> Result<(), Error> {
		let index = self.spread_index(field, size);
		slots(number, map, &mut self.slots);
		let Some(spread) = self.spreads.get_mut(index) else {
			return Ok(());
		};
		match &field.cells {
			Cells::Ints(_) => scatter(
				&mut spread.ints,
				record.store.ints(&field.cells),
				&self.slots,
			),
			Cells::Floats(_) => {
				scatter(
					&mut spread.floats,
					record.store.floats(&field.cells),
					&self.slots,
				);
			}
			Cells::Text(_) => {
				if remapped(map, size_alleles(number, size)) {
					return Err(Error::Incompatible(
						"string INFO tag with different alleles",
					));
				}
				spread.text.clear();
				spread
					.text
					.extend_from_slice(record.store.text(&field.cells));
			}
			Cells::Flag => {}
		}
		Ok(())
	}

	fn finish_spreads(&mut self, out: &mut Record) {
		for spread in self.spreads.iter().take(self.spread_count) {
			match spread.shape {
				Shape::Ints => out.push_info_ints(spread.key, &spread.ints),
				Shape::Floats => out.push_info_floats(spread.key, &spread.floats),
				Shape::Text => out.push_info_text(spread.key, &spread.text),
			}
		}
	}

	fn format(
		&mut self,
		merge: &Merge,
		taken: &[Taken],
		group: &Group,
		out: &mut Record,
		records: &[&Record],
		alleles: usize,
	) -> Result<(), Error> {
		self.keys.clear();
		self.ploidy.clear();
		let mut has_genotype = false;
		for record in records {
			for field in record.format() {
				if Some(field.key) == merge.genotype {
					has_genotype = true;
				} else if !self.keys.contains(&field.key) {
					self.keys.push(field.key);
				}
			}
		}
		if let (true, Some(key)) = (has_genotype, merge.genotype) {
			self.genotype(taken, group, key, out)?;
			if let Some(keys) = merge.counts
				&& out.count_alleles(keys.genotype, &mut self.counts)?
			{
				let total = self.counts.iter().sum();
				out.set_info_ints(keys.allele_number, &[total]);
				out.set_info_ints(keys.allele_count, self.counts.get(1..).unwrap_or_default());
			}
		}
		for index in 0..self.keys.len() {
			let key = *self.keys.get(index).ok_or(Error::Limits("tags"))?;
			self.field(merge, taken, group, key, alleles, out)?;
		}
		Ok(())
	}

	fn genotype(
		&mut self,
		taken: &[Taken],
		group: &Group,
		key: u32,
		out: &mut Record,
	) -> Result<(), Error> {
		let mut size = 1;
		for input in 0..self.samples.len() {
			if let Some((record, field)) = source(&self.present, taken, input, key) {
				for sample in record.store.ints(&field.cells).chunks(nonzero(field.width)) {
					size = size.max(live(sample).len());
				}
			}
		}
		self.ints.clear();
		self.ploidy.clear();
		for input in 0..self.samples.len() {
			let samples = self.samples.get(input).copied().unwrap_or(0);
			match source(&self.present, taken, input, key) {
				None => {
					self.ints.extend(repeat_n(0, samples * size));
					self.ploidy.extend(repeat_n(size, samples));
				}
				Some((record, field)) => {
					let values = record.store.ints(&field.cells);
					let chunks = values.chunks_exact(nonzero(field.width));
					if chunks.len() < samples {
						return Err(Error::Bcf("FORMAT field does not cover every sample"));
					}
					let map = group.map(input).unwrap_or_default();
					for sample in chunks.take(samples) {
						let kept = live(sample);
						self.ploidy.push(kept.len());
						self.ints.extend(kept.iter().map(|&value| {
							if value == MISSING_INT {
								0
							} else {
								genotype_allele(value, map)
							}
						}));
						self.ints
							.extend(repeat_n(END_INT, size.saturating_sub(kept.len())));
					}
				}
			}
		}
		out.push_format_ints(key, size, &self.ints);
		Ok(())
	}

	fn field(
		&mut self,
		merge: &Merge,
		taken: &[Taken],
		group: &Group,
		key: u32,
		alleles: usize,
		out: &mut Record,
	) -> Result<(), Error> {
		let number = merge
			.header
			.format_def(key)
			.map_or(Number::Variable, |definition| definition.number);
		let mut widest = 0;
		let mut cells = None;
		for input in 0..self.samples.len() {
			if let Some((_, field)) = source(&self.present, taken, input, key) {
				widest = widest.max(field.width);
				cells.get_or_insert(field.cells.clone());
			}
		}
		let size = match spread_size(number, alleles) {
			0 => nonzero(widest),
			size => size,
		};
		let layout = Layout {
			samples: &self.samples,
			ploidy: &self.ploidy,
			present: &self.present,
			taken,
			group,
			key,
			number,
			alleles,
			size,
		};
		match cells {
			Some(Cells::Ints(_)) => {
				self.ints.clear();
				gather_all(&mut self.ints, &mut self.slots, &layout, |record, field| {
					record.store.ints(&field.cells)
				})?;
				out.push_format_ints(key, size, &self.ints);
			}
			Some(Cells::Floats(_)) => {
				self.floats.clear();
				gather_all(
					&mut self.floats,
					&mut self.slots,
					&layout,
					|record, field| record.store.floats(&field.cells),
				)?;
				out.push_format_floats(key, size, &self.floats);
			}
			Some(Cells::Text(_)) => {
				let changed = group
					.members
					.iter()
					.any(|member| remapped(&member.map, alleles));
				if changed && spread_size(number, alleles) > 0 {
					return Err(Error::Incompatible(
						"string FORMAT tag with different alleles",
					));
				}
				let width = self.text_field(taken, key, number, size)?;
				out.push_format_text(key, width, &self.text);
			}
			Some(Cells::Flag) | None => return Err(Error::Bcf("FORMAT field without values")),
		}
		Ok(())
	}

	fn text_field(
		&mut self,
		taken: &[Taken],
		key: u32,
		number: Number,
		size: usize,
	) -> Result<usize, Error> {
		let default_len = match number {
			Number::PerAlt | Number::PerAllele | Number::PerGenotype => size * 2 - 1,
			Number::Fixed(_) | Number::Variable => 1,
		};
		let mut width = default_len;
		for input in 0..self.samples.len() {
			if let Some((_, field)) = source(&self.present, taken, input, key) {
				width = width.max(field.width);
			}
		}
		self.text.clear();
		for input in 0..self.samples.len() {
			let samples = self.samples.get(input).copied().unwrap_or(0);
			match source(&self.present, taken, input, key) {
				None => {
					for _ in 0..samples {
						let start = self.text.len();
						self.text.push(b'.');
						for _ in 1..default_len.div_ceil(2) {
							self.text.extend_from_slice(b",.");
						}
						self.text.resize(start + width, 0);
					}
				}
				Some((record, field)) => {
					let values = record.store.text(&field.cells);
					let chunks = values.chunks_exact(nonzero(field.width));
					if chunks.len() < samples {
						return Err(Error::Bcf("FORMAT field does not cover every sample"));
					}
					for sample in chunks.take(samples) {
						let start = self.text.len();
						self.text.extend_from_slice(sample);
						self.text.resize(start + width, 0);
					}
				}
			}
		}
		Ok(width)
	}
}
