use crate::Error;
use crate::header::Header;
use crate::value::{Cells, END_INT, MISSING_INT, Store};

#[derive(Debug, Clone, Copy)]
pub struct CountKeys {
	pub genotype: u32,
	pub allele_count: u32,
	pub allele_number: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Field {
	pub(crate) key: u32,
	pub(crate) width: usize,
	pub(crate) cells: Cells,
}

impl Field {
	pub fn key(&self) -> u32 {
		self.key
	}

	pub fn width(&self) -> usize {
		self.width
	}

	pub fn cells(&self) -> &Cells {
		&self.cells
	}
}

#[derive(Debug, Clone, Default)]
pub struct Record {
	pub(crate) contig: u32,
	pub(crate) start: u32,
	pub(crate) id: Vec<u8>,
	pub(crate) allele_bytes: Vec<u8>,
	pub(crate) allele_ends: Vec<usize>,
	pub(crate) qual: Option<f32>,
	pub(crate) filters: Vec<u32>,
	pub(crate) info: Vec<Field>,
	pub(crate) format: Vec<Field>,
	pub(crate) store: Store,
}

impl Record {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn clear(&mut self) {
		self.contig = 0;
		self.start = 0;
		self.id.clear();
		self.allele_bytes.clear();
		self.allele_ends.clear();
		self.qual = None;
		self.filters.clear();
		self.info.clear();
		self.format.clear();
		self.store.clear();
	}

	pub fn set_site(&mut self, contig: u32, start: u32) {
		self.contig = contig;
		self.start = start;
	}

	pub fn set_id(&mut self, id: &[u8]) {
		self.id.clear();
		self.id.extend_from_slice(id);
	}

	pub fn push_allele(&mut self, allele: &[u8]) {
		self.allele_bytes.extend_from_slice(allele);
		self.allele_ends.push(self.allele_bytes.len());
	}

	pub fn set_alleles<'a>(&mut self, alleles: impl IntoIterator<Item = &'a [u8]>) {
		self.allele_bytes.clear();
		self.allele_ends.clear();
		for allele in alleles {
			self.push_allele(allele);
		}
	}

	pub(crate) fn truncate_reference(&mut self, len: usize) {
		let Some(&first_end) = self.allele_ends.first() else {
			return;
		};
		let new_end = len.min(first_end);
		let removed = first_end - new_end;
		self.allele_bytes.drain(new_end..first_end);
		for end in &mut self.allele_ends {
			*end -= removed;
		}
	}

	pub(crate) fn remap(&mut self, contigs: &[u32], keys: &[Option<u32>]) -> Result<(), Error> {
		let lookup = |key: u32| {
			keys.get(key as usize)
				.copied()
				.flatten()
				.ok_or(Error::Bcf("tag is not in the output header"))
		};
		self.contig = *contigs
			.get(self.contig as usize)
			.ok_or(Error::Bcf("contig is not in the output header"))?;
		for filter in &mut self.filters {
			*filter = lookup(*filter)?;
		}
		for field in self.info.iter_mut().chain(&mut self.format) {
			field.key = lookup(field.key)?;
		}
		Ok(())
	}

	pub fn set_start(&mut self, start: u32) {
		self.start = start;
	}

	pub fn set_qual(&mut self, qual: Option<f32>) {
		self.qual = qual;
	}

	pub fn push_filter(&mut self, key: u32) {
		self.filters.push(key);
	}

	pub fn push_info_flag(&mut self, key: u32) {
		self.info.push(Field {
			key,
			width: 0,
			cells: Cells::Flag,
		});
	}

	pub fn push_info_ints(&mut self, key: u32, values: &[i32]) {
		let cells = self.store.push_ints(values.iter().copied());
		self.info.push(Field {
			key,
			width: values.len(),
			cells,
		});
	}

	pub fn push_info_floats(&mut self, key: u32, values: &[f32]) {
		let cells = self.store.push_floats(values.iter().copied());
		self.info.push(Field {
			key,
			width: values.len(),
			cells,
		});
	}

	pub fn push_info_text(&mut self, key: u32, value: &[u8]) {
		let cells = self.store.push_text(value);
		self.info.push(Field {
			key,
			width: value.len(),
			cells,
		});
	}

	pub fn push_format_ints(&mut self, key: u32, width: usize, values: &[i32]) {
		let cells = self.store.push_ints(values.iter().copied());
		self.format.push(Field { key, width, cells });
	}

	pub fn push_format_floats(&mut self, key: u32, width: usize, values: &[f32]) {
		let cells = self.store.push_floats(values.iter().copied());
		self.format.push(Field { key, width, cells });
	}

	pub fn push_format_text(&mut self, key: u32, width: usize, values: &[u8]) {
		let cells = self.store.push_text(values);
		self.format.push(Field { key, width, cells });
	}

	pub fn contig(&self) -> u32 {
		self.contig
	}

	pub fn start(&self) -> u32 {
		self.start
	}

	pub fn id(&self) -> &[u8] {
		&self.id
	}

	pub fn allele_count(&self) -> usize {
		self.allele_ends.len()
	}

	pub fn allele(&self, index: usize) -> Option<&[u8]> {
		let end = *self.allele_ends.get(index)?;
		let start = match index.checked_sub(1) {
			Some(previous) => *self.allele_ends.get(previous)?,
			None => 0,
		};
		self.allele_bytes.get(start..end)
	}

	pub fn alleles(&self) -> impl Iterator<Item = &[u8]> {
		(0..self.allele_count()).filter_map(|index| self.allele(index))
	}

	pub fn qual(&self) -> Option<f32> {
		self.qual
	}

	pub fn filters(&self) -> &[u32] {
		&self.filters
	}

	pub fn info(&self) -> &[Field] {
		&self.info
	}

	pub fn format(&self) -> &[Field] {
		&self.format
	}

	pub fn store(&self) -> &Store {
		&self.store
	}

	pub fn info_field(&self, key: u32) -> Option<&Field> {
		self.info.iter().find(|field| field.key == key)
	}

	pub fn format_field(&self, key: u32) -> Option<&Field> {
		self.format.iter().find(|field| field.key == key)
	}

	pub fn reference_len(&self) -> u32 {
		let len = self.allele(0).map_or(0, <[u8]>::len);
		u32::try_from(len).unwrap_or(u32::MAX)
	}

	pub fn end(&self, header: &Header) -> u32 {
		let reference_end = self.start.saturating_add(self.reference_len());
		let declared = header
			.key("END")
			.and_then(|key| self.info_field(key))
			.and_then(|field| self.store.ints(&field.cells).first().copied())
			.and_then(|end| u32::try_from(end).ok());
		declared.map_or(reference_end, |end| end.max(reference_end))
	}

	pub fn keep_samples(&mut self, keep: &[usize]) -> Result<(), Error> {
		for field in &mut self.format {
			let width = field.width;
			let (range, copied) = match &field.cells {
				Cells::Ints(range) => (range.clone(), self.store.ints.len()),
				Cells::Floats(range) => (range.clone(), self.store.floats.len()),
				Cells::Text(range) => (range.clone(), self.store.text.len()),
				Cells::Flag => return Err(Error::Bcf("FORMAT field without values")),
			};
			for &sample in keep {
				let from = sample
					.checked_mul(width)
					.and_then(|offset| range.start.checked_add(offset))
					.filter(|&from| from + width <= range.end)
					.ok_or(Error::Bcf("sample is outside the FORMAT field"))?;
				match &field.cells {
					Cells::Ints(_) => self.store.ints.extend_from_within(from..from + width),
					Cells::Floats(_) => self.store.floats.extend_from_within(from..from + width),
					_ => self.store.text.extend_from_within(from..from + width),
				}
			}
			field.cells = match &field.cells {
				Cells::Ints(_) => Cells::Ints(copied..self.store.ints.len()),
				Cells::Floats(_) => Cells::Floats(copied..self.store.floats.len()),
				_ => Cells::Text(copied..self.store.text.len()),
			};
		}
		Ok(())
	}

	pub fn set_info_ints(&mut self, key: u32, values: &[i32]) {
		let position = self.info.iter().position(|field| field.key == key);
		if values.is_empty() {
			if let Some(position) = position {
				self.info.remove(position);
			}
			return;
		}
		let field = Field {
			key,
			width: values.len(),
			cells: self.store.push_ints(values.iter().copied()),
		};
		match position.and_then(|position| self.info.get_mut(position)) {
			Some(existing) => *existing = field,
			None => self.info.push(field),
		}
	}

	pub(crate) fn count_alleles(
		&self,
		genotype: u32,
		counts: &mut Vec<i32>,
	) -> Result<bool, Error> {
		let Some(genotypes) = self.format_field(genotype) else {
			return Ok(false);
		};
		counts.clear();
		counts.resize(self.allele_count(), 0);
		for &value in self.store.ints(&genotypes.cells) {
			if value == END_INT || value == MISSING_INT {
				continue;
			}
			let Ok(index) = usize::try_from((value >> 1) - 1) else {
				continue;
			};
			*counts
				.get_mut(index)
				.ok_or(Error::Bcf("genotype allele is not in the record"))? += 1;
		}
		Ok(true)
	}

	pub fn update_counts(&mut self, keys: CountKeys, counts: &mut Vec<i32>) -> Result<(), Error> {
		if !self.count_alleles(keys.genotype, counts)? {
			return Ok(());
		}
		let total = counts.iter().sum();
		self.set_info_ints(keys.allele_count, counts.get(1..).unwrap_or_default());
		self.set_info_ints(keys.allele_number, &[total]);
		Ok(())
	}

	pub fn check(&self, header: &Header) -> Result<(), Error> {
		let samples = header.sample_count();
		for field in &self.format {
			let len = match &field.cells {
				Cells::Flag => return Err(Error::Bcf("FORMAT field without values")),
				Cells::Ints(range) | Cells::Floats(range) | Cells::Text(range) => range.len(),
			};
			if field.width == 0 || len != field.width.saturating_mul(samples) {
				return Err(Error::Bcf("FORMAT field does not cover every sample"));
			}
		}
		Ok(())
	}
}
