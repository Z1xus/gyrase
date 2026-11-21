use crate::Error;
use crate::header::Header;
use crate::record::{Field, Record};
use crate::value::{Cells, END_INT, MISSING_INT, Number, Numeric};

const fn genotype_index(low: usize, high: usize) -> usize {
	high * (high + 1) / 2 + low
}

fn nth_field(text: &[u8], index: usize) -> Option<&[u8]> {
	let mut parts = text.split(|byte| *byte == b',');
	let part = parts.nth(index)?;
	if part.is_empty() && parts.next().is_none() {
		return None;
	}
	Some(part)
}

fn push_fields(out: &mut Vec<u8>, text: &[u8], indices: &[usize]) {
	let start = out.len();
	for (position, &index) in indices.iter().enumerate() {
		let Some(part) = nth_field(text, index) else {
			out.truncate(start);
			out.push(b'.');
			return;
		};
		if position > 0 {
			out.push(b',');
		}
		out.extend_from_slice(part);
	}
}

fn pick_info<T: Numeric>(
	number: Number,
	src: &[T],
	alleles: usize,
	alt: usize,
	out: &mut Vec<T>,
) -> Option<()> {
	let count = src.iter().take_while(|value| !value.is_end()).count();
	let values = src.get(..count)?;
	let get = |index: usize| values.get(index).copied();
	match number {
		Number::PerAlt => {
			(count == alleles - 1).then_some(())?;
			out.push(get(alt)?);
		}
		Number::PerAllele => {
			(count == alleles).then_some(())?;
			out.extend([get(0)?, get(alt + 1)?]);
		}
		Number::PerGenotype => {
			(count == genotype_index(0, alleles)).then_some(())?;
			out.extend([
				get(0)?,
				get(genotype_index(0, alt + 1))?,
				get(genotype_index(alt + 1, alt + 1))?,
			]);
		}
		Number::Fixed(_) | Number::Variable => out.extend_from_slice(values),
	}
	Some(())
}

fn pick_format<T: Numeric>(
	number: Number,
	src: &[T],
	samples: usize,
	alleles: usize,
	alt: usize,
	out: &mut Vec<T>,
) -> Option<usize> {
	if src.len() == samples {
		out.extend_from_slice(src);
		return Some(1);
	}
	let width = src.len().checked_div(samples)?;
	let total = src.len();
	match number {
		Number::PerAlt => {
			(total == (alleles - 1) * samples).then_some(())?;
			for sample in src.chunks_exact(width) {
				let value = *sample.get(alt)?;
				out.push(if value.is_missing() || value.is_end() {
					T::MISSING
				} else {
					value
				});
			}
			Some(1)
		}
		Number::PerAllele => {
			(total == alleles * samples).then_some(())?;
			for sample in src.chunks_exact(width) {
				out.push(*sample.first()?);
				let value = *sample.get(alt + 1)?;
				out.push(if value.is_end() { T::MISSING } else { value });
			}
			Some(2)
		}
		Number::PerGenotype => {
			(total == genotype_index(0, alleles) * samples || total == alleles * samples)
				.then_some(())?;
			let all_haploid = width == alleles;
			for sample in src.chunks_exact(width) {
				out.push(*sample.first()?);
				if all_haploid || sample.iter().any(|value| value.is_end()) {
					out.push(*sample.get(alt + 1)?);
					if !all_haploid {
						out.push(T::END);
					}
				} else {
					out.push(*sample.get(genotype_index(0, alt + 1))?);
					out.push(*sample.get(genotype_index(alt + 1, alt + 1))?);
				}
			}
			Some(if all_haploid { 2 } else { 3 })
		}
		Number::Fixed(_) | Number::Variable => {
			out.extend_from_slice(src);
			Some(width)
		}
	}
}

fn text_indices(number: Number, alleles: usize, alt: usize, nfields: usize) -> Option<Vec<usize>> {
	match number {
		Number::PerAlt => Some(vec![alt]),
		Number::PerAllele => Some(vec![0, alt + 1]),
		Number::PerGenotype if nfields == alleles => Some(vec![0, alt + 1]),
		Number::PerGenotype => (nfields == genotype_index(0, alleles)).then(|| {
			vec![
				0,
				genotype_index(0, alt + 1),
				genotype_index(alt + 1, alt + 1),
			]
		}),
		Number::Fixed(_) | Number::Variable => None,
	}
}

fn remap_genotype(value: i32, alt: usize) -> i32 {
	let allele = value >> 1;
	let phase = value & 1;
	let alternate = i32::try_from(alt + 2).unwrap_or(i32::MAX);
	if value == MISSING_INT || allele == 0 || allele == 1 {
		value
	} else if allele == alternate {
		4 | phase
	} else {
		2 | phase
	}
}

#[derive(Debug, Default)]
pub(crate) struct Splitter {
	ints: Vec<i32>,
	floats: Vec<f32>,
	text: Vec<u8>,
	squeezed: Vec<u8>,
}

impl Splitter {
	pub(crate) fn split(
		&mut self,
		header: &Header,
		genotype: Option<u32>,
		src: &Record,
		alt: usize,
		dst: &mut Record,
	) -> Result<(), Error> {
		let (Some(reference), Some(alternate)) = (src.allele(0), src.allele(alt + 1)) else {
			return Err(Error::Vcf("record has no such allele"));
		};
		dst.clear();
		dst.set_site(src.contig, src.start);
		dst.set_id(&src.id);
		dst.set_qual(src.qual);
		dst.set_alleles([reference, alternate]);
		for &filter in &src.filters {
			dst.push_filter(filter);
		}
		for field in &src.info {
			self.info(header, src, field, alt, dst)?;
		}
		for field in &src.format {
			self.format(header, genotype, src, field, alt, dst)?;
		}
		Ok(())
	}

	fn info(
		&mut self,
		header: &Header,
		src: &Record,
		field: &Field,
		alt: usize,
		dst: &mut Record,
	) -> Result<(), Error> {
		let number = header
			.info_def(field.key)
			.map_or(Number::Variable, |definition| definition.number);
		let alleles = src.allele_count();
		let counted = match &field.cells {
			Cells::Flag => {
				dst.push_info_flag(field.key);
				Some(())
			}
			Cells::Ints(_) => {
				self.ints.clear();
				let picked = pick_info(
					number,
					src.store.ints(&field.cells),
					alleles,
					alt,
					&mut self.ints,
				);
				dst.push_info_ints(field.key, &self.ints);
				picked
			}
			Cells::Floats(_) => {
				self.floats.clear();
				let picked = pick_info(
					number,
					src.store.floats(&field.cells),
					alleles,
					alt,
					&mut self.floats,
				);
				dst.push_info_floats(field.key, &self.floats);
				picked
			}
			Cells::Text(_) => {
				self.text.clear();
				let text = src.store.text(&field.cells);
				match text_indices(number, alleles, alt, genotype_index(0, alleles)) {
					Some(indices) => push_fields(&mut self.text, text, &indices),
					None => self.text.extend_from_slice(text),
				}
				dst.push_info_text(field.key, &self.text);
				Some(())
			}
		};
		counted.ok_or_else(|| count_error(header, src, field.key))
	}

	fn format(
		&mut self,
		header: &Header,
		genotype: Option<u32>,
		src: &Record,
		field: &Field,
		alt: usize,
		dst: &mut Record,
	) -> Result<(), Error> {
		let number = header
			.format_def(field.key)
			.map_or(Number::Variable, |definition| definition.number);
		let samples = header.sample_count();
		let alleles = src.allele_count();
		let failed = || count_error(header, src, field.key);
		match &field.cells {
			Cells::Flag => return Err(Error::Bcf("FORMAT field without values")),
			Cells::Ints(_) if Some(field.key) == genotype => {
				self.ints.clear();
				let live = src.store.ints(&field.cells);
				for sample in live.chunks(field.width.max(1)) {
					let mut ended = false;
					for &value in sample {
						ended |= value == END_INT;
						self.ints.push(if ended {
							value
						} else {
							remap_genotype(value, alt)
						});
					}
				}
				dst.push_format_ints(field.key, field.width, &self.ints);
			}
			Cells::Ints(_) => {
				self.ints.clear();
				let width = pick_format(
					number,
					src.store.ints(&field.cells),
					samples,
					alleles,
					alt,
					&mut self.ints,
				)
				.ok_or_else(failed)?;
				dst.push_format_ints(field.key, width, &self.ints);
			}
			Cells::Floats(_) => {
				self.floats.clear();
				let width = pick_format(
					number,
					src.store.floats(&field.cells),
					samples,
					alleles,
					alt,
					&mut self.floats,
				)
				.ok_or_else(failed)?;
				dst.push_format_floats(field.key, width, &self.floats);
			}
			Cells::Text(_) => {
				let width = self
					.format_text(
						number,
						src.store.text(&field.cells),
						field.width,
						alleles,
						alt,
					)
					.ok_or_else(failed)?;
				dst.push_format_text(field.key, width, &self.text);
			}
		}
		Ok(())
	}

	fn format_text(
		&mut self,
		number: Number,
		src: &[u8],
		width: usize,
		alleles: usize,
		alt: usize,
	) -> Option<usize> {
		self.text.clear();
		if width == 0 || matches!(number, Number::Fixed(_) | Number::Variable) {
			self.text.extend_from_slice(src);
			return Some(width);
		}
		let mut longest = 0;
		for sample in src.chunks_exact(width) {
			let text = sample.split(|byte| *byte == 0).next().unwrap_or_default();
			let nfields = text.split(|byte| *byte == b',').count();
			let missing = text == b".";
			let start = self.text.len();
			match (missing, text_indices(number, alleles, alt, nfields)) {
				(true, _) => self.text.push(b'.'),
				(false, Some(indices)) => push_fields(&mut self.text, text, &indices),
				(false, None) => return None,
			}
			longest = longest.max(self.text.len() - start);
			self.text.resize(start + width, 0);
		}
		if longest < width {
			self.squeezed.clear();
			for sample in self.text.chunks_exact(width) {
				self.squeezed.extend_from_slice(sample.get(..longest)?);
			}
			std::mem::swap(&mut self.text, &mut self.squeezed);
		}
		Some(longest.min(width))
	}
}

fn count_error(header: &Header, record: &Record, key: u32) -> Error {
	Error::ValueCount {
		tag: header.name(key).unwrap_or_default().to_owned(),
		contig: header
			.contig(record.contig)
			.map(|contig| contig.name.clone())
			.unwrap_or_default(),
		position: u64::from(record.start) + 1,
	}
}
