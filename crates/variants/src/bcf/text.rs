use super::record::{BcfRecord, FormatEntry};
use super::typed::Typed;
use crate::Error;
use crate::header::Header;
use crate::number::{push_float, push_int};
use crate::value::{BcfType, is_missing_float};
use crate::vcf::{
	push_float_list, push_int_list, push_or_dot, push_text, tag_name, write_genotype,
};

const GENOTYPE_KEY: &str = "GT";
const MISSING_STRING_BYTE: u8 = 7;

impl BcfRecord<'_> {
	pub fn write_vcf(&self, header: &Header, out: &mut Vec<u8>) -> Result<(), Error> {
		let contig = header
			.contig(self.contig()?)
			.ok_or(Error::Bcf("contig id is not in the header"))?;
		out.extend_from_slice(contig.name.as_bytes());
		out.push(b'\t');
		push_int(out, i64::from(self.start()?) + 1);
		out.push(b'\t');
		let id_start = out.len();
		push_text(out, self.id()?);
		if out.len() == id_start {
			out.push(b'.');
		}
		out.push(b'\t');
		let mut alleles = self.alleles()?;
		push_or_dot(out, alleles.next().transpose()?.unwrap_or_default());
		out.push(b'\t');
		if self.allele_count() < 2 {
			out.push(b'.');
		}
		for (index, allele) in alleles.enumerate() {
			if index > 0 {
				out.push(b',');
			}
			out.extend_from_slice(allele?);
		}
		out.push(b'\t');
		match self.qual()? {
			Some(qual) if !is_missing_float(qual) => push_float(out, qual),
			_ => out.push(b'.'),
		}
		out.push(b'\t');
		self.write_filters(header, out)?;
		out.push(b'\t');
		self.write_info(header, out)?;
		if header.sample_count() > 0 {
			self.write_samples(header, out)?;
		}
		Ok(())
	}

	fn write_filters(&self, header: &Header, out: &mut Vec<u8>) -> Result<(), Error> {
		let filters = self.filters()?;
		if matches!(filters.ty, BcfType::Null) || filters.len == 0 {
			out.push(b'.');
			return Ok(());
		}
		for (index, key) in filters.ints().enumerate() {
			if index > 0 {
				out.push(b';');
			}
			let key = u32::try_from(key).map_err(|_| Error::Bcf("negative FILTER key"))?;
			out.extend_from_slice(tag_name(header, key)?.as_bytes());
		}
		Ok(())
	}

	fn write_info(&self, header: &Header, out: &mut Vec<u8>) -> Result<(), Error> {
		if self.info_count() == 0 {
			out.push(b'.');
			return Ok(());
		}
		for (index, entry) in self.info().enumerate() {
			let entry = entry?;
			if index > 0 {
				out.push(b';');
			}
			out.extend_from_slice(tag_name(header, entry.key)?.as_bytes());
			if entry.value.len == 0 || matches!(entry.value.ty, BcfType::Null) {
				continue;
			}
			out.push(b'=');
			push_typed(out, &entry.value);
		}
		Ok(())
	}

	fn write_samples(&self, header: &Header, out: &mut Vec<u8>) -> Result<(), Error> {
		let samples = header.sample_count();
		if self.format_count() == 0 {
			out.extend_from_slice(b"\t.");
			for _ in 0..samples {
				out.extend_from_slice(b"\t.");
			}
			return Ok(());
		}
		let mut fields: Vec<(FormatEntry<'_>, bool)> = Vec::with_capacity(self.format_count());
		for entry in self.format() {
			let entry = entry?;
			if matches!(entry.values.ty, BcfType::Null) {
				return Err(Error::Bcf("FORMAT field without a type"));
			}
			let genotype = header.name(entry.key) == Some(GENOTYPE_KEY);
			out.push(if fields.is_empty() { b'\t' } else { b':' });
			out.extend_from_slice(tag_name(header, entry.key)?.as_bytes());
			fields.push((entry, genotype));
		}
		for sample in 0..samples {
			out.push(b'\t');
			for (index, (entry, genotype)) in fields.iter().enumerate() {
				if index > 0 {
					out.push(b':');
				}
				write_sample_field(header, entry, *genotype, sample, out)?;
			}
		}
		Ok(())
	}
}

fn push_typed(out: &mut Vec<u8>, value: &Typed<'_>) {
	match value.ty {
		BcfType::Null => {}
		BcfType::Int8 | BcfType::Int16 | BcfType::Int32 => push_int_list(out, value.ints()),
		BcfType::Float => push_float_list(out, value.floats()),
		BcfType::Char => push_text(out, value.text()),
	}
}

fn write_sample_field(
	header: &Header,
	entry: &FormatEntry<'_>,
	genotype: bool,
	sample: usize,
	out: &mut Vec<u8>,
) -> Result<(), Error> {
	let width = entry.values.ty.width();
	let step = entry
		.values
		.len
		.checked_mul(width)
		.ok_or(Error::Bcf("FORMAT field is too large"))?;
	let start = sample
		.checked_mul(step)
		.ok_or(Error::Bcf("FORMAT field is too large"))?;
	let bytes = entry
		.values
		.bytes
		.get(start..start + step)
		.ok_or(Error::Bcf("FORMAT field is truncated"))?;
	let values = Typed {
		ty: entry.values.ty,
		len: entry.values.len,
		bytes,
	};
	match values.ty {
		BcfType::Int8 | BcfType::Int16 | BcfType::Int32 if genotype => {
			write_genotype(out, values.ints(), header.version());
		}
		BcfType::Char => match values.text() {
			[byte] => out.push(if *byte == MISSING_STRING_BYTE {
				b'.'
			} else {
				*byte
			}),
			text => push_text(out, text),
		},
		_ => push_typed(out, &values),
	}
	Ok(())
}
