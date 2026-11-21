use super::{Field, QueryOutput, Segment};
use crate::Error;
use crate::filter::variant_type;
use crate::header::Header;
use crate::number::{push_float, push_int};
use crate::site::Site;
use crate::value::{Cells, END_INT, Loaded, MISSING_INT, is_end_float, is_missing_float};
use crate::vcf::{push_float_list, push_int_list, push_sample_text, write_genotype, write_info};

const MISSING: u8 = b'.';

pub(super) fn render<S: Site>(
	segments: &[Segment],
	header: &Header,
	samples: &[usize],
	site: &S,
	out: &mut QueryOutput,
) -> Result<(), Error> {
	for segment in segments {
		match segment {
			Segment::Literal(bytes) => out.bytes.extend_from_slice(bytes),
			Segment::Field(field) => emit(*field, header, site, None, out)?,
			Segment::PerSample(inner) => {
				for &sample in samples {
					for segment in inner {
						match segment {
							Segment::Literal(bytes) => out.bytes.extend_from_slice(bytes),
							Segment::Field(field) => emit(*field, header, site, Some(sample), out)?,
							Segment::PerSample(_) => {}
						}
					}
				}
			}
		}
	}
	Ok(())
}

fn format_columns<S: Site>(header: &Header, site: &S, out: &mut QueryOutput) -> Result<(), Error> {
	site.read_record(header, &mut out.record)?;
	out.line.clear();
	out.record.write_vcf(header, &mut out.line)?;
	let columns = out.line.splitn(9, |&byte| byte == b'\t').nth(8);
	out.bytes.extend_from_slice(columns.unwrap_or(&[MISSING]));
	Ok(())
}

fn emit<S: Site>(
	field: Field,
	header: &Header,
	site: &S,
	sample: Option<usize>,
	out: &mut QueryOutput,
) -> Result<(), Error> {
	match field {
		Field::Chrom => {
			let id = site.contig(header)?;
			let contig = header
				.contig(id)
				.ok_or(Error::Bcf("contig id is not in the header"))?;
			out.bytes.extend_from_slice(contig.name.as_bytes());
		}
		Field::Pos => push_int(&mut out.bytes, i64::from(site.start()?) + 1),
		Field::Pos0 => push_int(&mut out.bytes, i64::from(site.start()?)),
		Field::End => push_int(&mut out.bytes, i64::from(site.end(header)?)),
		Field::End0 => push_int(&mut out.bytes, i64::from(site.end(header)?) - 1),
		Field::Id => {
			let id = site.id()?;
			out.bytes
				.extend_from_slice(if id.is_empty() { &[MISSING] } else { id });
		}
		Field::Ref => {
			site.alleles(&mut out.alleles)?;
			out.bytes
				.extend_from_slice(out.alleles.get(0).unwrap_or_default());
		}
		Field::Alt(subscript) => alts(site, subscript, out)?,
		Field::FirstAlt => alts(site, Some(0), out)?,
		Field::Qual => match site.qual()? {
			Some(qual) if !is_missing_float(qual) => push_float(&mut out.bytes, qual),
			_ => out.bytes.push(MISSING),
		},
		Field::Filter => {
			site.filters(header, &mut out.keys)?;
			if out.keys.is_empty() {
				out.bytes.push(MISSING);
			}
			for (index, &key) in out.keys.iter().enumerate() {
				if index > 0 {
					out.bytes.push(b';');
				}
				out.bytes
					.extend_from_slice(header.name(key).unwrap_or_default().as_bytes());
			}
		}
		Field::InfoAll => {
			site.read_record(header, &mut out.record)?;
			write_info(&out.record, header, &mut out.bytes)?;
		}
		Field::FormatKeys => format_columns(header, site, out)?,
		Field::Info { key, subscript } => {
			let found = site.info(header, key, &mut out.loaded)?;
			info(&out.loaded, found, subscript, &mut out.bytes);
		}
		Field::Format { key, subscript } => {
			let found = site.format(header, key, &mut out.loaded)?;
			format(&out.loaded, found, sample, subscript, &mut out.bytes);
		}
		Field::Genotype | Field::TargetGenotype => {
			genotype(field == Field::TargetGenotype, header, site, sample, out)?;
		}
		Field::Sample => {
			let name = sample.and_then(|index| header.samples().get(index));
			out.bytes
				.extend_from_slice(name.map_or(&[][..], |name| name.as_bytes()));
		}
		Field::Type => {
			site.alleles(&mut out.alleles)?;
			let reference = out.alleles.get(0).unwrap_or_default();
			let mut present = [false; 7];
			for allele in out.alleles.iter().skip(1) {
				let name = variant_type(reference, allele);
				let slot = TYPES.iter().position(|(candidate, _)| *candidate == name);
				if let Some(flag) = slot.and_then(|slot| present.get_mut(slot)) {
					*flag = true;
				}
			}
			if out.alleles.len() <= 1
				&& let Some(flag) = present.first_mut()
			{
				*flag = true;
			}
			let mut first = true;
			for ((_, upper), _) in TYPES.iter().zip(present).filter(|(_, flag)| *flag) {
				if !first {
					out.bytes.push(b',');
				}
				first = false;
				out.bytes.extend_from_slice(upper.as_bytes());
			}
		}
		Field::Line => {
			site.read_record(header, &mut out.record)?;
			out.line.clear();
			out.record.write_vcf(header, &mut out.line)?;
			out.bytes.extend_from_slice(&out.line);
		}
	}
	Ok(())
}

const TYPES: [(&str, &str); 7] = [
	("ref", "REF"),
	("snp", "SNP"),
	("mnp", "MNP"),
	("indel", "INDEL"),
	("other", "OTHER"),
	("bnd", "BND"),
	("overlap", "OVERLAP"),
];

fn alts<S: Site>(site: &S, subscript: Option<usize>, out: &mut QueryOutput) -> Result<(), Error> {
	site.alleles(&mut out.alleles)?;
	if out.alleles.len() <= 1 {
		out.bytes.push(MISSING);
		return Ok(());
	}
	if let Some(index) = subscript {
		let allele = index.checked_add(1).and_then(|at| out.alleles.get(at));
		out.bytes.extend_from_slice(allele.unwrap_or(&[MISSING]));
		return Ok(());
	}
	for (index, allele) in out.alleles.iter().skip(1).enumerate() {
		if index > 0 {
			out.bytes.push(b',');
		}
		out.bytes.extend_from_slice(allele);
	}
	Ok(())
}

fn info(loaded: &Loaded, found: bool, subscript: Option<usize>, out: &mut Vec<u8>) {
	if !found {
		out.push(MISSING);
		return;
	}
	let width = loaded.width();
	match loaded.cells() {
		Cells::Flag => out.push(b'1'),
		Cells::Ints(_) if width == 1 => scalar_int(loaded.ints().first().copied(), out),
		Cells::Floats(_) if width == 1 => scalar_float(loaded.floats().first().copied(), out),
		Cells::Text(_) if width == 1 => {
			out.push(loaded.text().first().copied().unwrap_or(MISSING));
		}
		Cells::Ints(_) => match subscript {
			Some(index) => scalar_int(loaded.ints().get(index).copied(), out),
			None => array_ints(loaded.ints(), out),
		},
		Cells::Floats(_) => match subscript {
			Some(index) => scalar_float(loaded.floats().get(index).copied(), out),
			None => array_floats(loaded.floats(), out),
		},
		Cells::Text(_) => match subscript {
			Some(index) => comma_field(loaded.text(), index, out),
			None => out.extend_from_slice(loaded.text()),
		},
	}
}

fn format(
	loaded: &Loaded,
	found: bool,
	sample: Option<usize>,
	subscript: Option<usize>,
	out: &mut Vec<u8>,
) {
	let Some(sample) = sample.filter(|_| found) else {
		out.push(MISSING);
		return;
	};
	match loaded.cells() {
		Cells::Flag => out.push(MISSING),
		Cells::Ints(_) => {
			let values = loaded.sample_ints(sample);
			match subscript {
				Some(index) => scalar_int(values.get(index).copied(), out),
				None => array_ints(values, out),
			}
		}
		Cells::Floats(_) => {
			let values = loaded.sample_floats(sample);
			match subscript {
				Some(index) => scalar_float(values.get(index).copied(), out),
				None => array_floats(values, out),
			}
		}
		Cells::Text(_) => {
			let text = loaded.sample_text(sample);
			match subscript {
				Some(index) => {
					let text = text.split(|&byte| byte == 0).next().unwrap_or_default();
					comma_field(text, index, out);
				}
				None => push_sample_text(out, text),
			}
		}
	}
}

fn genotype<S: Site>(
	target: bool,
	header: &Header,
	site: &S,
	sample: Option<usize>,
	out: &mut QueryOutput,
) -> Result<(), Error> {
	let key = header
		.key("GT")
		.ok_or_else(|| Error::UndefinedTag("GT".to_owned()))?;
	let found = site.format(header, key, &mut out.loaded)?;
	let Some(sample) = sample.filter(|_| found) else {
		out.bytes.push(MISSING);
		return Ok(());
	};
	let values = out.loaded.sample_ints(sample);
	if !target {
		let before = out.bytes.len();
		write_genotype(&mut out.bytes, values.iter().copied(), header.version());
		if out.bytes.len() == before {
			out.bytes.push(MISSING);
		}
		return Ok(());
	}
	site.alleles(&mut out.alleles)?;
	let mut written = 0;
	for (position, &value) in values.iter().enumerate() {
		if value == END_INT {
			break;
		}
		if position > 0 {
			out.bytes.push(if value & 1 == 1 { b'|' } else { b'/' });
		}
		match usize::try_from((value >> 1) - 1)
			.ok()
			.filter(|_| value != MISSING_INT)
		{
			Some(index) => out
				.bytes
				.extend_from_slice(out.alleles.get(index).unwrap_or(&[MISSING])),
			None => out.bytes.push(MISSING),
		}
		written += 1;
	}
	if written == 0 {
		out.bytes.push(MISSING);
	}
	Ok(())
}

fn scalar_int(value: Option<i32>, out: &mut Vec<u8>) {
	match value {
		Some(value) if value != MISSING_INT && value != END_INT => push_int(out, i64::from(value)),
		_ => out.push(MISSING),
	}
}

fn scalar_float(value: Option<f32>, out: &mut Vec<u8>) {
	match value {
		Some(value) if !is_missing_float(value) && !is_end_float(value) => push_float(out, value),
		_ => out.push(MISSING),
	}
}

fn array_ints(values: &[i32], out: &mut Vec<u8>) {
	let before = out.len();
	push_int_list(out, values.iter().copied());
	if out.len() == before && !values.is_empty() {
		out.push(MISSING);
	}
}

fn array_floats(values: &[f32], out: &mut Vec<u8>) {
	let before = out.len();
	push_float_list(out, values.iter().copied());
	if out.len() == before && !values.is_empty() {
		out.push(MISSING);
	}
}

fn comma_field(text: &[u8], index: usize, out: &mut Vec<u8>) {
	match text.split(|&byte| byte == b',').nth(index) {
		Some(field) if !field.is_empty() => out.extend_from_slice(field),
		_ => out.push(MISSING),
	}
}
