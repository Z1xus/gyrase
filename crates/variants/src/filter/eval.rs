use super::classify::variant_type;
use super::lex::Comparison;
use super::node::{Arithmetic, Function, GenotypePattern, Logic, Node, Select, SiteField};
use super::value::{Shape, Value};
use crate::Error;
use crate::header::Header;
use crate::site::{Alleles, Site};
use crate::value::{Cells, END_INT, Loaded, MISSING_INT, ValueType, is_end_float};

const MISSING_TEXT: &[u8] = b".";
const VARIANT_TYPES: [&str; 7] = ["ref", "snp", "mnp", "indel", "other", "bnd", "overlap"];
const NOT_PRESENT: Error = Error::Expression {
	position: 0,
	reason: "value is not available",
};

#[derive(Debug, Default)]
pub struct Scratch {
	pool: Vec<Value>,
	spare: Value,
	alleles: Alleles,
	keys: Vec<u32>,
	loaded: Loaded,
	text: Vec<u8>,
}

struct Context<'a, S> {
	site: &'a S,
	header: &'a Header,
}

pub(super) fn evaluate<S: Site>(
	node: &Node,
	site: &S,
	header: &Header,
	scratch: &mut Scratch,
) -> Result<bool, Error> {
	let context = Context { site, header };
	eval(&context, node, scratch, 0)?;
	Ok(scratch.pool.first().is_some_and(Value::truth))
}

fn eval<S: Site>(
	context: &Context<'_, S>,
	node: &Node,
	scratch: &mut Scratch,
	depth: usize,
) -> Result<(), Error> {
	while scratch.pool.len() < depth + 2 {
		scratch.pool.push(Value::default());
	}
	match node {
		Node::Number(_) | Node::Missing | Node::Text(_) => constant(node, scratch, depth)?,
		Node::Site(field) => load_site(context, *field, scratch, depth)?,
		Node::Info {
			key,
			value_type,
			elements,
		} => {
			load_info(context, *key, *value_type, *elements, scratch, depth)?;
		}
		Node::Format {
			key,
			value_type,
			samples,
			elements,
		} => {
			load_format(
				context,
				*key,
				*value_type,
				*samples,
				*elements,
				scratch,
				depth,
			)?;
		}
		Node::Genotype {
			pattern,
			negated,
			samples,
		} => {
			genotype(context, pattern, *negated, *samples, scratch, depth)?;
		}
		Node::IsMissing { node, negated } => {
			eval(context, node, scratch, depth)?;
			is_missing(*negated, scratch, depth)?;
		}
		Node::Compare {
			comparison,
			left,
			right,
		} => {
			eval(context, left, scratch, depth)?;
			eval(context, right, scratch, depth + 1)?;
			compare(
				*comparison,
				left.kind() == super::node::Kind::Text,
				scratch,
				depth,
			)?;
		}
		Node::Matches {
			node,
			regex,
			negated,
		} => {
			eval(context, node, scratch, depth)?;
			matches_regex(regex, *negated, scratch, depth)?;
		}
		Node::Not(node) => {
			eval(context, node, scratch, depth)?;
			for flag in &mut slot(scratch, depth)?.flags {
				*flag = !*flag;
			}
		}
		Node::Negate(node) => {
			eval(context, node, scratch, depth)?;
			for number in &mut slot(scratch, depth)?.numbers {
				*number = -*number;
			}
		}
		Node::Arithmetic {
			operator,
			left,
			right,
		} => {
			eval(context, left, scratch, depth)?;
			eval(context, right, scratch, depth + 1)?;
			arithmetic(*operator, scratch, depth)?;
		}
		Node::Logic {
			operator,
			left,
			right,
		} => {
			eval(context, left, scratch, depth)?;
			eval(context, right, scratch, depth + 1)?;
			logic(*operator, scratch, depth)?;
		}
		Node::Call { function, argument } => {
			eval(context, argument, scratch, depth)?;
			call(*function, scratch, depth)?;
		}
	}
	Ok(())
}

fn slot(scratch: &mut Scratch, depth: usize) -> Result<&mut Value, Error> {
	scratch.pool.get_mut(depth).ok_or(NOT_PRESENT)
}

fn finish(scratch: &mut Scratch, depth: usize) -> Result<(), Error> {
	let target = scratch.pool.get_mut(depth).ok_or(NOT_PRESENT)?;
	std::mem::swap(target, &mut scratch.spare);
	Ok(())
}

fn load_site<S: Site>(
	context: &Context<'_, S>,
	field: SiteField,
	scratch: &mut Scratch,
	depth: usize,
) -> Result<(), Error> {
	let Scratch {
		pool,
		alleles,
		keys,
		text,
		..
	} = scratch;
	let out = pool.get_mut(depth).ok_or(NOT_PRESENT)?;
	out.reset(Shape::Site);
	match field {
		SiteField::Chrom => {
			let id = context.site.contig(context.header)?;
			let contig = context
				.header
				.contig(id)
				.ok_or(Error::Bcf("contig id is not in the header"))?;
			out.push_text(contig.name.as_bytes());
		}
		SiteField::Pos => out.numbers.push(f64::from(context.site.start()?) + 1.0),
		SiteField::Id => {
			let id = context.site.id()?;
			out.push_text(if id.is_empty() { MISSING_TEXT } else { id });
		}
		SiteField::Ref => {
			context.site.alleles(alleles)?;
			out.push_text(alleles.get(0).unwrap_or_default());
		}
		SiteField::Alt => {
			context.site.alleles(alleles)?;
			for allele in alleles.iter().skip(1) {
				out.push_text(allele);
			}
			if out.ends.is_empty() {
				out.push_text(MISSING_TEXT);
			}
		}
		SiteField::Qual => out
			.numbers
			.push(context.site.qual()?.map_or(f64::NAN, f64::from)),
		SiteField::Filter => {
			context.site.filters(context.header, keys)?;
			for &key in keys.iter() {
				out.push_text(context.header.name(key).unwrap_or_default().as_bytes());
			}
			if out.ends.is_empty() {
				out.push_text(MISSING_TEXT);
			}
		}
		SiteField::AlternateCount => {
			context.site.alleles(alleles)?;
			out.numbers
				.push(count_as_float(alleles.len().saturating_sub(1)));
		}
		SiteField::SampleCount => out
			.numbers
			.push(count_as_float(context.header.sample_count())),
		SiteField::Type => {
			context.site.alleles(alleles)?;
			let reference = alleles.get(0).unwrap_or_default();
			let mut present = [false; VARIANT_TYPES.len()];
			for allele in alleles.iter().skip(1) {
				let name = variant_type(reference, allele);
				let slot = VARIANT_TYPES
					.iter()
					.position(|candidate| *candidate == name);
				if let Some(flag) = slot.and_then(|slot| present.get_mut(slot)) {
					*flag = true;
				}
			}
			if alleles.len() <= 1 {
				present[0] = true;
			}
			let combined = text;
			combined.clear();
			for (name, _) in VARIANT_TYPES.iter().zip(present).filter(|(_, flag)| *flag) {
				if !combined.is_empty() {
					combined.push(b',');
				}
				combined.extend_from_slice(name.as_bytes());
			}
			out.push_text(combined);
		}
	}
	Ok(())
}

fn push_numbers(loaded: &Loaded, selected: Select, out: &mut Vec<f64>, keep_width: usize) {
	let start = out.len();
	match loaded.cells() {
		Cells::Ints(_) => {
			for (index, &value) in loaded.ints().iter().enumerate() {
				if value == END_INT {
					continue;
				}
				out.push(
					if value == MISSING_INT || !selected.contains(index % keep_width) {
						f64::NAN
					} else {
						f64::from(value)
					},
				);
			}
		}
		Cells::Floats(_) => {
			for (index, &value) in loaded.floats().iter().enumerate() {
				if is_end_float(value) {
					continue;
				}
				out.push(if selected.contains(index % keep_width) {
					f64::from(value)
				} else {
					f64::NAN
				});
			}
		}
		_ => {}
	}
	let _ = start;
}

fn load_info<S: Site>(
	context: &Context<'_, S>,
	key: u32,
	value_type: ValueType,
	elements: Select,
	scratch: &mut Scratch,
	depth: usize,
) -> Result<(), Error> {
	let Scratch { pool, loaded, .. } = scratch;
	let out = pool.get_mut(depth).ok_or(NOT_PRESENT)?;
	out.reset(Shape::Site);
	if !context.site.info(context.header, key, loaded)? {
		return Ok(());
	}
	match value_type {
		ValueType::Flag => out.numbers.push(1.0),
		ValueType::Integer | ValueType::Float => {
			let width = loaded.width().max(1);
			push_numbers(
				loaded,
				elements,
				&mut out.numbers,
				width.max(loaded.ints().len().max(loaded.floats().len())),
			);
		}
		ValueType::String | ValueType::Character => {
			for (index, item) in loaded.text().split(|&byte| byte == b',').enumerate() {
				if elements.contains(index) {
					out.push_text(item);
				}
			}
		}
	}
	Ok(())
}

fn load_format<S: Site>(
	context: &Context<'_, S>,
	key: u32,
	value_type: ValueType,
	samples: Select,
	elements: Select,
	scratch: &mut Scratch,
	depth: usize,
) -> Result<(), Error> {
	let Scratch { pool, loaded, .. } = scratch;
	let out = pool.get_mut(depth).ok_or(NOT_PRESENT)?;
	out.reset(Shape::Samples);
	if !context.site.format(context.header, key, loaded)? {
		return Ok(());
	}
	let width = loaded.width().max(1);
	out.width = width;
	match value_type {
		ValueType::Integer | ValueType::Float | ValueType::Flag => {
			push_numbers(loaded, elements, &mut out.numbers, width);
			if out.numbers.len() != loaded.ints().len().max(loaded.floats().len()) {
				return Err(Error::Bcf("FORMAT vector has an end marker in the middle"));
			}
			for (index, number) in out.numbers.iter_mut().enumerate() {
				if !samples.contains(index / width) {
					*number = f64::NAN;
				}
			}
		}
		ValueType::String | ValueType::Character => {
			for sample in 0..loaded.samples() {
				let text = loaded.sample_text(sample);
				let trimmed = text.split(|&byte| byte == 0).next().unwrap_or_default();
				out.push_text(if samples.contains(sample) {
					trimmed
				} else {
					&[]
				});
			}
		}
	}
	Ok(())
}

fn genotype<S: Site>(
	context: &Context<'_, S>,
	pattern: &GenotypePattern,
	negated: bool,
	samples: Select,
	scratch: &mut Scratch,
	depth: usize,
) -> Result<(), Error> {
	let Scratch {
		pool, loaded, text, ..
	} = scratch;
	let out = pool.get_mut(depth).ok_or(NOT_PRESENT)?;
	out.reset(Shape::Samples);
	let key = context
		.header
		.key("GT")
		.ok_or_else(|| Error::UndefinedTag("GT".to_owned()))?;
	if !context.site.format(context.header, key, loaded)? {
		return Ok(());
	}
	for sample in 0..loaded.samples() {
		let alleles = loaded.sample_ints(sample);
		let matched = samples.contains(sample) && genotype_matches(alleles, pattern, text);
		out.flags.push(matched != negated);
	}
	Ok(())
}

fn genotype_matches(alleles: &[i32], pattern: &GenotypePattern, text: &mut Vec<u8>) -> bool {
	let called: Vec<i32> = alleles
		.iter()
		.copied()
		.take_while(|&value| value != END_INT)
		.collect();
	let indices = called.iter().map(|&value| {
		if value == MISSING_INT {
			-1
		} else {
			(value >> 1) - 1
		}
	});
	let any_missing = indices.clone().any(|index| index < 0);
	let first = indices.clone().next();
	let all_equal = indices.clone().all(|index| Some(index) == first);
	let ploidy = called.len();
	match pattern {
		GenotypePattern::Missing => any_missing || ploidy == 0,
		GenotypePattern::Reference => ploidy > 0 && indices.clone().all(|index| index == 0),
		GenotypePattern::Alternate => indices.clone().any(|index| index > 0),
		GenotypePattern::Homozygous => ploidy > 1 && !any_missing && all_equal,
		GenotypePattern::Heterozygous => ploidy > 1 && !any_missing && !all_equal,
		GenotypePattern::Haploid => ploidy == 1,
		GenotypePattern::Alleles(expected) => {
			text.clear();
			for (position, (&value, index)) in called.iter().zip(indices).enumerate() {
				if position > 0 {
					text.push(if value & 1 == 1 { b'|' } else { b'/' });
				}
				if index < 0 {
					text.push(b'.');
				} else {
					text.extend_from_slice(index.to_string().as_bytes());
				}
			}
			text == expected
		}
	}
}

fn is_missing(negated: bool, scratch: &mut Scratch, depth: usize) -> Result<(), Error> {
	let source = scratch.pool.get(depth).ok_or(NOT_PRESENT)?;
	let spare = &mut scratch.spare;
	spare.reset(source.shape);
	let missing = |value: &Value, sample: usize| {
		let numbers = value
			.number_items(sample)
			.iter()
			.all(|number| number.is_nan());
		let range = value.item_range(sample);
		let texts = range.clone().all(|index| {
			let text = value.item(index);
			text.is_empty() || text == MISSING_TEXT
		});
		if value.ends.is_empty() {
			numbers
		} else {
			texts
		}
	};
	match source.shape {
		Shape::Site => spare.flags.push(missing(source, 0) != negated),
		Shape::Samples => {
			for sample in 0..source.sample_count() {
				spare.flags.push(missing(source, sample) != negated);
			}
		}
	}
	finish(scratch, depth)
}

#[allow(clippy::float_cmp)]
fn number_test(comparison: Comparison, left: f64, right: f64) -> bool {
	match comparison {
		Comparison::Equal => left == right,
		Comparison::NotEqual => !left.is_nan() && !right.is_nan() && left != right,
		Comparison::Less => left < right,
		Comparison::LessOrEqual => left <= right,
		Comparison::Greater => left > right,
		Comparison::GreaterOrEqual => left >= right,
		Comparison::Match | Comparison::NotMatch => false,
	}
}

fn compare(
	comparison: Comparison,
	text: bool,
	scratch: &mut Scratch,
	depth: usize,
) -> Result<(), Error> {
	let (left, right) = pair(&scratch.pool, depth)?;
	let spare = &mut scratch.spare;
	let shared = left.sample_count().max(right.sample_count());
	let shape = if shared == 0 && left.shape == Shape::Site && right.shape == Shape::Site {
		Shape::Site
	} else {
		Shape::Samples
	};
	spare.reset(shape);
	let count = if shape == Shape::Site { 1 } else { shared };
	for sample in 0..count {
		let result = if text {
			text_test(comparison, left, right, sample)
		} else {
			let lefts = left.number_items(sample);
			let rights = right.number_items(sample);
			lefts
				.iter()
				.any(|&l| rights.iter().any(|&r| number_test(comparison, l, r)))
		};
		spare.flags.push(result);
	}
	finish(scratch, depth)
}

fn text_test(comparison: Comparison, left: &Value, right: &Value, sample: usize) -> bool {
	let lefts = left.item_range(sample);
	let rights = right.item_range(sample);
	if lefts.is_empty() || rights.is_empty() {
		return false;
	}
	let equal = |l: usize| rights.clone().any(|r| left.item(l) == right.item(r));
	match comparison {
		Comparison::Equal => lefts.clone().any(equal),
		Comparison::NotEqual => !lefts.clone().any(equal),
		_ => false,
	}
}

fn matches_regex(
	regex: &super::regex::Regex,
	negated: bool,
	scratch: &mut Scratch,
	depth: usize,
) -> Result<(), Error> {
	let source = scratch.pool.get(depth).ok_or(NOT_PRESENT)?;
	let spare = &mut scratch.spare;
	spare.reset(source.shape);
	let test = |sample: usize| {
		let range = source.item_range(sample);
		!range.is_empty()
			&& range
				.clone()
				.any(|index| regex.is_match(source.item(index)))
				!= negated
	};
	match source.shape {
		Shape::Site => spare.flags.push(test(0)),
		Shape::Samples => {
			for sample in 0..source.sample_count() {
				spare.flags.push(test(sample));
			}
		}
	}
	finish(scratch, depth)
}

fn pair(pool: &[Value], depth: usize) -> Result<(&Value, &Value), Error> {
	match pool.get(depth..=depth + 1) {
		Some([left, right]) => Ok((left, right)),
		_ => Err(NOT_PRESENT),
	}
}

fn arithmetic(operator: Arithmetic, scratch: &mut Scratch, depth: usize) -> Result<(), Error> {
	let (left, right) = pair(&scratch.pool, depth)?;
	let spare = &mut scratch.spare;
	let shape = if left.shape == Shape::Samples {
		left.shape
	} else {
		right.shape
	};
	spare.reset(shape);
	spare.width = if left.shape == Shape::Samples {
		left.width
	} else {
		right.width
	};
	let length = if left.numbers.is_empty() || right.numbers.is_empty() {
		0
	} else {
		left.numbers.len().max(right.numbers.len())
	};
	for index in 0..length {
		let a = left
			.numbers
			.get(index % left.numbers.len())
			.copied()
			.unwrap_or(f64::NAN);
		let b = right
			.numbers
			.get(index % right.numbers.len())
			.copied()
			.unwrap_or(f64::NAN);
		spare.numbers.push(match operator {
			Arithmetic::Add => a + b,
			Arithmetic::Subtract => a - b,
			Arithmetic::Multiply => a * b,
			Arithmetic::Divide => a / b,
		});
	}
	finish(scratch, depth)
}

fn logic(operator: Logic, scratch: &mut Scratch, depth: usize) -> Result<(), Error> {
	let (left, right) = pair(&scratch.pool, depth)?;
	let spare = &mut scratch.spare;
	let and = matches!(operator, Logic::AndSite | Logic::AndSample);
	let combine = |a: bool, b: bool| if and { a && b } else { a || b };
	let per_sample = matches!(operator, Logic::AndSample | Logic::OrSample)
		&& (left.shape == Shape::Samples || right.shape == Shape::Samples);
	if per_sample {
		spare.reset(Shape::Samples);
		let count = left.sample_count().max(right.sample_count());
		for sample in 0..count {
			let a = flag_at(left, sample);
			let b = flag_at(right, sample);
			spare.flags.push(combine(a, b));
		}
	} else {
		spare.reset(Shape::Site);
		spare.flags.push(combine(left.truth(), right.truth()));
	}
	finish(scratch, depth)
}

fn flag_at(value: &Value, sample: usize) -> bool {
	match value.shape {
		Shape::Site => value.truth(),
		Shape::Samples => value.flags.get(sample).copied().unwrap_or(false),
	}
}

fn call(function: Function, scratch: &mut Scratch, depth: usize) -> Result<(), Error> {
	let source = scratch.pool.get(depth).ok_or(NOT_PRESENT)?;
	let spare = &mut scratch.spare;
	spare.reset(Shape::Site);
	let values = source
		.numbers
		.iter()
		.copied()
		.filter(|number| !number.is_nan());
	match function {
		Function::Abs => {
			spare.shape = source.shape;
			spare.width = source.width;
			spare
				.numbers
				.extend(source.numbers.iter().map(|number| number.abs()));
		}
		Function::Count => spare.numbers.push(count_as_float(values.count())),
		Function::Max => spare.numbers.push(values.fold(f64::NAN, f64::max)),
		Function::Min => spare.numbers.push(values.fold(f64::NAN, f64::min)),
		Function::Sum => {
			let (sum, count) = values.fold((0.0, 0usize), |(sum, count), value| {
				(sum + value, count + 1)
			});
			spare.numbers.push(if count == 0 { f64::NAN } else { sum });
		}
		Function::Average => {
			let (sum, count) = values.fold((0.0, 0usize), |(sum, count), value| {
				(sum + value, count + 1)
			});
			spare.numbers.push(if count == 0 {
				f64::NAN
			} else {
				sum / count_as_float(count)
			});
		}
	}
	finish(scratch, depth)
}

fn count_as_float(count: usize) -> f64 {
	f64::from(u32::try_from(count).unwrap_or(u32::MAX))
}

fn constant(node: &Node, scratch: &mut Scratch, depth: usize) -> Result<(), Error> {
	let slot = slot(scratch, depth)?;
	slot.reset(Shape::Site);
	match node {
		Node::Number(number) => slot.numbers.push(*number),
		Node::Text(text) => slot.push_text(text),
		_ => slot.numbers.push(f64::NAN),
	}
	Ok(())
}
