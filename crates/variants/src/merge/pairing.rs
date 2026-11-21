use std::collections::HashMap;

use super::input::Taken;
use crate::record::Record;

struct Var {
	tokens: Vec<String>,
	text: String,
	members: Vec<(usize, usize)>,
	groups: Vec<usize>,
}

struct VarSet {
	vars: Vec<usize>,
	groups: Vec<usize>,
	count: usize,
}

fn tokens(record: &Record) -> Vec<String> {
	let reference = String::from_utf8_lossy(record.allele(0).unwrap_or_default()).into_owned();
	let mut out: Vec<String> = record
		.alleles()
		.skip(1)
		.map(|alt| format!("{reference}>{}", String::from_utf8_lossy(alt)))
		.collect();
	if out.is_empty() {
		out.push(format!("{reference}>."));
	}
	out
}

fn same_set(left: &Var, right: &Var) -> bool {
	left.text == right.text
		|| (left.tokens.len() == right.tokens.len()
			&& left.text.len() == right.text.len()
			&& left.tokens.iter().all(|token| {
				right
					.tokens
					.iter()
					.any(|other| other.eq_ignore_ascii_case(token))
			}))
}

fn overlaps(left: &[usize], right: &[usize]) -> bool {
	left.iter().any(|group| right.contains(group))
}

fn pairable(vars: &[Var], left: &VarSet, right: &VarSet) -> bool {
	left.vars.iter().any(|&a| {
		right
			.vars
			.iter()
			.any(|&b| match (vars.get(a), vars.get(b)) {
				(Some(a), Some(b)) => same_set(a, b),
				_ => false,
			})
	})
}

fn register(
	known: &mut HashMap<String, usize>,
	vars: &mut Vec<Var>,
	tokens: Vec<String>,
	position: usize,
) -> usize {
	let text = tokens.join(",");
	let mut key = text.clone();
	let mut counter = 0;
	while let Some(&index) = known.get(&key) {
		let last = vars
			.get(index)
			.and_then(|var| var.members.last())
			.map(|member| member.0);
		if last != Some(position) {
			return index;
		}
		key = format!("{text}{counter}");
		counter += 1;
	}
	known.insert(key, vars.len());
	vars.push(Var {
		tokens,
		text,
		members: Vec::new(),
		groups: Vec::new(),
	});
	vars.len() - 1
}

fn collect(taken: &[Taken], ready: &[usize]) -> Vec<Var> {
	let mut vars: Vec<Var> = Vec::new();
	let mut known: HashMap<String, usize> = HashMap::new();
	let mut keys: HashMap<Vec<String>, usize> = HashMap::new();
	for (position, &input) in ready.iter().enumerate() {
		let Some(source) = taken.get(input) else {
			continue;
		};
		let mut texts = Vec::new();
		let mut indices = Vec::new();
		for slot in 0..source.count() {
			let Some(record) = source.record(slot) else {
				continue;
			};
			let index = register(&mut known, &mut vars, tokens(record), position);
			if let Some(var) = vars.get_mut(index) {
				var.members.push((position, slot));
				texts.push(var.text.clone());
			}
			indices.push(index);
		}
		texts.sort();
		let next = keys.len();
		let group = *keys.entry(texts).or_insert(next);
		for index in indices {
			if let Some(var) = vars.get_mut(index)
				&& !var.groups.contains(&group)
			{
				var.groups.push(group);
			}
		}
	}
	vars
}

pub(super) fn first_row(taken: &[Taken], ready: &[usize]) -> Vec<Option<usize>> {
	let mut picked = vec![None; ready.len()];
	if ready.len() < 2 {
		return picked;
	}
	let vars = collect(taken, ready);
	let mut sets: Vec<VarSet> = vars
		.iter()
		.enumerate()
		.map(|(index, var)| VarSet {
			vars: vec![index],
			groups: var.groups.clone(),
			count: var.members.len(),
		})
		.collect();
	while !sets.is_empty() {
		let mut best = 0;
		for (index, set) in sets.iter().enumerate() {
			if sets
				.get(best)
				.is_some_and(|current| current.count < set.count)
			{
				best = index;
			}
		}
		let Some(lead) = sets.get(best) else {
			break;
		};
		let partner = sets.iter().position(|other| {
			!overlaps(&lead.groups, &other.groups) && pairable(&vars, lead, other)
		});
		let Some(partner) = partner else {
			for var in lead.vars.iter().filter_map(|&index| vars.get(index)) {
				for &(position, slot) in &var.members {
					if let Some(entry) = picked.get_mut(position) {
						*entry = Some(slot);
					}
				}
			}
			break;
		};
		let (low, high) = (best.min(partner), best.max(partner));
		let removed = sets.remove(high);
		if let Some(target) = sets.get_mut(low) {
			target.vars.extend(removed.vars);
			target.groups.extend(removed.groups);
			target.count += removed.count;
		}
	}
	picked
}
