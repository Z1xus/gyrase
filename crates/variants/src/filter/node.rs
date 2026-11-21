use super::regex::Regex;
use crate::value::ValueType;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Select {
	All,
	Range(usize, usize),
}

impl Select {
	pub(super) fn contains(self, index: usize) -> bool {
		match self {
			Self::All => true,
			Self::Range(first, last) => first <= index && index <= last,
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SiteField {
	Chrom,
	Pos,
	Id,
	Ref,
	Alt,
	Qual,
	Filter,
	Type,
	AlternateCount,
	SampleCount,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Arithmetic {
	Add,
	Subtract,
	Multiply,
	Divide,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Logic {
	AndSite,
	AndSample,
	OrSite,
	OrSample,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Function {
	Max,
	Min,
	Sum,
	Average,
	Abs,
	Count,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum GenotypePattern {
	Missing,
	Reference,
	Alternate,
	Homozygous,
	Heterozygous,
	Haploid,
	Alleles(Vec<u8>),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Kind {
	Numbers,
	Text,
	Truth,
}

#[derive(Debug, Clone)]
pub(super) enum Node {
	Number(f64),
	Missing,
	Text(Vec<u8>),
	Site(SiteField),
	Info {
		key: u32,
		value_type: ValueType,
		elements: Select,
	},
	Format {
		key: u32,
		value_type: ValueType,
		samples: Select,
		elements: Select,
	},
	Genotype {
		pattern: GenotypePattern,
		negated: bool,
		samples: Select,
	},
	IsMissing {
		node: Box<Node>,
		negated: bool,
	},
	Compare {
		comparison: super::lex::Comparison,
		left: Box<Node>,
		right: Box<Node>,
	},
	Matches {
		node: Box<Node>,
		regex: Regex,
		negated: bool,
	},
	Not(Box<Node>),
	Negate(Box<Node>),
	Arithmetic {
		operator: Arithmetic,
		left: Box<Node>,
		right: Box<Node>,
	},
	Logic {
		operator: Logic,
		left: Box<Node>,
		right: Box<Node>,
	},
	Call {
		function: Function,
		argument: Box<Node>,
	},
}

impl Node {
	pub(super) fn kind(&self) -> Kind {
		match self {
			Self::Number(_)
			| Self::Missing
			| Self::Negate(_)
			| Self::Arithmetic { .. }
			| Self::Call { .. } => Kind::Numbers,
			Self::Text(_) => Kind::Text,
			Self::Site(field) => match field {
				SiteField::Pos
				| SiteField::Qual
				| SiteField::AlternateCount
				| SiteField::SampleCount => Kind::Numbers,
				_ => Kind::Text,
			},
			Self::Info { value_type, .. } | Self::Format { value_type, .. } => match value_type {
				ValueType::Integer | ValueType::Float | ValueType::Flag => Kind::Numbers,
				ValueType::String | ValueType::Character => Kind::Text,
			},
			Self::Genotype { .. }
			| Self::IsMissing { .. }
			| Self::Compare { .. }
			| Self::Matches { .. }
			| Self::Not(_)
			| Self::Logic { .. } => Kind::Truth,
		}
	}
}
