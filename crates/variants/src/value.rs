use std::ops::Range;

pub const MISSING_INT: i32 = i32::MIN;
pub const END_INT: i32 = i32::MIN + 1;
pub const MISSING_FLOAT_BITS: u32 = 0x7F80_0001;
pub const END_FLOAT_BITS: u32 = 0x7F80_0002;

pub const fn missing_float() -> f32 {
	f32::from_bits(MISSING_FLOAT_BITS)
}

pub const fn end_float() -> f32 {
	f32::from_bits(END_FLOAT_BITS)
}

pub const fn is_missing_float(value: f32) -> bool {
	value.to_bits() == MISSING_FLOAT_BITS
}

pub const fn is_end_float(value: f32) -> bool {
	value.to_bits() == END_FLOAT_BITS
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Number {
	Fixed(u32),
	PerAlt,
	PerAllele,
	PerGenotype,
	Variable,
}

impl Number {
	pub fn parse(text: &str) -> Option<Self> {
		match text {
			"A" => Some(Self::PerAlt),
			"R" => Some(Self::PerAllele),
			"G" => Some(Self::PerGenotype),
			"." => Some(Self::Variable),
			digits => digits.parse().ok().map(Self::Fixed),
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValueType {
	Flag,
	Integer,
	Float,
	Character,
	String,
}

impl ValueType {
	pub fn parse(text: &str) -> Option<Self> {
		match text {
			"Flag" => Some(Self::Flag),
			"Integer" => Some(Self::Integer),
			"Float" => Some(Self::Float),
			"Character" => Some(Self::Character),
			"String" => Some(Self::String),
			_ => None,
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BcfType {
	Null,
	Int8,
	Int16,
	Int32,
	Float,
	Char,
}

impl BcfType {
	pub const fn from_code(code: u8) -> Option<Self> {
		match code {
			0 => Some(Self::Null),
			1 => Some(Self::Int8),
			2 => Some(Self::Int16),
			3 => Some(Self::Int32),
			5 => Some(Self::Float),
			7 => Some(Self::Char),
			_ => None,
		}
	}

	pub const fn code(self) -> u8 {
		match self {
			Self::Null => 0,
			Self::Int8 => 1,
			Self::Int16 => 2,
			Self::Int32 => 3,
			Self::Float => 5,
			Self::Char => 7,
		}
	}

	pub const fn width(self) -> usize {
		match self {
			Self::Null | Self::Int8 | Self::Char => 1,
			Self::Int16 => 2,
			Self::Int32 | Self::Float => 4,
		}
	}
}

#[derive(Debug, Clone, Default)]
pub struct Store {
	pub(crate) ints: Vec<i32>,
	pub(crate) floats: Vec<f32>,
	pub(crate) text: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Cells {
	Flag,
	Ints(Range<usize>),
	Floats(Range<usize>),
	Text(Range<usize>),
}

impl Store {
	pub fn clear(&mut self) {
		self.ints.clear();
		self.floats.clear();
		self.text.clear();
	}

	pub fn push_ints(&mut self, values: impl IntoIterator<Item = i32>) -> Cells {
		let start = self.ints.len();
		self.ints.extend(values);
		Cells::Ints(start..self.ints.len())
	}

	pub fn push_floats(&mut self, values: impl IntoIterator<Item = f32>) -> Cells {
		let start = self.floats.len();
		self.floats.extend(values);
		Cells::Floats(start..self.floats.len())
	}

	pub fn push_text(&mut self, bytes: &[u8]) -> Cells {
		let start = self.text.len();
		self.text.extend_from_slice(bytes);
		Cells::Text(start..self.text.len())
	}

	pub fn ints(&self, cells: &Cells) -> &[i32] {
		match cells {
			Cells::Ints(range) => self.ints.get(range.clone()).unwrap_or_default(),
			_ => &[],
		}
	}

	pub fn floats(&self, cells: &Cells) -> &[f32] {
		match cells {
			Cells::Floats(range) => self.floats.get(range.clone()).unwrap_or_default(),
			_ => &[],
		}
	}

	pub fn text(&self, cells: &Cells) -> &[u8] {
		match cells {
			Cells::Text(range) => self.text.get(range.clone()).unwrap_or_default(),
			_ => &[],
		}
	}
}

#[derive(Debug, Clone)]
pub struct Loaded {
	pub(crate) store: Store,
	pub(crate) cells: Cells,
	pub(crate) width: usize,
}

impl Default for Loaded {
	fn default() -> Self {
		Self {
			store: Store::default(),
			cells: Cells::Flag,
			width: 0,
		}
	}
}

impl Loaded {
	pub fn clear(&mut self) {
		self.store.clear();
		self.cells = Cells::Flag;
		self.width = 0;
	}

	pub fn width(&self) -> usize {
		self.width
	}

	pub fn cells(&self) -> &Cells {
		&self.cells
	}

	pub fn ints(&self) -> &[i32] {
		self.store.ints(&self.cells)
	}

	pub fn floats(&self) -> &[f32] {
		self.store.floats(&self.cells)
	}

	pub fn text(&self) -> &[u8] {
		self.store.text(&self.cells)
	}

	pub fn sample_ints(&self, sample: usize) -> &[i32] {
		slice_of(self.ints(), sample, self.width)
	}

	pub fn sample_floats(&self, sample: usize) -> &[f32] {
		slice_of(self.floats(), sample, self.width)
	}

	pub fn sample_text(&self, sample: usize) -> &[u8] {
		slice_of(self.text(), sample, self.width)
	}

	pub fn samples(&self) -> usize {
		let len = match &self.cells {
			Cells::Flag => return 0,
			Cells::Ints(range) | Cells::Floats(range) | Cells::Text(range) => range.len(),
		};
		len.checked_div(self.width).unwrap_or(0)
	}
}

fn slice_of<T>(all: &[T], sample: usize, width: usize) -> &[T] {
	sample
		.checked_mul(width)
		.and_then(|start| all.get(start..start.checked_add(width)?))
		.unwrap_or_default()
}

pub(crate) trait Numeric: Copy {
	const MISSING: Self;
	const END: Self;
	fn is_end(self) -> bool;
	fn is_missing(self) -> bool;
}

impl Numeric for i32 {
	const MISSING: Self = MISSING_INT;
	const END: Self = END_INT;

	fn is_end(self) -> bool {
		self == END_INT
	}

	fn is_missing(self) -> bool {
		self == MISSING_INT
	}
}

impl Numeric for f32 {
	const MISSING: Self = missing_float();
	const END: Self = end_float();

	fn is_end(self) -> bool {
		is_end_float(self)
	}

	fn is_missing(self) -> bool {
		is_missing_float(self)
	}
}
