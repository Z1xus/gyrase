use std::ops::Range;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(super) enum Shape {
	#[default]
	Site,
	Samples,
}

#[derive(Debug, Default, Clone)]
pub(super) struct Value {
	pub shape: Shape,
	pub width: usize,
	pub numbers: Vec<f64>,
	pub bytes: Vec<u8>,
	pub ends: Vec<usize>,
	pub flags: Vec<bool>,
}

impl Value {
	pub fn reset(&mut self, shape: Shape) {
		self.shape = shape;
		self.width = 1;
		self.numbers.clear();
		self.bytes.clear();
		self.ends.clear();
		self.flags.clear();
	}

	pub fn push_text(&mut self, text: &[u8]) {
		self.bytes.extend_from_slice(text);
		self.ends.push(self.bytes.len());
	}

	pub fn sample_count(&self) -> usize {
		match self.shape {
			Shape::Site => 0,
			Shape::Samples => self
				.flags
				.len()
				.max(self.ends.len())
				.max(self.numbers.len().checked_div(self.width).unwrap_or(0)),
		}
	}

	pub fn number_items(&self, sample: usize) -> &[f64] {
		match self.shape {
			Shape::Site => &self.numbers,
			Shape::Samples => sample
				.checked_mul(self.width)
				.and_then(|start| self.numbers.get(start..start.checked_add(self.width)?))
				.unwrap_or_default(),
		}
	}

	pub fn item_range(&self, sample: usize) -> Range<usize> {
		match self.shape {
			Shape::Site => 0..self.ends.len(),
			Shape::Samples if sample < self.ends.len() => sample..sample + 1,
			Shape::Samples => 0..0,
		}
	}

	pub fn item(&self, index: usize) -> &[u8] {
		let start = match index.checked_sub(1) {
			Some(previous) => self.ends.get(previous).copied().unwrap_or(0),
			None => 0,
		};
		let end = self.ends.get(index).copied().unwrap_or(0);
		self.bytes.get(start..end).unwrap_or_default()
	}

	pub fn truth(&self) -> bool {
		self.flags.iter().any(|&flag| flag)
	}
}
