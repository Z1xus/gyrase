#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Field {
	pub tag: [u8; 2],
	pub value: String,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Fields {
	entries: Vec<Field>,
}

impl Fields {
	pub fn new() -> Self {
		Self::default()
	}

	pub fn with(mut self, tag: [u8; 2], value: impl Into<String>) -> Self {
		self.entries.push(Field {
			tag,
			value: value.into(),
		});
		self
	}

	pub fn get(&self, tag: [u8; 2]) -> Option<&str> {
		self.entries
			.iter()
			.find(|field| field.tag == tag)
			.map(|field| field.value.as_str())
	}

	pub fn iter(&self) -> impl Iterator<Item = &Field> {
		self.entries.iter()
	}

	pub(crate) fn push(&mut self, field: Field) {
		self.entries.push(field);
	}
}
