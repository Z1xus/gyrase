#[derive(Debug)]
pub(super) struct Stream<'a> {
	data: &'a [u8],
	position: usize,
}

impl<'a> Stream<'a> {
	pub(super) fn new(data: &'a [u8]) -> Self {
		Self { data, position: 0 }
	}

	pub(super) fn byte(&mut self) -> u32 {
		let byte = self.data.get(self.position).copied().unwrap_or(0);
		self.position += 1;
		u32::from(byte)
	}

	pub(super) fn word(&mut self) -> u32 {
		let low = self.byte();
		low | (self.byte() << 8)
	}

	pub(super) fn is_overrun(&self) -> bool {
		self.position > self.data.len()
	}
}
