#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Chrom(u8);

const LABELS: [&str; Chrom::COUNT] = [
	"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
	"18", "19", "20", "21", "22", "X", "Y", "MT",
];

impl Chrom {
	pub const COUNT: usize = 25;
	pub const X: Self = Self(22);
	pub const Y: Self = Self(23);
	pub const MT: Self = Self(24);

	pub fn from_index(index: usize) -> Option<Self> {
		if index < Self::COUNT {
			u8::try_from(index).ok().map(Self)
		} else {
			None
		}
	}

	pub fn parse(name: &str) -> Option<Self> {
		let name = name.strip_prefix("chr").unwrap_or(name);
		let name = if name == "M" { "MT" } else { name };
		LABELS
			.iter()
			.position(|label| *label == name)
			.and_then(Self::from_index)
	}

	pub const fn index(self) -> usize {
		self.0 as usize
	}

	pub fn label(self) -> &'static str {
		LABELS.get(self.index()).copied().unwrap_or_default()
	}

	pub fn autosome_number(self) -> Option<u8> {
		(self.0 < Self::X.0).then_some(self.0 + 1)
	}
}
