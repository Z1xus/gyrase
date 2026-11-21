use std::collections::BTreeMap;

use super::layout::put_itf8;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Series {
	Flags,
	CompressionFlags,
	ReadLength,
	Position,
	ReadGroup,
	Name,
	MateFlags,
	MateReference,
	MatePosition,
	TemplateLength,
	TagLine,
	FeatureCount,
	FeatureCode,
	FeaturePosition,
	DeletionLength,
	Base,
	Substitution,
	Insertion,
	SoftClip,
	HardClip,
	Padding,
	ReferenceSkip,
	MappingQuality,
	Bases,
	Quality,
}

const COUNT: usize = 25;

impl Series {
	pub(super) const ALL: [Self; COUNT] = [
		Self::Flags,
		Self::CompressionFlags,
		Self::ReadLength,
		Self::Position,
		Self::ReadGroup,
		Self::Name,
		Self::MateFlags,
		Self::MateReference,
		Self::MatePosition,
		Self::TemplateLength,
		Self::TagLine,
		Self::FeatureCount,
		Self::FeatureCode,
		Self::FeaturePosition,
		Self::DeletionLength,
		Self::Base,
		Self::Substitution,
		Self::Insertion,
		Self::SoftClip,
		Self::HardClip,
		Self::Padding,
		Self::ReferenceSkip,
		Self::MappingQuality,
		Self::Bases,
		Self::Quality,
	];

	pub(super) fn key(self) -> [u8; 2] {
		match self {
			Self::Flags => *b"BF",
			Self::CompressionFlags => *b"CF",
			Self::ReadLength => *b"RL",
			Self::Position => *b"AP",
			Self::ReadGroup => *b"RG",
			Self::Name => *b"RN",
			Self::MateFlags => *b"MF",
			Self::MateReference => *b"NS",
			Self::MatePosition => *b"NP",
			Self::TemplateLength => *b"TS",
			Self::TagLine => *b"TL",
			Self::FeatureCount => *b"FN",
			Self::FeatureCode => *b"FC",
			Self::FeaturePosition => *b"FP",
			Self::DeletionLength => *b"DL",
			Self::Base => *b"BA",
			Self::Substitution => *b"BS",
			Self::Insertion => *b"IN",
			Self::SoftClip => *b"SC",
			Self::HardClip => *b"HC",
			Self::Padding => *b"PD",
			Self::ReferenceSkip => *b"RS",
			Self::MappingQuality => *b"MQ",
			Self::Bases => *b"BB",
			Self::Quality => *b"QS",
		}
	}

	pub(super) fn content_id(self) -> i32 {
		self as i32 + 1
	}

	pub(super) fn is_terminated(self) -> bool {
		matches!(
			self,
			Self::Name | Self::Insertion | Self::SoftClip | Self::Bases
		)
	}
}

#[derive(Debug)]
pub(super) struct Streams {
	series: [Vec<u8>; COUNT],
	tags: BTreeMap<i32, Vec<u8>>,
}

impl Streams {
	pub(super) fn new() -> Self {
		Self {
			series: std::array::from_fn(|_| Vec::new()),
			tags: BTreeMap::new(),
		}
	}

	fn stream(&mut self, series: Series) -> &mut Vec<u8> {
		&mut self.series[series as usize]
	}

	pub(super) fn int(&mut self, series: Series, value: i32) {
		put_itf8(self.stream(series), value);
	}

	pub(super) fn byte(&mut self, series: Series, value: u8) {
		self.stream(series).push(value);
	}

	pub(super) fn bytes(&mut self, series: Series, values: &[u8]) {
		self.stream(series).extend_from_slice(values);
	}

	pub(super) fn terminated(&mut self, series: Series, values: &[u8]) {
		let stream = self.stream(series);
		stream.extend_from_slice(values);
		stream.push(0);
	}

	pub(super) fn tag(&mut self, id: i32) -> &mut Vec<u8> {
		self.tags.entry(id).or_default()
	}

	pub(super) fn into_blocks(self) -> impl Iterator<Item = (i32, Vec<u8>)> {
		Series::ALL
			.into_iter()
			.zip(self.series)
			.map(|(series, bytes)| (series.content_id(), bytes))
			.chain(self.tags)
			.filter(|(_, bytes)| !bytes.is_empty())
	}
}
