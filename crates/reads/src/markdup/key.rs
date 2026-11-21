use crate::record::{Cigar, CigarKind, Flags, RecordRef};

const FF: u8 = 2;
const RR: u8 = 3;
const FR: u8 = 5;
const RF: u8 = 7;
const LEFT: u8 = 11;
const RIGHT: u8 = 13;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) struct SingleKey {
	contig: i32,
	coord: i64,
	orientation: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) struct PairKey {
	contig: i32,
	coord: i64,
	mate_contig: i32,
	mate_coord: i64,
	leftmost: u8,
	orientation: u8,
}

#[derive(Debug, Clone, Copy, Default)]
pub(super) struct Extent {
	pub start: i64,
	pub end: i64,
}

impl SingleKey {
	pub(super) fn coord(&self) -> i64 {
		self.coord
	}
}

pub(super) fn extent(record: RecordRef<'_>) -> Extent {
	let position = record.position().map_or(-1, i64::from);
	let cigar = record.cigar();
	let leading: i64 = cigar
		.ops()
		.take_while(|op| is_clip(op.kind()))
		.map(|op| i64::from(op.len()))
		.sum();
	let trailing = trailing_clips(cigar);
	let end = record
		.reference_end()
		.and_then(|end| i64::try_from(end).ok())
		.unwrap_or(0);
	Extent {
		start: position - leading + 1,
		end: end + trailing,
	}
}

pub(super) fn mate_extent(position: i64, cigar: &[u8]) -> Extent {
	let leading: i64 = text_ops(cigar)
		.take_while(|&(_, symbol)| matches!(symbol, b'S' | b'H'))
		.map(|(count, _)| count)
		.sum();
	let mut aligned = false;
	let mut span = 0;
	for (count, symbol) in text_ops(cigar) {
		match symbol {
			b'M' | b'D' | b'N' | b'=' | b'X' => {
				span += count;
				aligned = true;
			}
			b'S' | b'H' if aligned => span += count,
			_ => {}
		}
	}
	Extent {
		start: position - leading + 1,
		end: position + span,
	}
}

pub(super) fn single_key(contig: i32, flags: Flags, extent: Extent) -> SingleKey {
	if flags.is_reverse() {
		SingleKey {
			contig,
			coord: extent.end,
			orientation: RR,
		}
	} else {
		SingleKey {
			contig,
			coord: extent.start,
			orientation: FF,
		}
	}
}

pub(super) fn pair_key(
	contig: i32,
	mate_contig: i32,
	flags: Flags,
	this: Extent,
	mate: Extent,
) -> PairKey {
	let reverse = flags.is_reverse();
	let mate_reverse = flags.is_mate_reverse();
	let first = flags.is_first_in_pair();
	let leftmost = if contig != mate_contig {
		contig < mate_contig
	} else if reverse == mate_reverse {
		if reverse {
			this.end <= mate.end
		} else {
			this.start <= mate.start
		}
	} else if reverse {
		this.end <= mate.start
	} else {
		this.start <= mate.end
	};
	let (coord, mate_coord, orientation) = match (leftmost, reverse == mate_reverse) {
		(true, true) => (this.start, mate.end, if reverse == first { RR } else { FF }),
		(true, false) if !reverse => (this.start, mate.end, FR),
		(true, false) => (this.end, mate.start, RF),
		(false, true) => (this.end, mate.start, if reverse == first { FF } else { RR }),
		(false, false) if !reverse => (this.start, mate.end, RF),
		(false, false) => (this.end, mate.start, FR),
	};
	PairKey {
		contig,
		coord,
		mate_contig,
		mate_coord,
		leftmost: if leftmost { LEFT } else { RIGHT },
		orientation,
	}
}

fn is_clip(kind: CigarKind) -> bool {
	matches!(kind, CigarKind::SoftClip | CigarKind::HardClip)
}

fn trailing_clips(cigar: Cigar<'_>) -> i64 {
	cigar.ops().fold(0, |clipped, op| {
		if is_clip(op.kind()) {
			clipped + i64::from(op.len())
		} else {
			0
		}
	})
}

fn text_ops(cigar: &[u8]) -> impl Iterator<Item = (i64, u8)> + '_ {
	let mut rest = cigar;
	std::iter::from_fn(move || {
		if rest.first().is_none_or(|&byte| byte == b'*') {
			return None;
		}
		let digits = rest.iter().take_while(|byte| byte.is_ascii_digit()).count();
		let count = if digits == 0 {
			1
		} else {
			rest[..digits].iter().fold(0i64, |number, &digit| {
				number
					.saturating_mul(10)
					.saturating_add(i64::from(digit - b'0'))
			})
		};
		let symbol = *rest.get(digits)?;
		rest = &rest[digits + 1..];
		Some((count, symbol))
	})
}
