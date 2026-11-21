use crate::alignment::{CigarKind, CigarOp};
use crate::nucleotide::ascii;
use crate::num::u32_of;

pub(crate) struct Alignment<'a> {
	pub(crate) operations: &'a [CigarOp],
	pub(crate) query: &'a [u8],
	pub(crate) reference: &'a [u8],
	pub(crate) left_clip: u32,
	pub(crate) right_clip: u32,
}

pub(crate) fn describe(
	alignment: &Alignment<'_>,
	cigar: &mut Vec<CigarOp>,
	md: &mut Vec<u8>,
	with_md: bool,
) -> u32 {
	cigar.clear();
	md.clear();
	if alignment.left_clip > 0 {
		cigar.push(CigarOp::new(CigarKind::SoftClip, alignment.left_clip));
	}
	cigar.extend_from_slice(alignment.operations);
	if alignment.right_clip > 0 {
		cigar.push(CigarOp::new(CigarKind::SoftClip, alignment.right_clip));
	}
	let (mut q, mut r) = (0usize, 0usize);
	let mut run = 0u32;
	let mut edits = 0u32;
	for operation in alignment.operations {
		let length = operation.length() as usize;
		match operation.kind() {
			CigarKind::Match => {
				let pairs = alignment
					.query
					.iter()
					.skip(q)
					.zip(alignment.reference.iter().skip(r))
					.take(length);
				for (&query, &reference) in pairs {
					if query == reference && reference < 4 {
						run += 1;
					} else {
						if with_md {
							push_number(md, run);
							md.push(ascii(reference));
						}
						run = 0;
						edits += 1;
					}
				}
				q += length;
				r += length;
			}
			CigarKind::Insertion => {
				edits += u32_of(length);
				q += length;
			}
			CigarKind::Deletion => {
				if with_md {
					push_number(md, run);
					md.push(b'^');
					md.extend(
						alignment
							.reference
							.iter()
							.skip(r)
							.take(length)
							.map(|&code| ascii(code)),
					);
				}
				run = 0;
				edits += u32_of(length);
				r += length;
			}
			CigarKind::SoftClip => q += length,
		}
	}
	if with_md {
		push_number(md, run);
	}
	edits
}

fn push_number(out: &mut Vec<u8>, value: u32) {
	let mut digits = [0u8; 10];
	let mut remaining = value;
	let mut count = 0;
	loop {
		if let Some(slot) = digits.get_mut(count) {
			*slot = b'0' + u8::try_from(remaining % 10).unwrap_or(0);
		}
		count += 1;
		remaining /= 10;
		if remaining == 0 {
			break;
		}
	}
	out.extend(digits.iter().take(count).rev());
}
