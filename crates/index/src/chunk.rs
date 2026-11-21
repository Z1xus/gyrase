use gyrase_io::VirtualOffset;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Chunk {
	pub begin: VirtualOffset,
	pub end: VirtualOffset,
}
