mod bgzf;
mod error;
mod format;
mod workers;

pub use bgzf::{
	BgzfReader, BgzfWriter, Block, BlockPosition, BlockReader, BlockSpan, Finished, Input, Level,
	MAX_BLOCK_LEN, VirtualOffset,
};
pub use error::Error;
pub use format::{Compression, Decoded, decode};
pub use workers::Workers;
