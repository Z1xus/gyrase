mod bgzf;
mod error;
mod format;
mod workers;

pub use bgzf::{
	BgzfReader, BgzfWriter, Block, BlockPosition, BlockReader, BlockSpan, Chunk, Finished, Input,
	Level, MAX_BLOCK_LEN, Readahead, VirtualOffset,
};
pub use error::Error;
pub use format::{Compression, Decoded, decode, peek};
pub use workers::Workers;
