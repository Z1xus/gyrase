mod block;
mod level;
mod offset;
mod reader;
mod stream;
mod writer;

pub(crate) use block::is_bgzf;
pub use block::{Block, MAX_BLOCK_LEN};
pub use level::Level;
pub use offset::{BlockPosition, VirtualOffset};
pub use reader::{BlockReader, BlockSpan, Input};
pub use stream::BgzfReader;
pub use writer::{BgzfWriter, Finished};
