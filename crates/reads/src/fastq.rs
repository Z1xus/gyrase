mod batch;
mod pairs;
mod reader;

pub use batch::{FastqBatch, ReadRef};
pub use pairs::{InterleavedFastq, PairedBatch, PairedFastq, ReadPair, interleaved_pairs};
pub use reader::FastqReader;
