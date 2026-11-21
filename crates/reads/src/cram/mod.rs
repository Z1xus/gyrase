mod codec;
mod compression;
mod cursor;
mod decode;
mod emit;
mod encoding;
mod error;
mod reader;
mod records;
mod reference;
mod regenerate;
mod slice;
mod structure;
mod write;

pub use error::{Codec, Error};
pub use reader::{CramBatch, CramReader};
pub use write::{CramWriter, Finished};
