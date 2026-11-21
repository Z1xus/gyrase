mod decode;
mod encode;
mod record;
pub(crate) mod typed;

pub use encode::Encoder;
pub use record::{BcfRecord, FormatEntry, InfoEntry};
pub use typed::Typed;
