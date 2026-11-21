mod bin;
mod cigar;
mod fixed;
mod flags;
mod owned;
mod sequence;
mod tags;
mod view;

pub use cigar::{Cigar, CigarKind, CigarOp};
pub use flags::Flags;
pub use owned::{ArrayWriter, QualityInput, Record};
pub use sequence::{Qualities, Sequence};
pub use tags::{Array, Element, Elements, Field, TagIter, Tags, Value};
pub use view::RecordRef;
