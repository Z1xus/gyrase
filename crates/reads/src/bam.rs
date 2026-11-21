mod batch;
mod frame;
mod index;
mod indexed;
mod query;
mod reader;
mod writer;

pub use batch::{Located, RecordBatch};
pub(crate) use index::fits_bai;
pub use index::{csi_geometry, index_bam};
pub(crate) use indexed::IndexedWriter;
pub use query::Query;
pub use reader::BamReader;
pub use writer::{BamWriter, Written};
