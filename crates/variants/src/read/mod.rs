mod any;
mod bcf;
mod chunks;
mod header;
mod regions;
mod vcf;

pub use any::Reader;
pub use bcf::{BcfReader, RecordBatch};
pub use chunks::ChunkReader;
pub use header::read_header_text;
pub use regions::open_regions;
pub use vcf::{LineBatch, VcfReader};
