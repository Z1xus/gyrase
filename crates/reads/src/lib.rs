mod bam;
mod bytes;
pub mod cram;
mod error;
mod fastq;
mod header;
mod markdup;
mod record;
mod sam;
mod sort;
mod stats;

pub use bam::{
	BamReader, BamWriter, Located, Query, RecordBatch, Written, csi_geometry, index_bam,
	index_layout,
};
pub use error::{
	BamError, Error, FastqError, HeaderError, MarkError, RecordError, SamError, SamField, SortError,
};
pub use fastq::{
	FastqBatch, FastqReader, InterleavedFastq, PairedBatch, PairedFastq, ReadPair, ReadRef,
	interleaved, interleaved_pairs,
};
pub use header::{Field, Fields, Header, Line, Reference, SortOrder};
pub use markdup::{MarkSettings, MarkStats};
pub use record::{
	Array, ArrayWriter, Cigar, CigarKind, CigarOp, Element, Elements, Field as TagField, Flags,
	Qualities, QualityInput, Record, RecordRef, Sequence, TagIter, Tags, Value,
};
pub use sam::{SamReader, SamWriter};
pub use sort::{Sorted, Sorter};
pub use stats::{Flagstat, Idxstats, flagstat, idxstats};
