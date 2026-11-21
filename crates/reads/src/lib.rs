mod bam;
mod bytes;
pub mod cram;
mod error;
mod fastq;
mod header;
mod record;
mod sam;
mod stats;

pub use bam::{BamReader, BamWriter, Chunk, Located, Query, RecordBatch, Written};
pub use error::{
	BamError, Error, FastqError, HeaderError, RecordError, SamError, SamField, StatsError,
};
pub use fastq::{
	FastqBatch, FastqReader, InterleavedFastq, PairedBatch, PairedFastq, ReadPair, ReadRef,
	interleaved_pairs,
};
pub use header::{Field, Fields, Header, Line, Reference};
pub use record::{
	Array, ArrayWriter, Cigar, CigarKind, CigarOp, Element, Elements, Field as TagField, Flags,
	Qualities, QualityInput, Record, RecordRef, Sequence, TagIter, Tags, Value,
};
pub use sam::{SamReader, SamWriter};
pub use stats::{Flagstat, Idxstats, flagstat, idxstats};
