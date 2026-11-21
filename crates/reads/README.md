# gyrase-reads

SAM, BAM, CRAM and FASTQ readers, writers, sorting, duplicate marking and read statistics.

Part of [gyrase](https://github.com/Z1xus/gyrase).

## Features

- Parallel BAM and CRAM readers and writers (`BamReader`, `BamWriter`, `cram::CramReader`, `cram::CramWriter`).
- SAM and FASTQ readers and writers (`SamReader`, `SamWriter`, `FastqReader`, `PairedFastq`).
- Coordinate sorting (`Sorter`), duplicate marking (`MarkSettings`, `MarkStats`) and BAM indexing (`index_bam`).
- Alignment statistics (`flagstat`, `idxstats`).
- Safe Rust (`#![deny(unsafe_code)]`).

## License

Dual-licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at your option.
