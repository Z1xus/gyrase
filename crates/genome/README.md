# gyrase-genome

Core genomics types and indexed FASTA reference access.

Part of [gyrase](https://github.com/Z1xus/gyrase).

## Features

- Nucleotide base (`Base`) and SAM 4-bit code (`Code4`) types.
- Contig tables (`Contigs`, `Contig`, `ContigId`) with chromosome alias resolution (`chr1`, `1`, `NC_000001.11`, `chrM`, `MT`).
- Half-open genomic intervals (`Interval`) and region string parsing.
- Human reference build detection (`Build`: GRCh37, GRCh38, CHM13).
- Random-access FASTA reader (`Reference`) for plain and BGZF-compressed FASTA files, with automatic `.fai` indexing.
- Safe Rust (`#![deny(unsafe_code)]`).

## License

Dual-licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at your option.
