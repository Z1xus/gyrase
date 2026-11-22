# synapton

Short-read DNA aligner in safe Rust, built for [gyrase](https://github.com/Z1xus/gyrase).

## Usage

Open a reference FASTA with `gyrase_genome::Reference`, build or load an index with `Index::load_or_build`, create an `Aligner` with `AlignConfig::default()` and `gyrase_io::Workers::available()`, then call `align_reads` or `align_pairs`. Use `BamEmitter` to convert an `AlignedBatch` into `gyrase-reads` BAM records.

## License

Dual-licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at your option.
