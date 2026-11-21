# synapton

Short-read DNA aligner in safe Rust, built for [gyrase](https://github.com/Z1xus/gyrase).

On chr20 (2.5M read pairs, 12 threads, Ryzen 5 7600X), synapton is about 5 times faster than bwa-mem2 (14.91 s vs 76.33 s), 3 times faster than minimap2 (45.06 s) and slightly faster than strobealign (15.92 s), at 316 MiB peak RSS and 0.029% misplaced reads at MAPQ 20 or more.

## Usage

Open a reference FASTA with `gyrase_genome::Reference`, build or load an index with `Index::load_or_build`, create an `Aligner` with `AlignConfig::default()` and `gyrase_io::Workers::available()`, then call `align_reads` or `align_pairs`. Use `BamEmitter` to convert an `AlignedBatch` into `gyrase-reads` BAM records.

## License

Dual-licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at your option.
