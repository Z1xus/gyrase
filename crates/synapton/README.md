# synapton

Fast short-read DNA aligner written in safe Rust (`#![deny(unsafe_code)]`).

`synapton` is the aligner used inside [gyrase](https://github.com/Z1xus/gyrase). It indexes a reference genome with open-syncmer randstrobes, chains seed hits, extends alignments with SIMD-dispatched banded affine-gap dynamic programming (`multiversion`), pairs mates and computes MAPQ scores. Output is deterministic for any thread count.

## Usage

Add `synapton`, `gyrase-genome` and `gyrase-io` to `Cargo.toml`:

```rust
use std::path::Path;
use std::sync::atomic::AtomicBool;

use gyrase_genome::Reference;
use gyrase_io::Workers;
use synapton::{AlignConfig, Aligner, Index, Placement, Read, ReadPair, SeedConfig};

fn main() -> Result<(), synapton::Error> {
    let workers = Workers::available()?;
    let reference = Reference::open(Path::new("reference.fa"))?;
    let stop = AtomicBool::new(false);
    let index = Index::load_or_build(
        Path::new("reference.syni"),
        &reference,
        SeedConfig::short_read(),
        &workers,
        &stop,
    )?;

    let aligner = Aligner::new(&index, AlignConfig::default(), &workers);
    let pairs = [ReadPair {
        first: Read {
            name: b"read1",
            sequence: b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            qualities: b"",
        },
        second: Read {
            name: b"read1",
            sequence: b"TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
            qualities: b"",
        },
    }];

    let batch = aligner.align_pairs(&pairs)?;
    for record in batch.records() {
        if let Placement::Mapped(alignment) = &record.placement {
            let contig = index.contigs().get(alignment.contig);
            println!(
                "{:?} pos={} mapq={}",
                contig.map(|c| c.name.as_str()),
                alignment.position,
                alignment.mapq
            );
        }
    }
    Ok(())
}
```

Use `BamEmitter` to convert an `AlignedBatch` into `gyrase-reads` BAM records with standard `NM`, `MD`, `AS`, `XS`, `SA`, `MC` and `ms` tags.

## Benchmark (chr20)

Conditions: 2,500,000 paired-end reads (5,000,000 reads, 151 bp) from HG002 NovaSeq PCR-free 35x on GRCh38 chromosome 20 only, prebuilt index, full SAM output to disk, 12 threads, AMD Ryzen 5 7600X (6 cores, 12 threads). Faster of two runs per tool. Agreement is the share of reads placed within 10 bp of the original HG002 bwa-mem BAM inside GIAB v4.2.1 high-confidence regions where the original MAPQ is at least 20.

| Aligner | Wall time (chr20) | Speed (chr20) | Peak RSS (chr20) | Agreement in GIAB BED (chr20) |
| --- | --- | --- | --- | --- |
| synapton | 15.97 s | 313,087 reads/s | 315 MiB | 99.950% |
| strobealign 0.20.0 | 17.29 s | 289,184 reads/s | 595 MiB | 99.858% |
| minimap2 2.31 | 53.78 s | 92,971 reads/s | 962 MiB | 99.499% |
| bwa-mem2 2.2.1 | 85.67 s | 58,363 reads/s | 3,788 MiB | 99.995% |

These numbers are for chromosome 20 only. `bwa-mem2` has higher placement agreement with the original bwa-mem BAM (99.995% vs 99.950% in GIAB high-confidence regions), and all three baseline tools have fewer discordant high-MAPQ reads outside GIAB high-confidence regions (0.180% for `synapton` vs 0.158% for `strobealign`, 0.077% for `minimap2` and 0.023% for `bwa-mem2`).

## License

Dual-licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at your option.
