# synapton

Fast short-read DNA aligner written in safe Rust, with no `unsafe` code.

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

Use `BamEmitter` to convert an `AlignedBatch` into `gyrase-reads` BAM records with standard `NM`, `AS`, `XS`, `SA`, `MC` and `ms` tags. Set `AlignConfig.md` to `true` to also emit the `MD` tag (off by default). Enable the optional `gpu` Cargo feature for `wgpu`-accelerated seeding on indexes that fit in one storage buffer binding (larger indexes fall back to CPU automatically).

## Benchmark (chr20)

Conditions: 2,500,000 paired-end reads (5,000,000 reads, 151 bp) from HG002 NovaSeq PCR-free 35x on GRCh38 chromosome 20 only, prebuilt index, full SAM output to disk, 12 threads, AMD Ryzen 5 7600X (6 cores, 12 threads) and NVIDIA GeForce RTX 4060 (for the optional `gpu` feature build). Median of three runs per tool. Agreement is the share of reads placed within 10 bp of the original HG002 bwa-mem BAM inside GIAB v4.2.1 high-confidence regions where the original MAPQ is at least 20. Misplaced at MAPQ >= 20 is the share of reads with original MAPQ at least 20 placed more than 50 bp away or on the opposite strand with MAPQ at least 20.

| Aligner | Wall time (chr20) | Speed (chr20) | Peak RSS (chr20) | Agreement in GIAB BED (chr20) | Misplaced at MAPQ >= 20 |
| --- | --- | --- | --- | --- | --- |
| synapton (CPU) | 14.91 s | 335,345 reads/s | 316 MiB | 99.950% | 0.029% |
| synapton (GPU, `gpu` feature) | 12.06 s | 414,594 reads/s | 710 MiB | 99.950% | 0.029% |
| strobealign 0.20.0 | 15.92 s | 314,070 reads/s | 530 MiB | 99.858% | 0.158% |
| minimap2 2.31 | 45.06 s | 110,963 reads/s | 928 MiB | 99.499% | 0.077% |
| bwa-mem2 2.2.1 | 76.33 s | 65,505 reads/s | 3,766 MiB | 99.995% | 0.023% |

These numbers are for chromosome 20 only. On this chromosome 20 test, `synapton` (CPU) is about 1.07 times faster than `strobealign`, 3.0 times faster than `minimap2` and 5.1 times faster than `bwa-mem2`, and uses less memory than all three. `bwa-mem2` has higher placement agreement with the original bwa-mem BAM (99.995% vs 99.950% in GIAB high-confidence regions). `synapton` has fewer misplaced reads at MAPQ 20 or more than `minimap2` and `strobealign` (0.029% vs 0.077% and 0.158%) and slightly more than `bwa-mem2` (0.023%), while assigning MAPQ 20 or more to 92.2% of mapped reads (vs 94.0% for `bwa-mem2`).

## License

Dual-licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at your option.
