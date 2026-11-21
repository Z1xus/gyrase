# gyrase-variants

VCF and BCF readers, writers, indexing, normalization, merging, concatenation, querying and filtering.

Part of [gyrase](https://github.com/Z1xus/gyrase).

## Features

- Read and write VCF, bgzipped VCF and BCF files (`Reader`, `VcfReader`, `BcfReader`, `Writer`).
- Indexed region queries (`open_regions`) with tabix and CSI index collection (`IndexCollector`).
- Variant normalization (`Norm`), merging (`Merge`), concatenation (`Concat`), formatting queries (`Query`) and filtering (`Filter`, `View`).
- Safe Rust (`#![deny(unsafe_code)]`).

## License

Licensed under [GPL-3.0-or-later](https://www.gnu.org/licenses/gpl-3.0.html).
