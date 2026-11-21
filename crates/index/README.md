# gyrase-index

Genomics index readers, writers and query planners for BAI, CSI, TBI and CRAI files.

Part of [gyrase](https://github.com/Z1xus/gyrase).

## Features

- Read, write and build BAI, CSI and TBI indexes (`Index`, `Builder`, `Geometry`, `TabixHeader`).
- Read and write CRAM `.crai` indexes (`Crai`).
- Plan indexed region queries into coalesced BGZF block spans (`Plan`).
- Safe Rust (`#![deny(unsafe_code)]`).

## License

Dual-licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at your option.
