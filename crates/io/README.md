# gyrase-io

Parallel BGZF, gzip and stream I/O primitives for genomics tools.

Part of [gyrase](https://github.com/Z1xus/gyrase).

## Features

- Parallel BGZF reader (`BgzfReader`) and writer (`BgzfWriter`) using `libdeflater`.
- Compression sniffing and transparent decoding (`decode`, `peek`) for plain, gzip and BGZF streams.
- Virtual offset (`VirtualOffset`), chunk (`Chunk`) and block span (`BlockSpan`) types for indexed seeking.
- Explicit worker pool (`Workers`) with no global state.
- Safe Rust (`#![deny(unsafe_code)]`).

## License

Dual-licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at your option.
