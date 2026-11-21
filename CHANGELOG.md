# Changelog

## Unreleased
- Speed up `synapton` with a mappability certificate, unique-pair fast path, single-pass seeding and direct single-diagonal chaining
- Improve `synapton` accuracy by limiting ungapped acceptance to at most two interior mismatches and scaling down MAPQ gaps on near-tied alignments
- Make the `MD` tag optional in `synapton` and off by default (`AlignConfig.md`)
- Update the `synapton` index format to v3 (existing v2 index files rebuild automatically on the next run)
- Skip writing a BAM file by default for FASTQ input by keeping sorted alignments in memory and spilling to a temporary BAM above budget (use `gyr --bam` or `-b` to keep the BAM file)
- Encode BAM records in parallel across aligner worker threads
- Add an optional `gpu` Cargo feature (`wgpu` seeding front end, off by default) and `gyr --no-gpu`

## 0.1.0
- Initial release
