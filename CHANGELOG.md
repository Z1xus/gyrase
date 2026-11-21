# Changelog

## [0.1.0]

First release.

- Turn FASTQ, SAM, BAM or CRAM files into a VCF and raw data files for 23andMe, AncestryDNA, FTDNA, MyHeritage and LivingDNA.
- Turn a VCF or gVCF into raw data files.
- Align reads to GRCh38 with the built-in synapton aligner, then sort and mark duplicates.
- Sort and index unsorted or unindexed inputs automatically.
- Download and check the reference genome when it is missing.
- Use `gyr samtools` and `gyr bcftools` for common samtools and bcftools commands.
- Use the `gyrase` app to do all of this without a terminal.
