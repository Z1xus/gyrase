# gyrase-cli

Command-line interface (`gyr`) for [gyrase](https://github.com/Z1xus/gyrase).

Turns FASTQ, SAM, BAM, CRAM and VCF files into a VCF file and raw genotype files for 23andMe, AncestryDNA, FTDNA, MyHeritage and LivingDNA. Also provides `gyr samtools` and `gyr bcftools` subcommands.

## Install

```sh
cargo install gyrase-cli
```

## Usage

```sh
gyr sample.bam
gyr sample_R1.fastq.gz sample_R2.fastq.gz -o out
gyr sample_R1.fastq.gz sample_R2.fastq.gz --bam
gyr sample.cram -f 23andme,ancestry
gyr samtools flagstat sample.bam
gyr bcftools view sample.vcf.gz
```

## License

Licensed under [GPL-3.0-or-later](https://www.gnu.org/licenses/gpl-3.0.html).
