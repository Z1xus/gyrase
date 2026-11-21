# gyrase

Desktop GUI app for [gyrase](https://github.com/Z1xus/gyrase).

Turns whole-genome FASTQ, SAM, BAM, CRAM and VCF files into a VCF file and raw genotype files for 23andMe, AncestryDNA, FTDNA, MyHeritage and LivingDNA.

## Install

On Debian or Ubuntu, install the GUI libraries first:

```sh
sudo apt-get install -y libfontconfig-dev libxkbcommon-dev
```

Then install `gyrase`:

```sh
cargo install gyrase
```

## License

Licensed under [GPL-3.0-or-later](https://www.gnu.org/licenses/gpl-3.0.html).
