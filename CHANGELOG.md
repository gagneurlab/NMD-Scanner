# Changelog

## [0.3.0](https://github.com/gagneurlab/NMD-Scanner/compare/v0.2.0...v0.3.0) (2026-10-05)


### ⚠ BREAKING CHANGES

* read_vcf no longer returns the Qual, Filter and Info columns. Nothing in NMD-Scanner reads them, and polars-bio does not keep their text. A VCF now needs its header, at least the ##fileformat and #CHROM lines; without it read_vcf raises a ValueError that says so.
* take the stop codon from the annotation's stop_codon rows ([#31](https://github.com/gagneurlab/NMD-Scanner/issues/31))
* return the same columns and dtypes for every result, also an empty one ([#29](https://github.com/gagneurlab/NMD-Scanner/issues/29))

### Features

* add annotate() to get the result table without writing it ([#37](https://github.com/gagneurlab/NMD-Scanner/issues/37)) ([95d9f10](https://github.com/gagneurlab/NMD-Scanner/commit/95d9f10d43a81e9290a3e114364bf5877c6d50a4))
* read gene annotations from GFF3, not just GTF ([#35](https://github.com/gagneurlab/NMD-Scanner/issues/35)) ([d00904f](https://github.com/gagneurlab/NMD-Scanner/commit/d00904f426cbaf3228a49585da3362ecedeac845))
* read VCF and GFF3 with polars-bio, drop GTF input and pyranges ([#40](https://github.com/gagneurlab/NMD-Scanner/issues/40)) ([f405e57](https://github.com/gagneurlab/NMD-Scanner/commit/f405e5733eb0a7de55efe120aa5250034f860544))


### Bug Fixes

* compute exon numbers as ints under pandas 3 ([#30](https://github.com/gagneurlab/NMD-Scanner/issues/30)) ([d7b942c](https://github.com/gagneurlab/NMD-Scanner/commit/d7b942c5646ab417de7e0ae7d89e47976e1e675c))
* locate the CDS in the transcript by its coordinates ([#34](https://github.com/gagneurlab/NMD-Scanner/issues/34)) ([961340d](https://github.com/gagneurlab/NMD-Scanner/commit/961340db69977d7d9aab0c3e676283e582546ac8))
* measure the 50 nt rule and ptc_to_intron to the last exon junction ([#33](https://github.com/gagneurlab/NMD-Scanner/issues/33)) ([79fb2ff](https://github.com/gagneurlab/NMD-Scanner/commit/79fb2ff9482ec884fe68b2f93ceb10dc58eb6346))
* read VCF fields as text ([#36](https://github.com/gagneurlab/NMD-Scanner/issues/36)) ([053aa6a](https://github.com/gagneurlab/NMD-Scanner/commit/053aa6ac71e2b8599bba67e0862afaf564809e27))
* return the same columns and dtypes for every result, also an empty one ([#29](https://github.com/gagneurlab/NMD-Scanner/issues/29)) ([5e29ae2](https://github.com/gagneurlab/NMD-Scanner/commit/5e29ae25c0dae2e0f2f155872f24001b04f3d8f4))
* stop swapping UTR lengths on minus strand ([#32](https://github.com/gagneurlab/NMD-Scanner/issues/32)) ([5c6f5ae](https://github.com/gagneurlab/NMD-Scanner/commit/5c6f5ae646e3b21ae4fa767a1c12829e40240cdc))
* take the stop codon from the annotation's stop_codon rows ([#31](https://github.com/gagneurlab/NMD-Scanner/issues/31)) ([f33872d](https://github.com/gagneurlab/NMD-Scanner/commit/f33872d54350897ac91e6a6a16d56940483f5b42))
* write parquet output with a fixed schema and int exon numbers ([#28](https://github.com/gagneurlab/NMD-Scanner/issues/28)) ([6ec5c6e](https://github.com/gagneurlab/NMD-Scanner/commit/6ec5c6e58330a23b8c9a744bf5b4af6f97818b52))


### Documentation

* fix README and CONTRIBUTING inaccuracies ([#20](https://github.com/gagneurlab/NMD-Scanner/issues/20)) ([86559d9](https://github.com/gagneurlab/NMD-Scanner/commit/86559d97f6b3082f7274c8db9e822140178bdf3f))
* update the README for the 0.3.0 release ([#43](https://github.com/gagneurlab/NMD-Scanner/issues/43)) ([0098330](https://github.com/gagneurlab/NMD-Scanner/commit/0098330a025b4fdbd0f322c2c90795563eb7310a))

## [0.2.0](https://github.com/gagneurlab/NMD-Scanner/compare/v0.1.1...v0.2.0) (2026-05-14)


### Features

* cli single file output ([#17](https://github.com/gagneurlab/NMD-Scanner/issues/17)) ([8b6a1e5](https://github.com/gagneurlab/NMD-Scanner/commit/8b6a1e56ec9ebcf1cb428cb60a493e76e636ef7c))


### Bug Fixes

* nmd rule test and exports ([#7](https://github.com/gagneurlab/NMD-Scanner/issues/7)) ([4adc2ef](https://github.com/gagneurlab/NMD-Scanner/commit/4adc2ef9f6d18d7128a358af9d6f939f00b0a15e))
* reject multi-allelic VCF records in read_vcf ([#8](https://github.com/gagneurlab/NMD-Scanner/issues/8)) ([5811fb5](https://github.com/gagneurlab/NMD-Scanner/commit/5811fb5035a95bdc0eb65d738cff1c7ecad755c6))


### Documentation

* move misplaced module docstring into add_nmd_features ([#13](https://github.com/gagneurlab/NMD-Scanner/issues/13)) ([0adce83](https://github.com/gagneurlab/NMD-Scanner/commit/0adce8398a73e729c38b8f60f7f9a7e4184fef0d))

## [0.1.1](https://github.com/gagneurlab/NMD-Scanner/compare/v0.1.0...v0.1.1) (2026-05-13)


### Bug Fixes

* exon_number, ptc_exon_length bug, long and last exon rule, implementation of PTC_to_intron ([dd3edff](https://github.com/gagneurlab/NMD-Scanner/commit/dd3edff2c00c98b75162319aefb9b840ca10208a))
