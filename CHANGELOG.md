# Changelog

## [0.4.0](https://github.com/gagneurlab/NMD-Scanner/compare/v0.3.0...v0.4.0) (2026-10-07)


### ⚠ BREAKING CHANGES

* raise a KeyError in add_nmd_features and evaluate_nmd_escape_rules for a row without a column they read
* rename the model inputs ptc_to_intron and stop_codon_distance
* rename the output columns to one naming scheme
* drop the output columns that carry no information of their own
* make the NMD rules null on rows without a PTC or with a null rule input
* type the list columns as named structs and keep their shape through Parquet and CSV
* make the closed value sets categorical and give variant_id null without an ID
* start the output with the columns that identify a row
* remove cli.parquet_schema and cli.to_parquet_safe

### Features

* add nmd_model_status and MODEL_INPUTS for the NMD efficiency model ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* add ptc_pos_in_alt_transcript, ptc_exon_number and stop_classification ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* add sequences=False and --no-sequences to drop the 4 sequence columns ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* add the output columns unknown_reason, alt_cds_start_in_transcript, cds_frame, has_start_codon and alt_transcript_exon_info ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* add to_arrow to convert the result table to a typed Arrow table ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* drop the output columns that carry no information of their own ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* make the closed value sets categorical and give variant_id null without an ID ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* remove cli.parquet_schema and cli.to_parquet_safe ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* rename the model inputs ptc_to_intron and stop_codon_distance ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* rename the output columns to one naming scheme ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* start the output with the columns that identify a row ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* type the list columns as named structs and keep their shape through Parquet and CSV ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))


### Bug Fixes

* apply 5'UTR changes at the CDS start to the alt transcript ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* classify the first in-frame stop codon of the alt transcript ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* classify the rescued ORF after a start loss ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* compute the exon features of a PTC from the exons of the alt transcript ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* give a row with null transcript columns if no joined transcript has exon rows ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* give each VCF record its own row ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* give the length of the PTC exon in the alt transcript as ptc_exon_length ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* handle variants at exon boundaries by their splice sites ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* judge start loss on the annotated start codon ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* make the NMD rules null on rows without a PTC or with a null rule input ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* measure ptc_to_intron of a last-exon PTC to the transcript end ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* measure stop_codon_distance to the stop codon in the alt transcript ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* measure the start-proximal rule from the rescued ATG after a start loss ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* number the exons of the alt transcript scan by the alt exon lengths ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* raise a ValueError for a CDS row outside the exon rows of its transcript ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* raise a ValueError that names the broken input rule ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* read codons in the frame of the annotated CDS ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* require a leading ATG to flag start loss ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* skip variants with a symbolic ALT allele or a breakend ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* skip VCF records with ALT "." or "*" ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* start the transcript codon scan at the CDS start in the transcript ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* use the annotated start codon as the start codon position ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))


### Performance Improvements

* build the alt CDS of each variant without pandas copies ([ae86a07](https://github.com/gagneurlab/NMD-Scanner/commit/ae86a0728871e1a9f38f5d70d8b7b2c9f1292788))
* collect the codon scan results in lists instead of DataFrame cells ([ae86a07](https://github.com/gagneurlab/NMD-Scanner/commit/ae86a0728871e1a9f38f5d70d8b7b2c9f1292788))
* compute the exon boundaries once per transcript ([ae86a07](https://github.com/gagneurlab/NMD-Scanner/commit/ae86a0728871e1a9f38f5d70d8b7b2c9f1292788))
* scan the whole alt transcript for stop codons only on a stop loss ([ae86a07](https://github.com/gagneurlab/NMD-Scanner/commit/ae86a0728871e1a9f38f5d70d8b7b2c9f1292788))


### Documentation

* define the transcript position and the exon number, and drop the old pipeline description ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* describe short alt transcripts, chromosome ends and REF ends at the coding region edge ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* describe the meaning, dtype and null cases of every output column ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* draw the transcript layout of the ptc_to_intron and NMD rule cases ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* drop the GFF3 details from the new pipeline description ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* drop the history from the docs and the test comments ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* format the markdown files with prettier ([5526ad0](https://github.com/gagneurlab/NMD-Scanner/commit/5526ad0f001a12385c77fd7519926f511135df90))
* give ptc_to_start_codon a null clause for a stop codon as start codon ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* move the GFF3 defects and input checks into Input Defects.md ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* say that the output for CDS rows that share bases can be wrong ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))
* say why likely_misannotated does not flag a stop codon as start codon ([92c2ac4](https://github.com/gagneurlab/NMD-Scanner/commit/92c2ac41a99ef1da6ebaced23a74262b1b26c52e))


### Code Refactoring

* raise a KeyError in add_nmd_features and evaluate_nmd_escape_rules for a row without a column they read ([ae86a07](https://github.com/gagneurlab/NMD-Scanner/commit/ae86a0728871e1a9f38f5d70d8b7b2c9f1292788))

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
