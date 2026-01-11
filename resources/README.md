# Resources

This directory contains datasets, reference files, and test resources used for the development, benchmarking, and validation of the NMD-Scanner.
The files in this directory are not required for normal usage of the NMD-Scanner package, but are provided for reproducibility, benchmarking, testing and example analyses.

---

## Benchmark datasets

### `MMRF_benchmark/`

- **`MMRF_TARGET_dataset.vcf`** : VCF file generated from the MMRF / TARGET dataset used to benchmark NMD efficiency prediction. The original dataset and processing scripts are provided by the original NMDEff project: <https://github.com/hjkng/nmdeff>

### `TCGA_benchmark/`

- **`tcga_dataset.vcf`** : VCF file generated from the TCGA dataset as described in the original study (<https://github.com/hjkng/nmdeff>)
  used for NMD efficiency benchmarking and for model training.

---

## Files for Testing

### `test_files/`

VCF files created for testing the NMD-Scanner functionality.

- **`test_variants.vcf`** : Test VCF containing variants on the plus strand only.

- **`test_variants_minus.vcf`** : Test VCF containing variants on the minus strand only.

- **`variants.vcf`** : Combined test VCF containing both plus- and minus-strand variants.

All test VCFs were generated using:

```text
/scripts/create_test_VCF.ipynb
```
and include edge-case variants such as start and stop-loss mutations, frameshifts and nonsense mutations.

### `test_output_files/`

Contains intermediate and final output files used for testing of the computation of NMD features (see `../tests/test_rules.py`).

---

## Example input files

These files are provided to allow testing on a small genomic region (chromosome 18).

- **`chr18.fa.gz`** : Reference FASTA for chromosome 18.

- **`chr18.fa.gz.fai`** : FASTA index file.

- **`chr18.fa.gz.gzi`** : FASTA gzip index.

- **`chr18.gtf.gz`** : GTF annotation file for chromosome 18.

- **`part-00241-61a0abbf-fbf9-444f-8287-4e46ad4b9b7b-c000.vcf`** : Example VCF file containing some variants from chromosome 18, used as a realistic input example for testing and demonstration.

--- 

## Notes
- These resources are intended for testing, benchmarking and reproducibility.
- Users of the NMD-Scanner do not need to download or use these files to run the tool on their own data.
- Dataset usage of the TCGA and MMRF/TARGET datasets follows the terms of the original source (<https://github.com/hjkng/nmdeff>).