# Scripts

These scripts are not part of the package. They run from the repo root with `uv run --group scripts python scripts/<name>.py`, and they use the `nmd_scanner` of the checkout.

## make_benchmark_vcfs.py

Builds the benchmark VCFs `resources/TCGA_benchmark/tcga_dataset.vcf` and `resources/MMRF_benchmark/MMRF_TARGET_dataset.vcf` from the NMDEff study tables. They have REF and ALT on the forward strand of GRCh38, and they equal the committed files byte for byte.

```bash
uv run --group scripts python scripts/make_benchmark_vcfs.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa
```

- Inputs: a GENCODE GFF3 of GRCh38 and a GRCh38 FASTA with chr names and a `.fai` index.
- Outputs: the two VCFs above, one row per table row, in table order.
- Data: `tcga_dataset.csv` and `MMRF_TARGET_dataset.csv` of <https://github.com/hjkng/nmdeff>, downloaded at a pinned commit and checked against their sha256.
- It parses the substitution from HGVSc, which has the alleles on the transcript strand. It takes the strand from the GFF3 and complements REF and ALT on the minus strand. For transcripts missing from the GFF3, it uses the orientation whose REF matches the FASTA.
- It fails if a row is not a single-base substitution, or if a REF does not match the FASTA.
- `train_model.py` and `validate_model.py` import these steps from it.

## train_model.py

Trains the NMD efficiency model that the package ships as `src/nmd_scanner/nmd_efficiency_rf.onnx`. NMD-Scanner adds its prediction to the output by default (see the main README).

```bash
uv run --group scripts python scripts/train_model.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa --out-dir out/
```

- Inputs: the GFF3 and the FASTA, as above.
- Data: `tcga_dataset.csv` of NMDEff, at the pinned commit and with the sha256 check of `make_benchmark_vcfs.py`. Its column `NMD_efficiency` is the target. The variants are those of `resources/TCGA_benchmark/tcga_dataset.vcf`.
- It tunes a random forest with a grid search. The folds are grouped by chromosome.
- In nested cross-validation with the same folds, it predicts each row out of fold.
- It saves the random forest tuned on all rows as ONNX. It fails if the ONNX model predicts the training rows differently from the random forest.
- Outputs in `--out-dir`:
  - `nmd_efficiency_rf.onnx`, the model.
  - `oof_predictions.parquet`, the out-of-fold predictions, which `validate_model.py` reads.
  - `cv_results.json`, the R2 of each fold and the hyperparameters.
  - `variants.vcf`, the input of `nmd_scanner.annotate`.

## validate_model.py

Computes the numbers of the section "Model validation" in the main README. It tests `src/nmd_scanner/nmd_efficiency_rf.onnx` without retraining on TCGA out of fold, Geuvadis, GTEx and MMRF-TARGET.

```bash
uv run --group scripts python scripts/validate_model.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa \
    --train-dir train/ --out-dir out/
```

- Inputs: the GFF3 and the FASTA, as above. `--train-dir` is the `--out-dir` of the `train_model.py` run that trained the model. The TCGA row uses its out-of-fold predictions. The script fails if the model of that run differs from `nmd_efficiency_rf.onnx`.
- Data, with a sha256 check for each file:
  - `MMRF_TARGET_dataset.csv` of NMDEff, at the commit that `make_benchmark_vcfs.py` pins.
  - Supplementary Data 2 of Kim et al. 2024 (GTEx).
  - 2 tables of Zenodo record 16666299 (Geuvadis). The record is one zip of 3.2 GB, so the script reads only those 2 members with HTTP range requests. `--geuvadis-dir` reads them from a local directory instead, with the same sha256 check.
- It builds one VCF of the variants and runs `nmd_scanner.annotate`. For GTEx and MMRF-TARGET, it converts HGVSc to genomic alleles with `make_benchmark_vcfs.py`.
- It keeps the rows the model can predict that have a variant allele fraction above 0 in the RNA. For the germline variants, the NMD efficiency is -log2(VAF_RNA / 0.5).
- The 95% CIs come from 2,000 bootstrap resamples of the rows per cohort.
- Output: `out/validation.json`, and a short summary in the log. Per cohort, the file has the Spearman correlations of the model and of the rules, the gain of the model over the rules, and the mean measured NMD efficiency of the rows with `nmd_escape` False. For GTEx, it also has that mean over the variant-tissue rows, in total and per tissue. The README table is written by hand from these numbers.
