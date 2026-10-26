# Scripts

These scripts are not part of the package. They run from the repo root with `uv run --group scripts python scripts/<name>.py`, and they use the `nmd_scanner` of the checkout.

`train_model.py` and `validate_model.py` are also notebooks. They are jupytext percent scripts, and each is paired with a `.ipynb` of the same name that git ignores. As a notebook, a script takes its arguments from an environment variable, `TRAIN_MODEL_ARGS` or `VALIDATE_MODEL_ARGS`, in the syntax of the command line. The kernel runs in `scripts/`, so the paths in the variable must be absolute. The notebook commands below write the `.ipynb` with jupytext, execute it with nbconvert, and render it to HTML with the figures inline. `TQDM_DISABLE=1` keeps progress bars out of the HTML.

Both ways draw each figure with plotnine and save it at 600 dpi as `figures/<name>.png` and `figures/<name>.pdf` in `--out-dir`. `figures.py` has their theme.

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
- `join_features` joins the rows of a study table to the `nmd_scanner.annotate` output, keeps the rows with `nmd_efficiency_prediction_status` ok, and logs the rows that it drops.
- `read_genes` returns the gene of each transcript, for the gene bootstrap of `validate_model.py`.

## train_model.py

Trains the NMD efficiency model that the package ships as `src/nmd_scanner/nmd_efficiency_rf.onnx`. NMD-Scanner adds its prediction to the output by default (see the main README).

As a script:

```bash
uv run --group scripts python scripts/train_model.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa --out-dir train/
```

As a notebook, rendered to `train/train_model.html`:

```bash
uv run --group scripts jupytext --to ipynb scripts/train_model.py &&
    TQDM_DISABLE=1 TRAIN_MODEL_ARGS="--gff3 $PWD/gencode.v42.annotation.gff3.gz --fasta $PWD/GRCh38.fa --out-dir $PWD/train" \
    uv run --group scripts jupyter nbconvert --to html --execute scripts/train_model.ipynb --output-dir train/
```

- Inputs: the GFF3 and the FASTA, as above.
- Data: `tcga_dataset.csv` of NMDEff, at the pinned commit and with the sha256 check of `make_benchmark_vcfs.py`. Its column `NMD_efficiency` is the target. The variants are those of `resources/TCGA_benchmark/tcga_dataset.vcf`.
- It logs the rows that it drops, and why.
- It tunes a random forest with a grid search. The folds are grouped by chromosome.
- In nested cross-validation with the same folds, it predicts each row out of fold.
- It saves the random forest tuned on all rows as ONNX. It fails if the ONNX model predicts the training rows differently from the random forest.
- Outputs in `--out-dir`:
  - `nmd_efficiency_rf.onnx`, the model.
  - `nmd_efficiency_rf.pkl`, the same random forest as a joblib pickle. It loads with the scikit-learn version of `uv.lock`. The repo root has a copy. The wheel and the sdist do not contain it.
  - `oof_predictions.parquet`, the out-of-fold predictions, which `validate_model.py` reads.
  - `cv_results.json`, the R2 of each fold and the hyperparameters.
  - `variants.vcf`, the input of `nmd_scanner.annotate`.
- Figures:
  - `target_histogram`: the histogram of `NMD_efficiency`.
  - `tcga_measured_vs_predicted`: the measured NMD efficiency over the out-of-fold prediction, with Spearman and R2.
  - `calibration_oof`: per decile of the out-of-fold prediction, the mean measured over the mean predicted NMD efficiency, with the 95% CI of the measured mean. It adds the least-squares line of measured on predicted over all rows, with its slope and intercept.
  - `feature_importance`: the impurity importance and the permutation importance of the 19 inputs in the random forest tuned on all rows.
  - `learning_curve`: the R2 on the training rows and on the held-out fold over the number of training rows, for the hyperparameters of that random forest and the chromosome folds.
  - `oob_error`: 1 - out-of-bag R2 over the number of trees, for the same hyperparameters on all rows.
  - The two curves come from extra fits for the figures. They do not change the outputs above.

## validate_model.py

Computes the numbers of the section "Model validation" in the main README. It tests `src/nmd_scanner/nmd_efficiency_rf.onnx` without retraining on TCGA out of fold, Geuvadis, GTEx and MMRF-TARGET.

As a script:

```bash
uv run --group scripts python scripts/validate_model.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa \
    --train-dir train/ --out-dir validate/
```

As a notebook, rendered to `validate/validate_model.html`:

```bash
uv run --group scripts jupytext --to ipynb scripts/validate_model.py &&
    TQDM_DISABLE=1 VALIDATE_MODEL_ARGS="--gff3 $PWD/gencode.v42.annotation.gff3.gz --fasta $PWD/GRCh38.fa --train-dir $PWD/train --out-dir $PWD/validate" \
    uv run --group scripts jupyter nbconvert --to html --execute scripts/validate_model.ipynb --output-dir validate/
```

- Inputs: the GFF3 and the FASTA, as above. `--train-dir` is the `--out-dir` of the `train_model.py` run that trained the model. The TCGA row uses its out-of-fold predictions. The script fails if the model of that run differs from `nmd_efficiency_rf.onnx`.
- Data, with a sha256 check for each file:
  - `MMRF_TARGET_dataset.csv` of NMDEff, at the commit that `make_benchmark_vcfs.py` pins.
  - Supplementary Data 2 of Kim et al. 2024 (GTEx).
  - 2 tables of Zenodo record 16666299 (Geuvadis). The record is one zip of 3.2 GB, so the script reads only those 2 members with HTTP range requests. `--geuvadis-dir` reads them from a local directory instead, with the same sha256 check.
- It builds one VCF of the variants and runs `nmd_scanner.annotate`. For GTEx and MMRF-TARGET, it converts HGVSc to genomic alleles with `make_benchmark_vcfs.py`.
- It keeps the rows the model can predict that have a variant allele fraction above 0 in the RNA, and it logs the rows it drops. For Geuvadis, the NMD efficiency is -log2(VAF_RNA / 0.5). GTEx and MMRF-TARGET use the column `NMD_efficiency` of their tables.
- The 95% CIs come from 2,000 bootstrap resamples of the genes per cohort. A resample draws genes with replacement and takes all rows of each drawn gene, because rows of one gene share a transcript structure.
- Outputs: `validation.json` in `--out-dir`, a short summary in the log, and the figures. Per cohort, `validation.json` has the numbers of rows and genes, the Spearman correlations of the model and of the rules, the gain of the model over the rules, and the mean measured NMD efficiency of the rows with `nmd_escape` False. For GTEx, it also has that mean over the variant-tissue rows, in total and per tissue. For MMRF-TARGET, `mmrf_target_by_transcript` has the correlations and the gain also for the rows on transcripts of the TCGA training rows, and for the other rows. The README table is written by hand from these numbers.
- Figures:
  - `measured_vs_predicted`: per cohort, the measured NMD efficiency over the prediction, with the Spearman correlation.
  - `calibration`: per cohort, the calibration plot of `calibration_oof`. Its least-squares line is the linear recalibration that the main README suggests.
  - `efficiency_by_nmd_escape`: per cohort, the measured NMD efficiency of the rows with `nmd_escape` False and True.
  - `spearman_per_cohort`: per cohort, the Spearman correlations of the model and of the rules with their 95% CIs, which are the numbers of the README table.
  - `gtex_tissues`: per GTEx tissue with at least 100 variant-tissue rows, the mean measured NMD efficiency of the rows with `nmd_escape` False, with its 95% CI.
  - `mmrf_target_by_transcript`: the Spearman correlations and the gain in MMRF-TARGET, for the rows on transcripts of the TCGA training rows and for the other rows. These are the numbers of `mmrf_target_by_transcript` in `validation.json`.
