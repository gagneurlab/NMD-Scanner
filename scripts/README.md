# Scripts and Notebooks

This folder contains Jupyter notebooks used for development, testing, exploratory
analysis, and model training for the NMD-Scanner project. These notebooks are
**not part of the Python package** but document the code used during the research and
implementation process.

## Contents

### 1. `create_test_VCF.ipynb`

Creates synthetic VCF files containing possible edge-case variants used to test the robustness
of the function incorporating variants into the sequence.

The notebook generates variants such as:

- start codon substitutions, deletions, insertions, duplications
- stop codon substitutions, length-changing variants at the end of the CDS
- nonsense mutations
- frameshift and in-frame indels
- variants before and after the CDS
- complex or multi-base substitutions

---

### 2. `train.ipynb`

Notebook used to train predictive models to evaluate NMD efficiency. The workflow includes:

#### **(A) Benchmarking against NMDEff**

- Runs the original NMDEff implementation (<https://github.com/hjkng/nmdeff>)
- Extracts NMD efficiency scores using their official scripts
- Prepares the same variants using the NMD-Scanner to compare feature sets and predictions

#### **(B) Model exploration**

Multiple machine-learning models were tested, including:

- Random Forest
- XGBoost
- LightGBM
- Gradient Boosting
- Ridge
- Lasso
- SVR

#### **(C) Final model selection**

- After evaluation, the Random Forest Regressor gave the best performance based on cross-validation metrics.
- The trained model is saved (joblib) for later prediction.

---

### 3. `validation_MMRF_TARGET.ipynb`

Evaluation of the trained model using independent datasets (MMRF/TARGET, <https://github.com/hjkng/nmdeff>).
It includes:

- loading validation datasets
- generating NMD-related features
- comparing predictions with reference measurements
- producing evaluation plots and performance metrics

---

### 4. `nmd-vep.ipynb`

This notebook contains the full NMD-Scanner implementation split into individual cells, allowing users to:

- follow the workflow step-by-step
- inspect intermediate outputs
- debug or modify specific steps of the pipeline
- reuse, adapt or extend the code for custom analyses.

---

### 5. `train_model.py`

Retrains the NMD efficiency model on the features of the installed NMD-Scanner release. It replaces the
notebooks above for that purpose: it runs end to end and pins its dependencies in a PEP 723 header.

```bash
uv run scripts/train_model.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa --out-dir out/
```

- The target is `NMD_efficiency` of the NMDEff TCGA benchmark, downloaded at a pinned commit.
- The variants are `resources/TCGA_benchmark/tcga_dataset.vcf`.
- It compares a random forest with the hyperparameters of the former `best_model.pkl`, a tuned random forest
  and LightGBM with default hyperparameters in nested cross-validation. The folds are grouped by chromosome, and
  for comparison by variant.
- It saves the tuned random forest and the LightGBM model, both fit on all usable rows, to `out/models/`.
- It also saves the tuned random forest as ONNX, to `out/models/nmd_efficiency_rf.onnx`, and checks that the ONNX
  model predicts the training rows like the random forest. This file is `nmd_efficiency_rf.onnx` at the root of the
  repository. The main README shows how to load it.

---

### 6. `make_tcga_vcf.py`

Builds `resources/TCGA_benchmark/tcga_dataset.vcf` from the NMDEff study table, so that the VCF is reproducible.
The output equals the committed file byte for byte.

```bash
uv run scripts/make_tcga_vcf.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa
```

- It downloads `tcga_dataset.csv` at the commit that `train_model.py` pins and checks its sha256. `--csv` reads a
  local copy instead, with the same check.
- It parses the substitution from HGVSc, which has the alleles on the transcript strand. It takes the strand from
  the GENCODE GFF3 and complements REF and ALT on the minus strand. For transcripts missing from the GFF3, it uses
  the orientation whose REF matches the FASTA.
- It fails if a row is not a single-base substitution or if a REF does not match the FASTA.

---

## Notes

- None of the notebooks are required for end users of the package.
- They are included for transparency and reproducibility of the development process.
- The recommended entry point for users is the installed nmd_scanner package or the CLI interface described in the main README.
