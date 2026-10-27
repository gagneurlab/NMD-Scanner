# Scripts

The scripts build the benchmarks, train the NMD efficiency model and validate it. They run on the checkout and need the extra `scripts`:

```bash
uv run --extra scripts python scripts/build_benchmarks.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa --out-dir benchmarks/
uv run --extra scripts python scripts/train_model.py --benchmark-dir benchmarks/ --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa --out-dir train/
uv run --extra scripts python scripts/validate_model.py --benchmark-dir benchmarks/ --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa --train-dir train/ --out-dir validate/
```

- `build_benchmarks.py` downloads the study tables at pinned versions and writes one table per dataset of `nmd_scanner.benchmark.DATASETS`. A new dataset needs one module in `nmd_scanner.benchmark.datasets` and one `DATASETS` entry.
- `train_model.py` trains the random forest that the package ships as `src/nmd_scanner/nmd_efficiency_rf.onnx`. It also saves the forest as a joblib pickle, a copy of which is `nmd_efficiency_rf.pkl` at the repo root.
- `validate_model.py` writes the numbers of the "Model validation" table of the main README to `validation.json`.
- Each script saves its figures as PNG and PDF in `<out-dir>/figures/`.

`train_model.py` and `validate_model.py` are also jupytext notebooks. To render one to HTML, pass its arguments with absolute paths in `TRAIN_MODEL_ARGS` or `VALIDATE_MODEL_ARGS`:

```bash
uv run --extra scripts jupytext --to ipynb scripts/train_model.py &&
    TQDM_DISABLE=1 TRAIN_MODEL_ARGS="--benchmark-dir $PWD/benchmarks --gff3 $PWD/gencode.v42.annotation.gff3.gz --fasta $PWD/GRCh38.fa --out-dir $PWD/train" \
    uv run --extra scripts jupyter nbconvert --to html --execute scripts/train_model.ipynb --output-dir train/
```

The study tables follow the terms of their sources: NMDEff (Kim et al. 2024, <https://github.com/hjkng/nmdeff>) for TCGA and MMRF-TARGET, Supplementary Data 2 of Kim et al. 2024 for GTEx v8, and Zenodo record 16666299 (Iha et al. 2025) for Geuvadis.
