# Scripts

The scripts build the benchmarks, train the NMD efficiency model and validate it. They need the extra `scripts`, and they run from the repo root:

```bash
uv run --extra scripts python scripts/build_benchmarks.py
uv run --extra scripts python scripts/train_model.py
uv run --extra scripts python scripts/validate_model.py.py
```

They read and write `data/` by default: the GFF3 `data/gencode.v42.annotation.gff3.gz` and the FASTA `data/GRCh38.fa` with its `.fai` index go in, and `data/benchmarks/`, `data/train/` and `data/validate/` come out.

- `build_benchmarks.py` downloads the study tables at pinned versions and writes one table per dataset of `nmd_scanner.benchmark.DATASETS`. A new dataset needs one module in `nmd_scanner.benchmark.datasets` and one `DATASETS` entry.
- `train_model.py` trains the random forest and writes it to `data/train/nmd_efficiency_rf.onnx`, and as a joblib pickle to `data/train/nmd_efficiency_rf.pkl`. To ship the model, copy the ONNX to `src/nmd_scanner/` and the pickle to the repo root. `validate_model.py.py` fails if the ONNX in `data/train/` differs from the one in `src/nmd_scanner/`.
- `validate_model.py.py` writes `validation.json`, and it draws the figures of the training and the validation into `data/validate/figures/`. It sets its paths in its first code cell. The main README shows one of these figures, and after a rerun `git add -f` updates that PNG.

`validate_model.py.py` is also a jupytext notebook. To render it to HTML:

```bash
TQDM_DISABLE=1 uv run --extra scripts jupytext --to ipynb --execute --run-path . scripts/validate_model.py.py
uv run --extra scripts jupyter nbconvert --to html scripts/validate_model.py.ipynb --output-dir data/validate
```

The study tables follow the terms of their sources: NMDEff (Kim et al. 2024, <https://github.com/hjkng/nmdeff>) for TCGA and MMRF-TARGET, Supplementary Data 2 of Kim et al. 2024 for GTEx v8, and Zenodo record 16666299 (Iha et al. 2025) for Geuvadis.
