# /// script
# requires-python = "==3.12.*"
# dependencies = [
#     "nmd-scanner==0.4.0",
#     "scikit-learn==1.9.1",
#     "lightgbm==4.7.0",
#     "pandas==3.0.6",
#     "numpy==2.5.3",
#     "scipy==1.18.1",
#     "joblib==1.6.0",
#     "pyarrow==24.0.0",
#     "pyfaidx==0.9.0.4",
#     "skl2onnx==1.20.0",
#     "onnx==1.23.2",
#     "onnxruntime==1.30.0",
# ]
# [tool.uv]
# exclude-newer = "2026-10-07T16:00:00Z"
# ///
"""
Train the NMD efficiency model on NMD-Scanner features of the NMDEff TCGA benchmark, and compare two random forests
with LightGBM in nested cross-validation, grouped by chromosome.

Run it with:

    uv run scripts/train_model.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa --out-dir out/

Steps:
1. Download tcga_dataset.csv of NMDEff (https://github.com/hjkng/nmdeff) at a pinned commit and check its sha256.
   Its column NMD_efficiency is the target.
2. Read resources/TCGA_benchmark/tcga_dataset.vcf. It has one row per TCGA row, in the same order, with the genomic
   alleles. Check that each REF matches the FASTA, and write the unique variants to <out-dir>/tcga_variants.vcf.
3. Run nmd_scanner.annotate on that VCF and keep the rows of the TCGA transcripts. NMD-Scanner's start is 0-based,
   so start + 1 is the TCGA start. Join on the transcript ID without version, chrom, start, end, ref and alt.
4. Keep the rows with nmd_model_status "ok". The inputs are nmd_scanner.schema.MODEL_INPUTS as float64.
5. Nested cross-validation of 3 models: a random forest with the hyperparameters of the former best_model.pkl, a
   random forest tuned with the grid of scripts/train_new.ipynb, and LightGBM with its default hyperparameters. The
   outer loop has 5 folds. The inner loop, which tunes the random forest, splits the training rows of an outer fold
   the same way.
   - Main grouping by chromosome: StratifiedGroupKFold over the chromosomes, stratified by 5 quantile bins of the
     target. The folds then have similar sizes and target distributions. Genes can overlap (antisense, nested), so a
     grouping by gene could still leak.
   - Comparison grouped by variant (chrom, pos, ref, alt): GroupKFold. It is close to the ungrouped folds of the
     notebook, but keeps the TCGA rows of one variant in one fold.
6. Tune the random forest on all usable rows with folds grouped by chromosome, and refit it. Fit LightGBM with its
   default hyperparameters on all usable rows. Save both to <out-dir>/models/. Also save the random forest as ONNX
   (see save_onnx), and check that its predictions equal those of the random forest, rounded to float32.

Outputs in --out-dir: tcga_dataset.csv, tcga_variants.vcf, tcga_features.parquet (the annotate result without
sequences), training_rows.parquet, oof_predictions.parquet, cv_results.json, cv_results.md and models/.
"""

import argparse
import dataclasses
import hashlib
import importlib.metadata
import json
import logging
import os
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import joblib
import numpy as np
import onnx
import onnxruntime
import pandas as pd
import pyarrow.parquet as pq
from lightgbm import LGBMRegressor
from onnx import TensorProto, helper, numpy_helper
from pyfaidx import Fasta
from scipy.stats import spearmanr
from skl2onnx.common.tree_ensemble import add_tree_to_attribute_pairs, get_default_tree_regressor_attribute_pairs
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, root_mean_squared_error
from sklearn.model_selection import GridSearchCV, GroupKFold, StratifiedGroupKFold

import nmd_scanner
from nmd_scanner.schema import MODEL_INPUTS

logger = logging.getLogger("train_model")

# tcga_dataset.csv was last changed in this commit, which is also the tag v1.0 of NMDEff
TCGA_URL = "https://raw.githubusercontent.com/hjkng/nmdeff/08c92768fcb689236a833db6a2f2d9bcbe919f12/tcga_dataset.csv"
TCGA_SHA256 = "7458af15a699a643de3638114af99a17d7eb95b3f80bce187e45347b9c1dc3ca"
DEFAULT_VCF = Path(__file__).resolve().parent.parent / "resources" / "TCGA_benchmark" / "tcga_dataset.vcf"

TARGET = "NMD_efficiency"
SEED = 42
N_FOLDS = 5
# The chromosome folds are stratified by this many quantile bins of the target
N_TARGET_BINS = 5

# Hyperparameters of best_model.pkl, the model before nmd_efficiency_rf.onnx: the best estimator of the grid search in
# scripts/train_new.ipynb
OLD_RF_PARAMS = {
    "n_estimators": 300,
    "max_depth": 11,
    "max_features": 0.3,
    "min_samples_leaf": 5,
    "min_samples_split": 15,
}

# refined_grid5 of scripts/train_new.ipynb, the grid that found OLD_RF_PARAMS: 324 candidates
RF_GRID = {
    "n_estimators": [200, 300, 400],
    "max_depth": [5, 7, 9, 11],
    "max_features": ["sqrt", 0.2, 0.3],
    "min_samples_split": [5, 10, 15],
    "min_samples_leaf": [2, 3, 5],
}

# LightGBM with its default hyperparameters. deterministic makes a fit give the same model on every run.
LGBM_PARAMS = {"random_state": SEED, "deterministic": True, "verbosity": -1}

MODELS = ("rf_old_params", "rf_tuned", "lgbm_default")
GROUPINGS = ("chromosome", "variant")

# Train and test row indices of one fold
Split = tuple[np.ndarray, np.ndarray]


@dataclass(frozen=True)
class StepCount:
    step: str
    rows: int


@dataclass(frozen=True)
class TrainingData:
    rows: pd.DataFrame
    x: pd.DataFrame
    y: np.ndarray
    groups: dict[str, np.ndarray]


@dataclass(frozen=True)
class Fit:
    estimator: object
    params: dict
    inner_r2: float | None


@dataclass(frozen=True)
class FoldInfo:
    grouping: str
    fold: int
    n_rows: int
    target_mean: float
    # The test chromosomes of a chromosome fold, None for a variant fold
    chromosomes: list[str] | None


@dataclass(frozen=True)
class Metrics:
    r2: float
    spearman: float
    rmse: float


@dataclass(frozen=True)
class FoldResult:
    grouping: str
    model: str
    fold: int
    n_train: int
    n_test: int
    metrics: Metrics
    inner_r2: float | None
    params: dict


@dataclass(frozen=True)
class Summary:
    grouping: str
    model: str
    r2_mean: float
    r2_sd: float
    spearman_mean: float
    spearman_sd: float
    rmse_mean: float
    rmse_sd: float


@dataclass(frozen=True)
class FinalModel:
    model: str
    path: str
    metadata_path: str
    # The ONNX file of the random forest, None for LightGBM
    onnx_path: str | None
    params: dict
    inner_r2: float | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, type=Path, help="GENCODE GFF3 of GRCh38, optionally gzip-compressed")
    parser.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA with chr names and a .fai index")
    parser.add_argument("--out-dir", required=True, type=Path, help="directory for the outputs")
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF, help="TCGA benchmark VCF (default: %(default)s)")
    parser.add_argument(
        "--features",
        type=Path,
        help="tcga_features.parquet of an earlier run with the same inputs; skips nmd_scanner.annotate",
    )
    parser.add_argument("--n-jobs", type=int, default=os.cpu_count(), help="parallel fits (default: all cores)")
    return parser.parse_args()


def download_tcga(out_dir: Path) -> pd.DataFrame:
    path = out_dir / "tcga_dataset.csv"
    if not path.exists():
        logger.info("Downloading %s", TCGA_URL)
        with urllib.request.urlopen(TCGA_URL, timeout=60) as response:
            path.write_bytes(response.read())
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != TCGA_SHA256:
        raise ValueError(f"{path} has sha256 {digest}, expected {TCGA_SHA256}")
    return pd.read_csv(path)


def read_vcf_alleles(vcf: Path) -> pd.DataFrame:
    records = []
    with open(vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            records.append(fields[:5])
    return pd.DataFrame(records, columns=["chrom", "pos", "id", "ref", "alt"]).astype({"pos": "int64"})


def add_alleles(tcga: pd.DataFrame, vcf: pd.DataFrame, fasta_path: Path) -> pd.DataFrame:
    """
    Add ref and alt of the VCF row of the same index to each TCGA row. Fail if a VCF row differs in chrom or pos from
    its TCGA row, or if a REF does not match the FASTA.
    """
    if len(vcf) != len(tcga):
        raise ValueError(f"the VCF has {len(vcf)} rows, tcga_dataset.csv {len(tcga)}")
    same_site = (vcf["chrom"].to_numpy() == tcga["chromosome"].to_numpy()) & (
        vcf["pos"].to_numpy() == tcga["start"].to_numpy()
    )
    if not same_site.all():
        raise ValueError(f"{(~same_site).sum()} VCF rows differ in chrom or pos from the TCGA row of the same index")

    fasta = Fasta(str(fasta_path))
    fasta_ref = [
        str(fasta[chrom][pos - 1 : pos - 1 + len(ref)]).upper()
        for chrom, pos, ref in zip(vcf["chrom"], vcf["pos"], vcf["ref"], strict=True)
    ]
    mismatches = vcf["ref"].to_numpy() != np.array(fasta_ref)
    if mismatches.any():
        raise ValueError(f"{mismatches.sum()} VCF rows have a REF that does not match the FASTA")

    rows = tcga.copy()
    rows["tcga_row"] = np.arange(len(tcga))
    rows["ref"] = vcf["ref"].to_numpy()
    rows["alt"] = vcf["alt"].to_numpy()
    return rows


def write_vcf(rows: pd.DataFrame, path: Path) -> int:
    variants = rows[["chromosome", "start", "ref", "alt"]].drop_duplicates()
    with open(path, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in variants.itertuples(index=False):
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    return len(variants)


def annotate(vcf: Path, gff3: Path, fasta: Path, out_path: Path) -> pd.DataFrame:
    started = time.time()
    features = nmd_scanner.annotate(vcf, gff3, fasta, sequences=False)
    logger.info("nmd_scanner.annotate: %d rows in %.0f s", len(features), time.time() - started)
    pq.write_table(nmd_scanner.to_arrow(features), out_path)
    return features


def strip_version(ids: pd.Series) -> pd.Series:
    return ids.astype("string").str.split(".").str[0]


def join_features(rows: pd.DataFrame, features: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """
    Join each TCGA row to the NMD-Scanner row of its transcript and variant. NMD-Scanner's start is 0-based and its
    end is exclusive, so start + 1 is the TCGA start and end is the TCGA end.
    """
    scanner = features.assign(
        transcript=strip_version(features["transcript_id"]),
        gene=strip_version(features["gene_id"]),
        start_1based=features["start"].astype("int64") + 1,
        end=features["end"].astype("int64"),
        chrom=features["chrom"].astype(str),
        ref=features["ref"].astype(str),
        alt=features["alt"].astype(str),
        nmd_model_status=features["nmd_model_status"].astype(str),
        strand=features["strand"].astype(str),
    )
    key = ["transcript", "chrom", "start_1based", "end", "ref", "alt"]
    if scanner.duplicated(key).any():
        raise ValueError("NMD-Scanner returned more than one row for a transcript and variant")
    tcga = rows.rename(columns={"Transcript_ID": "transcript", "chromosome": "chrom", "start": "start_1based"})[
        ["tcga_row", "transcript", "chrom", "start_1based", "end", "ref", "alt", "Hugo_Symbol", TARGET]
    ]
    joined = tcga.merge(scanner.drop(columns=["start"]), on=key, how="inner", validate="many_to_one")

    # The join of scripts/train_new.ipynb: transcript, start and end only. A site with two alt alleles in one
    # transcript then gives two rows for each of its TCGA rows.
    old_style = tcga[["transcript", "start_1based", "end"]].merge(
        scanner[["transcript", "start_1based", "end"]], on=["transcript", "start_1based", "end"], how="inner"
    )
    unshifted = tcga[["transcript", "start_1based", "end"]].merge(
        scanner[["transcript", "start", "end"]].astype({"start": "int64"}),
        left_on=["transcript", "start_1based", "end"],
        right_on=["transcript", "start", "end"],
        how="inner",
    )
    transcripts_found = tcga["transcript"].isin(set(scanner["transcript"]))
    checks = {
        "join_on_transcript_start_end_only": len(old_style),
        "join_without_start_shift": len(unshifted),
        "tcga_rows_whose_transcript_has_no_scanner_row": int((~transcripts_found).sum()),
    }
    return joined, checks


def training_data(joined: pd.DataFrame) -> TrainingData:
    rows = joined[joined["nmd_model_status"] == "ok"].reset_index(drop=True)
    x = rows[MODEL_INPUTS].astype("float64")
    variant = rows["chrom"] + ":" + rows["start_1based"].astype(str) + ":" + rows["ref"] + ":" + rows["alt"]
    groups = {"chromosome": rows["chrom"].to_numpy(dtype=str), "variant": variant.to_numpy()}
    return TrainingData(rows=rows, x=x, y=rows[TARGET].to_numpy(dtype="float64"), groups=groups)


def make_splits(grouping: str, y: np.ndarray, groups: np.ndarray) -> list[Split]:
    """
    Split the rows into N_FOLDS folds. chromosome: StratifiedGroupKFold over the chromosomes, stratified by
    N_TARGET_BINS quantile bins of y. variant: GroupKFold over the variants.
    """
    placeholder = np.zeros(len(y))
    if grouping == "chromosome":
        bins = pd.qcut(y, N_TARGET_BINS, labels=False)
        splitter = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
        return list(splitter.split(placeholder, bins, groups))
    if grouping == "variant":
        splitter = GroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
        return list(splitter.split(placeholder, y, groups))
    raise ValueError(f"unknown grouping {grouping}")


def fit_model(model: str, x: pd.DataFrame, y: np.ndarray, groups: np.ndarray, grouping: str, n_jobs: int) -> Fit:
    """
    Fit one model. rf_tuned runs its grid search on folds of the given grouping, which make_splits builds from
    groups. The other models ignore grouping and groups.
    """
    if model == "rf_old_params":
        estimator = RandomForestRegressor(**OLD_RF_PARAMS, random_state=SEED, n_jobs=n_jobs)
        estimator.fit(x, y)
        return Fit(estimator=estimator, params=OLD_RF_PARAMS, inner_r2=None)
    if model == "lgbm_default":
        estimator = LGBMRegressor(**LGBM_PARAMS)
        estimator.fit(x, y)
        return Fit(estimator=estimator, params=LGBM_PARAMS, inner_r2=None)
    if model == "rf_tuned":
        search = GridSearchCV(
            RandomForestRegressor(random_state=SEED, n_jobs=1),
            RF_GRID,
            scoring="r2",
            cv=make_splits(grouping, y, groups),
            n_jobs=n_jobs,
        )
        search.fit(x, y)
        return Fit(estimator=search.best_estimator_, params=search.best_params_, inner_r2=float(search.best_score_))
    raise ValueError(f"unknown model {model}")


def metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Metrics:
    return Metrics(
        r2=float(r2_score(y_true, y_pred)),
        spearman=float(spearmanr(y_true, y_pred).statistic),
        rmse=float(root_mean_squared_error(y_true, y_pred)),
    )


def fold_infos(data: TrainingData, grouping: str, splits: list[Split]) -> list[FoldInfo]:
    chromosomes = data.groups["chromosome"]
    return [
        FoldInfo(
            grouping=grouping,
            fold=fold,
            n_rows=len(test),
            target_mean=float(data.y[test].mean()),
            chromosomes=sorted(set(chromosomes[test]), key=chromosome_order) if grouping == "chromosome" else None,
        )
        for fold, (_, test) in enumerate(splits)
    ]


def chromosome_order(chrom: str) -> tuple[int, str]:
    name = chrom.removeprefix("chr")
    if name.isdigit():
        return (int(name), "")
    return (100, name)


def nested_cv(
    data: TrainingData, grouping: str, splits: list[Split], model: str, n_jobs: int
) -> tuple[list[FoldResult], pd.DataFrame]:
    groups = data.groups[grouping]
    folds = []
    predictions = []
    for fold, (train, test) in enumerate(splits):
        started = time.time()
        fit = fit_model(model, data.x.iloc[train], data.y[train], groups[train], grouping, n_jobs)
        y_pred = fit.estimator.predict(data.x.iloc[test])
        result = FoldResult(
            grouping=grouping,
            model=model,
            fold=fold,
            n_train=len(train),
            n_test=len(test),
            metrics=metrics(data.y[test], y_pred),
            inner_r2=fit.inner_r2,
            params=fit.params,
        )
        logger.info(
            "%s by %s, fold %d: R2 %.4f, %.0f s", model, grouping, fold, result.metrics.r2, time.time() - started
        )
        folds.append(result)
        predictions.append(
            pd.DataFrame(
                {
                    "tcga_row": data.rows["tcga_row"].to_numpy()[test],
                    "grouping": grouping,
                    "model": model,
                    "fold": fold,
                    "y": data.y[test],
                    "prediction": y_pred,
                }
            )
        )
    return folds, pd.concat(predictions, ignore_index=True)


def summarize(folds: list[FoldResult]) -> Summary:
    r2 = np.array([f.metrics.r2 for f in folds])
    spearman = np.array([f.metrics.spearman for f in folds])
    rmse = np.array([f.metrics.rmse for f in folds])
    return Summary(
        grouping=folds[0].grouping,
        model=folds[0].model,
        r2_mean=float(r2.mean()),
        r2_sd=float(r2.std(ddof=1)),
        spearman_mean=float(spearman.mean()),
        spearman_sd=float(spearman.std(ddof=1)),
        rmse_mean=float(rmse.mean()),
        rmse_sd=float(rmse.std(ddof=1)),
    )


def save_onnx(estimator: RandomForestRegressor, path: Path) -> None:
    """
    Write the random forest as one ONNX TreeEnsembleRegressor node (ai.onnx.ml opset 3). The input "input" has type
    double and shape [N, 19], with the columns of MODEL_INPUTS in that order. By the ONNX specification, the output
    "variable" has type float and shape [N, 1]. The metadata key feature_names holds MODEL_INPUTS as a JSON list.

    The RandomForestRegressor converter of skl2onnx would store the split thresholds and leaf values as float32. This
    function uses the tree helpers of the converter, but keeps them as float64.
    """
    attrs = get_default_tree_regressor_attribute_pairs()
    attrs["n_targets"] = 1
    n_trees = len(estimator.estimators_)
    for tree_id, tree in enumerate(estimator.estimators_):
        add_tree_to_attribute_pairs(
            attrs,
            False,
            tree.tree_,
            tree_id,
            1 / n_trees,
            0,
            False,
            adjust_threshold_for_sklearn=True,
            dtype=np.float64,
        )
    attrs["nodes_values_as_tensor"] = numpy_helper.from_array(np.array(attrs.pop("nodes_values"), dtype=np.float64))
    attrs["target_weights_as_tensor"] = numpy_helper.from_array(np.array(attrs.pop("target_weights"), dtype=np.float64))
    del attrs["nodes_hitrates"]
    node = helper.make_node("TreeEnsembleRegressor", ["input"], ["variable"], domain="ai.onnx.ml", **attrs)
    graph = helper.make_graph(
        [node],
        "nmd_efficiency_rf",
        [helper.make_tensor_value_info("input", TensorProto.DOUBLE, [None, len(MODEL_INPUTS)])],
        [helper.make_tensor_value_info("variable", TensorProto.FLOAT, [None, 1])],
    )
    model = helper.make_model(
        graph, opset_imports=[helper.make_opsetid("", 15), helper.make_opsetid("ai.onnx.ml", 3)], ir_version=8
    )
    helper.set_model_props(model, {"feature_names": json.dumps(list(MODEL_INPUTS))})
    onnx.checker.check_model(model)
    onnx.save(model, path)


def check_onnx(estimator: RandomForestRegressor, path: Path, x: pd.DataFrame) -> float:
    """
    Fail if the ONNX model predicts other values than the random forest, rounded to float32. Return the largest
    absolute difference to the float64 predictions of the random forest.
    """
    session = onnxruntime.InferenceSession(path, providers=["CPUExecutionProvider"])
    onnx_pred = session.run(None, {"input": x.to_numpy(dtype=np.float64)})[0].ravel()
    rf_pred = estimator.predict(x)
    differs = onnx_pred != rf_pred.astype(np.float32)
    if differs.any():
        raise ValueError(f"the ONNX model differs from the random forest in {differs.sum()} of {len(x)} rows")
    return float(np.abs(onnx_pred.astype(np.float64) - rf_pred).max())


def save_final_model(model: str, fit: Fit, data: TrainingData, models_dir: Path) -> FinalModel:
    metadata = {
        "feature_names": list(MODEL_INPUTS),
        "target": TARGET,
        "params": fit.params,
        "inner_cv_r2_grouped_by_chromosome": fit.inner_r2,
        "n_training_rows": len(data.y),
        "nmd_scanner_version": importlib.metadata.version("nmd-scanner"),
    }
    onnx_path = None
    if model == "rf_tuned":
        path = models_dir / "nmd_efficiency_rf.joblib"
        joblib.dump(fit.estimator, path)
        metadata["sklearn_version"] = importlib.metadata.version("scikit-learn")
        onnx_path = models_dir / "nmd_efficiency_rf.onnx"
        save_onnx(fit.estimator, onnx_path)
        max_diff = check_onnx(fit.estimator, onnx_path, data.x)
        logger.info("%s predicts the training rows as the random forest, max abs diff %.3g", onnx_path, max_diff)
    elif model == "lgbm_default":
        path = models_dir / "nmd_efficiency_lgbm.txt"
        fit.estimator.booster_.save_model(path)
        metadata["lightgbm_version"] = importlib.metadata.version("lightgbm")
    else:
        raise ValueError(f"no final model for {model}")
    metadata_path = path.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return FinalModel(
        model=model,
        path=f"models/{path.name}",
        metadata_path=f"models/{metadata_path.name}",
        onnx_path=None if onnx_path is None else f"models/{onnx_path.name}",
        params=fit.params,
        inner_r2=fit.inner_r2,
    )


def start_loss_stats(data: TrainingData) -> dict:
    start_loss = data.rows["start_loss"].astype(bool).to_numpy()
    return {
        "ok_rows_with_start_loss": int(start_loss.sum()),
        "target_median_start_loss": float(np.median(data.y[start_loss])) if start_loss.any() else None,
        "target_median_other": float(np.median(data.y[~start_loss])),
    }


def markdown(
    row_counts: list[StepCount], infos: list[FoldInfo], summaries: list[Summary], finals: list[FinalModel]
) -> str:
    lines = ["# NMD efficiency model: nested grouped CV", "", "| step | rows |", "| --- | ---: |"]
    lines += [f"| {c.step} | {c.rows} |" for c in row_counts]
    lines += ["", "| chromosome fold | rows | target mean | chromosomes |", "| ---: | ---: | ---: | --- |"]
    lines += [
        f"| {i.fold} | {i.n_rows} | {i.target_mean:.3f} | {' '.join(i.chromosomes)} |"
        for i in infos
        if i.grouping == "chromosome"
    ]
    lines += [
        "",
        f"{N_FOLDS} outer folds; the inner loop splits the same way. Seed {SEED}. Mean and sd (ddof 1) over the outer"
        " folds.",
        "",
        "| grouping | model | R2 | Spearman | RMSE |",
        "| --- | --- | --- | --- | --- |",
    ]
    lines += [
        f"| {s.grouping} | {s.model} | {s.r2_mean:.3f} ± {s.r2_sd:.3f} | {s.spearman_mean:.3f} ± {s.spearman_sd:.3f}"
        f" | {s.rmse_mean:.3f} ± {s.rmse_sd:.3f} |"
        for s in summaries
    ]
    lines += ["", "| final model | inner CV R2 | params |", "| --- | --- | --- |"]
    lines += [
        f"| {f.model} | {'none' if f.inner_r2 is None else f'{f.inner_r2:.3f}'} | `{json.dumps(f.params)}` |"
        for f in finals
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    logger.setLevel(logging.INFO)
    out_dir = args.out_dir
    models_dir = out_dir / "models"
    models_dir.mkdir(parents=True, exist_ok=True)

    tcga = download_tcga(out_dir)
    vcf = read_vcf_alleles(args.vcf)
    rows = add_alleles(tcga, vcf, args.fasta)
    vcf_path = out_dir / "tcga_variants.vcf"
    n_variants = write_vcf(rows, vcf_path)

    if args.features is None:
        features = annotate(vcf_path, args.gff3, args.fasta, out_dir / "tcga_features.parquet")
    else:
        features = pd.read_parquet(args.features)
    joined, join_checks = join_features(rows, features)
    data = training_data(joined)
    key_columns = ["tcga_row", "chrom", "start_1based", "end", "ref", "alt", "transcript", "gene"]
    data.rows[[*key_columns, "strand", "Hugo_Symbol", TARGET, *MODEL_INPUTS]].to_parquet(
        out_dir / "training_rows.parquet"
    )

    row_counts = [
        StepCount("tcga_dataset.csv rows", len(tcga)),
        StepCount("unique variants in tcga_variants.vcf", n_variants),
        StepCount("NMD-Scanner rows, all transcripts", len(features)),
        StepCount("rows joined to the TCGA transcript", len(joined)),
        StepCount("rows with nmd_model_status ok", len(data.y)),
    ]
    status_counts = joined.groupby("nmd_model_status").size().to_dict()
    start_loss_by_status = joined[joined["start_loss"].astype("boolean").fillna(False)].groupby("nmd_model_status")
    for count in row_counts:
        logger.info("%s: %d", count.step, count.rows)
    logger.info("nmd_model_status of the joined rows: %s", status_counts)

    infos = []
    folds = []
    summaries = []
    oof = []
    for grouping in GROUPINGS:
        splits = make_splits(grouping, data.y, data.groups[grouping])
        infos += fold_infos(data, grouping, splits)
        for model in MODELS:
            model_folds, predictions = nested_cv(data, grouping, splits, model, args.n_jobs)
            folds += model_folds
            summaries.append(summarize(model_folds))
            oof.append(predictions)
    for info in infos:
        logger.info("%s fold %d: %d rows, target mean %.3f", info.grouping, info.fold, info.n_rows, info.target_mean)
    oof_predictions = pd.concat(oof, ignore_index=True)
    oof_predictions.to_parquet(out_dir / "oof_predictions.parquet")

    finals = []
    for model in ("rf_tuned", "lgbm_default"):
        fit = fit_model(model, data.x, data.y, data.groups["chromosome"], "chromosome", args.n_jobs)
        finals.append(save_final_model(model, fit, data, models_dir))
        logger.info("final %s: inner CV R2 %s, %s", model, fit.inner_r2, fit.params)

    results = {
        "versions": {
            name: importlib.metadata.version(name)
            for name in ("nmd-scanner", "scikit-learn", "lightgbm", "pandas", "numpy", "scipy")
        },
        "inputs": {
            "tcga_url": TCGA_URL,
            "tcga_sha256": TCGA_SHA256,
            "vcf": args.vcf.name,
            "gff3": args.gff3.name,
            "fasta": args.fasta.name,
        },
        "row_counts": [dataclasses.asdict(c) for c in row_counts],
        "join_checks": join_checks,
        "status_counts": {str(k): int(v) for k, v in status_counts.items()},
        "start_loss_rows_by_status": {str(k): int(v) for k, v in start_loss_by_status.size().items()},
        "start_loss": start_loss_stats(data),
        "cv": {
            "splitters": {
                "chromosome": f"StratifiedGroupKFold(n_splits={N_FOLDS}, shuffle=True, random_state={SEED}), y in"
                f" {N_TARGET_BINS} quantile bins, groups chrom",
                "variant": f"GroupKFold(n_splits={N_FOLDS}, shuffle=True, random_state={SEED}), groups variant",
            },
            "fold_assignment": {
                chrom: info.fold for info in infos if info.grouping == "chromosome" for chrom in info.chromosomes
            },
            "folds_info": [dataclasses.asdict(i) for i in infos],
            "rf_old_params": OLD_RF_PARAMS,
            "rf_grid": RF_GRID,
            "lgbm_params": LGBM_PARAMS,
            "summaries": [dataclasses.asdict(s) for s in summaries],
            "folds": [dataclasses.asdict(f) for f in folds],
        },
        "final_models": [dataclasses.asdict(f) for f in finals],
    }
    (out_dir / "cv_results.json").write_text(json.dumps(results, indent=2) + "\n")
    (out_dir / "cv_results.md").write_text(markdown(row_counts, infos, summaries, finals))
    logger.info("Wrote %s", out_dir / "cv_results.md")


if __name__ == "__main__":
    main()
