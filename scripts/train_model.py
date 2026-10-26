# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
"""
Train the NMD efficiency model nmd_efficiency_rf.onnx on NMD-Scanner features of the NMDEff TCGA benchmark, and
predict each training row out of fold.

Run it as a script with:

    uv run --group scripts python scripts/train_model.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa \
        --out-dir train/

It is also a notebook (see scripts/README.md). In a notebook, the arguments come from the environment variable
TRAIN_MODEL_ARGS.

Steps:
1. Download tcga_dataset.csv of NMDEff and add the genomic alleles, with the steps of make_benchmark_vcfs.py. Its
   column NMD_efficiency is the target.
2. Run nmd_scanner.annotate on the unique variants. Join each TCGA row to the NMD-Scanner row of its transcript and
   variant, and keep the rows with nmd_efficiency_prediction_status "ok". The inputs are
   nmd_scanner.schema.MODEL_INPUTS as float64.
3. Split the rows into 5 folds with StratifiedGroupKFold over the chromosomes, stratified by 5 quantile bins of the
   target. The folds then have similar sizes and target distributions. Genes can overlap (antisense, nested), so a
   grouping by gene could still leak.
4. For each fold, tune a random forest on the other folds, and predict the rows of the fold. The grid search of the
   tuning splits its rows the same way. validate_model.py reads these out-of-fold predictions.
5. Tune the random forest the same way on all rows, and save it as ONNX (see save_onnx). Fail if the ONNX model
   predicts other values than the random forest, rounded to float32.

Outputs in --out-dir: variants.vcf (the input of nmd_scanner.annotate), oof_predictions.parquet,
nmd_efficiency_rf.onnx, cv_results.json (the R2 of each fold and the hyperparameters), and figures/ with each figure
as PNG and PDF.
"""

import argparse
import json
import logging
import os
import shlex
from pathlib import Path

import numpy as np
import onnx
import onnxruntime
import pandas as pd
from figures import (
    AXIS,
    BLUE,
    MARKER_STYLE,
    MUTED_INK,
    SECONDARY_INK,
    SURFACE,
    in_notebook,
    new_figure,
    plot_calibration,
    plot_measured_vs_predicted,
    save_figure,
)
from make_benchmark_vcfs import KEY, add_alleles, add_key, read_nmdeff, read_strands, write_vcf
from onnx import TensorProto, helper, numpy_helper
from pyfaidx import Fasta
from scipy.stats import spearmanr
from skl2onnx.common.tree_ensemble import add_tree_to_attribute_pairs, get_default_tree_regressor_attribute_pairs
from sklearn.base import clone
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, learning_curve

import nmd_scanner
from nmd_scanner.schema import MODEL_INPUTS

logger = logging.getLogger("train_model")

# %%
SEED = 42
N_FOLDS = 5
# The folds are stratified by this many quantile bins of the target
N_TARGET_BINS = 5
# The grid of the random forest search: 324 candidates
RF_GRID = {
    "n_estimators": [200, 300, 400],
    "max_depth": [5, 7, 9, 11],
    "max_features": ["sqrt", 0.2, 0.3],
    "min_samples_split": [5, 10, 15],
    "min_samples_leaf": [2, 3, 5],
}
# The permutation importance permutes each input this many times
N_PERMUTATIONS = 10
# The learning curve fits the final hyperparameters on these shares of the training rows of the fold with the fewest,
# in every fold
LEARNING_CURVE_SIZES = np.linspace(0.1, 1.0, 8)
# The out-of-bag error curve fits the final hyperparameters with these numbers of trees
OOB_TREES = [50, 100, 200, 300, 400, 500]


def chromosome_folds(y: np.ndarray, chrom: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split the rows into N_FOLDS folds of whole chromosomes, stratified by N_TARGET_BINS quantile bins of y."""
    splitter = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    return list(splitter.split(np.zeros(len(y)), pd.qcut(y, N_TARGET_BINS, labels=False), chrom))


def tune(x: pd.DataFrame, y: np.ndarray, chrom: np.ndarray) -> GridSearchCV:
    """Run a grid search over RF_GRID on chromosome_folds, and refit the best random forest on all given rows."""
    forest = RandomForestRegressor(random_state=SEED, n_jobs=1)
    search = GridSearchCV(forest, RF_GRID, scoring="r2", cv=chromosome_folds(y, chrom), n_jobs=-1)
    return search.fit(x, y)


def save_onnx(estimator: RandomForestRegressor, path: Path) -> None:
    """
    Write the random forest as one ONNX TreeEnsembleRegressor node (ai.onnx.ml opset 3). The input "input" has type
    double and shape [N, len(MODEL_INPUTS)], with the columns of MODEL_INPUTS in that order. By the ONNX specification, the output
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


# %%
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--gff3", required=True, type=Path, help="GENCODE GFF3 of GRCh38, optionally gzip-compressed")
parser.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA with chr names and a .fai index")
parser.add_argument("--out-dir", required=True, type=Path, help="directory for the outputs")
# In a notebook, sys.argv holds the arguments of the kernel
args = parser.parse_args(shlex.split(os.environ.get("TRAIN_MODEL_ARGS", "")) if in_notebook() else None)
logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s")
logger.setLevel(logging.INFO)
args.out_dir.mkdir(parents=True, exist_ok=True)

# %%
tcga = add_alleles(read_nmdeff("tcga_dataset.csv"), read_strands(args.gff3), Fasta(str(args.fasta)))
vcf = args.out_dir / "variants.vcf"
write_vcf(tcga.drop_duplicates(["chrom", "pos", "ref", "alt"]), vcf)
features = add_key(nmd_scanner.annotate(vcf, args.gff3, args.fasta, sequences=False, prediction=False))
# tcga_row is the index of the row in tcga_dataset.csv
tcga = tcga[[*KEY, "NMD_efficiency"]].assign(tcga_row=np.arange(len(tcga)))
rows = tcga.merge(features, on=KEY, how="inner", validate="many_to_one")
rows = rows[rows["nmd_efficiency_prediction_status"] == "ok"].reset_index(drop=True)
logger.info("%d TCGA rows, %d of them with nmd_efficiency_prediction_status ok", len(tcga), len(rows))
x = rows[MODEL_INPUTS].astype("float64")
y = rows["NMD_efficiency"].to_numpy(dtype="float64")
chrom = rows["chrom"].to_numpy(dtype=str)

# %% [markdown]
# The histogram shows the target `NMD_efficiency` of the training rows. 0 means no NMD, and 1 means that the variant
# allele fraction in the RNA is half of that in the DNA.

# %%
fig, (ax,) = new_figure(6, 3.5)
ax.hist(y, bins=60, color=BLUE, edgecolor=SURFACE, linewidth=0.5)
ax.grid(axis="x", visible=False)
ax.set_title(f"NMD_efficiency of {len(y):,} TCGA rows")
ax.set_xlabel("NMD_efficiency")
ax.set_ylabel("Rows")
save_figure(fig, args.out_dir, "target_histogram")

# %%
# The rows of each fold in a block, in fold order, as validate_model.py expects them
predictions = []
folds = []
for fold, (train, test) in enumerate(chromosome_folds(y, chrom)):
    search = tune(x.iloc[train], y[train], chrom[train])
    prediction = search.predict(x.iloc[test])
    columns = ["tcga_row", *KEY, "nmd_escape", "NMD_efficiency"]
    predictions.append(rows.iloc[test][columns].assign(fold=fold, prediction=prediction))
    chromosomes = sorted(set(chrom[test].tolist()))
    r2 = r2_score(y[test], prediction)
    folds.append({"chromosomes": chromosomes, "rows": len(test), "r2": r2, "params": search.best_params_})
    logger.info("fold %d (%s): R2 %.3f, %s", fold, " ".join(chromosomes), r2, search.best_params_)
oof = pd.concat(predictions)
oof.to_parquet(args.out_dir / "oof_predictions.parquet", index=False)

# %% [markdown]
# The scatter shows the measured NMD efficiency of each TCGA row over its out-of-fold prediction, with Spearman and R2
# over all rows. The gray line marks equal values.

# %%
fig, (ax,) = new_figure(5, 4.5)
oof_spearman = spearmanr(oof["prediction"], oof["NMD_efficiency"]).statistic
oof_r2 = r2_score(oof["NMD_efficiency"], oof["prediction"])
plot_measured_vs_predicted(
    ax,
    oof["prediction"].to_numpy(),
    oof["NMD_efficiency"].to_numpy(),
    f"TCGA, out of fold, {len(oof):,} rows\nSpearman {oof_spearman:.2f}, R2 {oof_r2:.2f}",
)
save_figure(fig, args.out_dir, "tcga_measured_vs_predicted")

# %% [markdown]
# The dots show the mean measured over the mean predicted NMD efficiency in each quantile bin of the out-of-fold
# prediction, with the 95% CI of the measured mean. The blue line is the least-squares fit of measured on predicted over all rows,
# and the gray line marks equal values.

# %%
fig, (ax,) = new_figure(5, 6)
plot_calibration(
    ax, oof["prediction"].to_numpy(), oof["NMD_efficiency"].to_numpy(), f"TCGA, out of fold, {len(oof):,} rows"
)
fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center")
save_figure(fig, args.out_dir, "calibration_oof")

# %%
final = tune(x, y, chrom)
onnx_path = args.out_dir / "nmd_efficiency_rf.onnx"
save_onnx(final.best_estimator_, onnx_path)
session = onnxruntime.InferenceSession(onnx_path, providers=["CPUExecutionProvider"])
onnx_prediction = session.run(None, {"input": x.to_numpy(dtype=np.float64)})[0].ravel()
if (onnx_prediction != final.predict(x).astype(np.float32)).any():
    raise ValueError(f"{onnx_path} predicts other values than the random forest")

results = {"rows": len(rows), "folds": folds, "inner_cv_r2": final.best_score_, "params": final.best_params_}
(args.out_dir / "cv_results.json").write_text(json.dumps(results, indent=2) + "\n")
logger.info("final model: inner CV R2 %.3f, %s", final.best_score_, final.best_params_)

# %% [markdown]
# The bars show the impurity importance and the permutation importance of the model inputs in the final random forest.
# The permutation importance is the loss in R2 on the training rows when one input is permuted, with the standard
# deviation over the permutations.

# %%
impurity = pd.Series(final.best_estimator_.feature_importances_, index=MODEL_INPUTS).sort_values()
permutation = permutation_importance(
    final.best_estimator_, x, y, scoring="r2", n_repeats=N_PERMUTATIONS, random_state=SEED, n_jobs=-1
)
permutation_order = np.argsort(permutation.importances_mean)
fig, (impurity_ax, permutation_ax) = new_figure(9, 5, ncols=2)
impurity_ax.barh(impurity.index, impurity.to_numpy(), height=0.6, color=BLUE)
impurity_ax.set_title("Impurity importance")
impurity_ax.set_xlabel("Mean decrease in impurity")
permutation_ax.barh(
    np.array(MODEL_INPUTS)[permutation_order],
    permutation.importances_mean[permutation_order],
    xerr=permutation.importances_std[permutation_order],
    height=0.6,
    color=BLUE,
    error_kw={"ecolor": SECONDARY_INK, "elinewidth": 0.75},
)
permutation_ax.set_title("Permutation importance")
permutation_ax.set_xlabel(f"Decrease in R2 on the training rows,\nmean of {N_PERMUTATIONS} permutations")
for ax in (impurity_ax, permutation_ax):
    ax.grid(axis="y", visible=False)
    ax.margins(y=0.01)
save_figure(fig, args.out_dir, "feature_importance")

# %% [markdown]
# The learning curve shows the R2 of the final hyperparameters on the training rows and on the held-out fold, for
# growing numbers of training rows, with bands of one standard deviation over the folds. Each chromosome fold trains on
# the same numbers of rows, up to all training rows of the fold with the fewest.

# %%
learning_folds = chromosome_folds(y, chrom)
# learning_curve would scale shares by the training rows of the first fold only
learning_sizes = (LEARNING_CURVE_SIZES * min(len(train) for train, _ in learning_folds)).astype(int)
# A copy of the final random forest with its hyperparameters, unfitted. The final model stays as it is.
train_sizes, train_scores, test_scores = learning_curve(
    clone(final.best_estimator_),
    x,
    y,
    cv=learning_folds,
    train_sizes=learning_sizes,
    scoring="r2",
    shuffle=True,
    random_state=SEED,
    n_jobs=-1,
)
fig, (ax,) = new_figure(6, 4)
for scores, label, color in ((train_scores, "Training rows", MUTED_INK), (test_scores, "Held-out fold", BLUE)):
    mean = scores.mean(axis=1)
    sd = scores.std(axis=1)
    ax.fill_between(train_sizes, mean - sd, mean + sd, color=color, alpha=0.1, linewidth=0)
    ax.plot(train_sizes, mean, marker="o", color=color, **MARKER_STYLE)
    ax.annotate(
        label, (train_sizes[-1], mean[-1]), xytext=(8, 0), textcoords="offset points", va="center", color=SECONDARY_INK
    )
ax.set_xlim(right=train_sizes[-1] * 1.2)
ax.set_title("Learning curve of the final hyperparameters")
ax.set_xlabel("Training rows")
ax.set_ylabel("R2")
save_figure(fig, args.out_dir, "learning_curve")

# %% [markdown]
# The curve shows the out-of-bag error of the final hyperparameters, fitted on all rows with a growing number of trees.
# The out-of-bag prediction of a row is the mean prediction of the trees that did not draw it, and the out-of-bag error
# is 1 - R2 of these predictions.

# %%
oob_errors = [
    1 - clone(final.best_estimator_).set_params(n_estimators=trees, oob_score=True, n_jobs=-1).fit(x, y).oob_score_
    for trees in OOB_TREES
]
fig, (ax,) = new_figure(6, 4)
ax.plot(OOB_TREES, oob_errors, marker="o", color=BLUE, **MARKER_STYLE)
final_trees = final.best_params_["n_estimators"]
ax.axvline(final_trees, color=AXIS, linewidth=1)
ax.annotate(
    f"Final model, {final_trees} trees",
    (final_trees, 1),
    xycoords=("data", "axes fraction"),
    xytext=(4, -2),
    textcoords="offset points",
    va="top",
    color=SECONDARY_INK,
)
ax.set_title("Out-of-bag error of the final hyperparameters")
ax.set_xlabel("Trees")
ax.set_ylabel("1 - out-of-bag R2")
save_figure(fig, args.out_dir, "oob_error")
