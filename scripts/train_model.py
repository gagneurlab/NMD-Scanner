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
Train the NMD efficiency model nmd_efficiency_rf.onnx on NMD-Scanner features of the TCGA benchmark, and predict each
training row out of fold.

Run it as a script with:

    uv run --extra scripts python scripts/train_model.py --benchmark-dir benchmarks/ \
        --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa --out-dir train/

--benchmark-dir is the --out-dir of scripts/build_benchmarks.py, with the same GFF3 and FASTA.

It is also a notebook (see scripts/README.md). In a notebook, the arguments come from the environment variable
TRAIN_MODEL_ARGS.

Steps:
1. Read tcga.parquet of --benchmark-dir (see nmd_scanner.benchmark.datasets.tcga). Its column NMD_efficiency is the
   target.
2. Run nmd_scanner.annotate on the unique variants. Join each TCGA row to the NMD-Scanner row of its transcript and
   variant, and keep the rows with nmd_efficiency_prediction_status "ok". Log the rows that this drops, and why. The
   inputs are nmd_scanner.schema.MODEL_INPUTS as float64.
3. Split the rows into 5 folds with StratifiedGroupKFold over the chromosomes, stratified by 5 quantile bins of the
   target. The folds then have similar sizes and target distributions. Genes can overlap (antisense, nested), so a
   grouping by gene could still leak.
4. For each fold, tune a random forest on the other folds, and predict the rows of the fold. The grid search of the
   tuning splits its rows the same way. validate_model.py reads these out-of-fold predictions.
5. Tune the random forest the same way on all rows, and save it as ONNX (see save_onnx). Fail if the ONNX model
   predicts other values than the random forest, rounded to float32. Save the same random forest as a joblib pickle.

Outputs in --out-dir: variants.vcf (the input of nmd_scanner.annotate), oof_predictions.parquet,
nmd_efficiency_rf.onnx, nmd_efficiency_rf.pkl, cv_results.json (the R2 of each fold and the hyperparameters), and
figures/ with each figure as PNG and PDF.
"""

import argparse
import json
import logging
import os
import shlex
from pathlib import Path

import joblib
import numpy as np
import onnx
import onnxruntime
import plotnine as pn
import polars as pl
from onnx import TensorProto, helper, numpy_helper
from scipy.stats import spearmanr
from skl2onnx.common.tree_ensemble import add_tree_to_attribute_pairs, get_default_tree_regressor_attribute_pairs
from sklearn.base import clone
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, learning_curve

from nmd_scanner.benchmark import KEY, join_features, run_nmd_scanner
from nmd_scanner.benchmark.figures import (
    AXIS,
    BLUE,
    MUTED_INK,
    SECONDARY_INK,
    SURFACE,
    calibration,
    calibration_plot,
    in_notebook,
    save_figure,
)
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
    target_bins = np.digitize(y, np.quantile(y, np.linspace(0, 1, N_TARGET_BINS + 1)[1:-1]), right=True)
    splitter = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    return list(splitter.split(np.zeros(len(y)), target_bins, chrom))


def tune(x_df: pl.DataFrame, y: np.ndarray, chrom: np.ndarray) -> GridSearchCV:
    """Run a grid search over RF_GRID on chromosome_folds, and refit the best random forest on all given rows."""
    forest = RandomForestRegressor(random_state=SEED, n_jobs=1)
    search = GridSearchCV(forest, RF_GRID, scoring="r2", cv=chromosome_folds(y, chrom), n_jobs=-1)
    return search.fit(x_df, y)


def save_onnx(estimator: RandomForestRegressor, path: Path) -> None:
    """
    Write the random forest as one ONNX TreeEnsembleRegressor node (ai.onnx.ml opset 3). The input "input" has type
    double and shape [N, len(MODEL_INPUTS)], with the columns of MODEL_INPUTS in that order. By the ONNX specification,
    the output "variable" has type float and shape [N, 1]. The metadata key feature_names holds MODEL_INPUTS as a JSON
    list.

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
parser.add_argument(
    "--benchmark-dir", required=True, type=Path, help="--out-dir of build_benchmarks.py, with tcga.parquet"
)
parser.add_argument("--gff3", required=True, type=Path, help="GENCODE GFF3 of GRCh38, optionally gzip-compressed")
parser.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA with chr names and a .fai index")
parser.add_argument("--out-dir", required=True, type=Path, help="directory for the outputs")
# In a notebook, sys.argv holds the arguments of the kernel
args = parser.parse_args(shlex.split(os.environ.get("TRAIN_MODEL_ARGS", "")) if in_notebook() else None)
logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s")
for name in ("train_model", "nmd_scanner.benchmark"):
    logging.getLogger(name).setLevel(logging.INFO)
args.out_dir.mkdir(parents=True, exist_ok=True)

# %%
tcga_df = pl.read_parquet(args.benchmark_dir / "tcga.parquet")
features_df = run_nmd_scanner(tcga_df, args.out_dir / "variants.vcf", args.gff3, args.fasta, prediction=False)
rows_df = join_features(tcga_df, features_df, "TCGA").sort("source_row")
logger.info("%d TCGA rows, %d of them with nmd_efficiency_prediction_status ok", len(tcga_df), len(rows_df))
# A DataFrame, so that the random forest keeps the input names as feature_names_in_
x_df = rows_df.select(pl.col(MODEL_INPUTS).cast(pl.Float64))
y = rows_df["NMD_efficiency"].to_numpy()
chrom = rows_df["chrom"].to_numpy()

# %% [markdown]
# The histogram shows the target `NMD_efficiency` of the training rows. 0 means no NMD, and 1 means that the variant
# allele fraction in the RNA is half of that in the DNA.

# %%
plot = (
    pn.ggplot(rows_df, pn.aes("NMD_efficiency"))
    + pn.geom_histogram(bins=60, fill=BLUE, color=SURFACE, size=0.25)
    + pn.labs(title=f"NMD_efficiency of {len(y):,} TCGA rows", x="NMD_efficiency", y="Rows")
    + pn.theme(panel_grid_major_x=pn.element_blank())
)
save_figure(plot, args.out_dir, "target_histogram", 6, 3.5)

# %%
# The rows of each fold in a block, in fold order, as validate_model.py expects them
predictions = []
folds = []
for fold, (train, test) in enumerate(chromosome_folds(y, chrom)):
    search = tune(x_df[train], y[train], chrom[train])
    prediction = search.predict(x_df[test])
    columns = ["source_row", *KEY, "gene", "nmd_escape", "NMD_efficiency"]
    predictions.append(
        rows_df[test].select(columns).with_columns(fold=pl.lit(fold, dtype=pl.Int64), prediction=pl.Series(prediction))
    )
    chromosomes = sorted(set(chrom[test].tolist()))
    r2 = r2_score(y[test], prediction)
    folds.append({"chromosomes": chromosomes, "rows": len(test), "r2": r2, "params": search.best_params_})
    logger.info("fold %d (%s): R2 %.3f, %s", fold, " ".join(chromosomes), r2, search.best_params_)
oof_df = pl.concat(predictions)
oof_df.write_parquet(args.out_dir / "oof_predictions.parquet")

# %% [markdown]
# The scatter shows the measured NMD efficiency of each TCGA row over its out-of-fold prediction, with Spearman and R2
# over all rows. The gray line marks equal values.

# %%
oof_measured = oof_df["NMD_efficiency"].to_numpy()
oof_spearman = spearmanr(oof_df["prediction"].to_numpy(), oof_measured).statistic
oof_r2 = r2_score(oof_measured, oof_df["prediction"].to_numpy())
plot = (
    pn.ggplot(oof_df, pn.aes("prediction", "NMD_efficiency"))
    + pn.geom_point(color=BLUE, alpha=0.3, size=1.5, stroke=0)
    + pn.geom_abline(color=AXIS)
    + pn.labs(
        title=f"TCGA, out of fold, {len(oof_df):,} rows\nSpearman {oof_spearman:.2f}, R2 {oof_r2:.2f}",
        x="Predicted NMD efficiency",
        y="Measured NMD efficiency",
    )
)
save_figure(plot, args.out_dir, "tcga_measured_vs_predicted", 5, 4.5)

# %% [markdown]
# The dots show the mean measured over the mean predicted NMD efficiency in each quantile bin of the out-of-fold
# prediction, with the 95% CI of the measured mean. The blue line is the least-squares fit of measured on predicted over
# all rows, and the gray line marks equal values.

# %%
bins_df, slope, intercept = calibration(oof_df)
plot = calibration_plot(bins_df.with_columns(slope=slope, intercept=intercept)) + pn.labs(
    title=f"TCGA, out of fold, {len(oof_df):,} rows\nFit: slope {slope:.2f}, intercept {intercept:.2f}"
)
save_figure(plot, args.out_dir, "calibration_oof", 5, 6)

# %%
final = tune(x_df, y, chrom)
onnx_path = args.out_dir / "nmd_efficiency_rf.onnx"
save_onnx(final.best_estimator_, onnx_path)
session = onnxruntime.InferenceSession(onnx_path, providers=["CPUExecutionProvider"])
onnx_prediction = session.run(None, {"input": x_df.to_numpy()})[0].ravel()
if (onnx_prediction != final.predict(x_df).astype(np.float32)).any():
    raise ValueError(f"{onnx_path} predicts other values than the random forest")

results = {"rows": len(rows_df), "folds": folds, "inner_cv_r2": final.best_score_, "params": final.best_params_}
(args.out_dir / "cv_results.json").write_text(json.dumps(results, indent=2) + "\n")
joblib.dump(final.best_estimator_, args.out_dir / "nmd_efficiency_rf.pkl")
logger.info("final model: inner CV R2 %.3f, %s", final.best_score_, final.best_params_)

# %% [markdown]
# The bars show the impurity importance and the permutation importance of the model inputs in the final random forest.
# The permutation importance is the loss in R2 on the training rows when one input is permuted, with the standard
# deviation over the permutations.

# %%
# The random forest checks the input names. permutation_importance keeps them only in a pandas DataFrame.
permutation = permutation_importance(
    final.best_estimator_, x_df.to_pandas(), y, scoring="r2", n_repeats=N_PERMUTATIONS, random_state=SEED, n_jobs=-1
)
importance_df = pl.concat(
    [
        pl.DataFrame(
            {
                "kind": "Impurity importance\nMean decrease in impurity",
                "input": MODEL_INPUTS,
                "value": final.best_estimator_.feature_importances_,
            }
        ),
        pl.DataFrame(
            {
                "kind": (
                    "Permutation importance\nDecrease in R2 on the training rows,"
                    f"\nmean and SD of {N_PERMUTATIONS} permutations"
                ),
                "input": MODEL_INPUTS,
                "value": permutation.importances_mean,
                "sd": permutation.importances_std,
            }
        ),
    ],
    how="diagonal",
).sort("value", "input")
# Each panel orders its bars by importance, so a bar is one kind and input. Bars of equal importance are ordered
# by input name. The axis shows only the input.
bars = importance_df["kind"] + "|" + importance_df["input"]
importance_df = importance_df.with_columns(bar=bars.cast(pl.Enum(bars)))
plot = (
    pn.ggplot(importance_df, pn.aes("bar", "value"))
    + pn.geom_col(fill=BLUE, width=0.6)
    + pn.geom_linerange(
        pn.aes(ymin="value - sd", ymax="value + sd"), data=importance_df.drop_nulls(), color=SECONDARY_INK
    )
    + pn.scale_x_discrete(labels=lambda bars: [bar.split("|")[1] for bar in bars])
    + pn.coord_flip()
    + pn.facet_wrap("kind", scales="free")
    + pn.labs(x="", y="Importance")
    + pn.theme(panel_grid_major_y=pn.element_blank())
)
save_figure(plot, args.out_dir, "feature_importance", 9, 5)

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
    x_df,
    y,
    cv=learning_folds,
    train_sizes=learning_sizes,
    scoring="r2",
    shuffle=True,
    random_state=SEED,
    n_jobs=-1,
)
curves_df = pl.concat(
    [
        pl.DataFrame(
            {
                "rows": train_sizes,
                "mean": scores.mean(axis=1),
                "sd": scores.std(axis=1),
                "label": label,
                "color": color,
            }
        )
        for scores, label, color in ((train_scores, "Training rows", MUTED_INK), (test_scores, "Held-out fold", BLUE))
    ]
)
plot = (
    pn.ggplot(curves_df, pn.aes("rows", "mean"))
    + pn.geom_ribbon(pn.aes(ymin="mean - sd", ymax="mean + sd", fill="color"), alpha=0.1)
    + pn.geom_line(pn.aes(color="color"))
    + pn.geom_point(pn.aes(color="color"), size=2.5)
    + pn.geom_text(
        pn.aes(label="label"),
        data=curves_df.filter(pl.col("rows") == train_sizes[-1]),
        ha="left",
        nudge_x=0.03 * train_sizes[-1],
        color=SECONDARY_INK,
        size=9,
    )
    + pn.scale_color_identity()
    + pn.scale_fill_identity()
    + pn.scale_x_continuous(expand=(0.05, 0, 0.2, 0))
    + pn.labs(title="Learning curve of the final hyperparameters", x="Training rows", y="R2")
)
save_figure(plot, args.out_dir, "learning_curve", 6, 4)

# %% [markdown]
# The curve shows the out-of-bag error of the final hyperparameters, fitted on all rows with a growing number of trees.
# The out-of-bag prediction of a row is the mean prediction of the trees that did not draw it, and the out-of-bag error
# is 1 - R2 of these predictions.

# %%
oob_errors = [
    1 - clone(final.best_estimator_).set_params(n_estimators=trees, oob_score=True, n_jobs=-1).fit(x_df, y).oob_score_
    for trees in OOB_TREES
]
final_trees = final.best_params_["n_estimators"]
plot = (
    pn.ggplot(pl.DataFrame({"trees": OOB_TREES, "error": oob_errors}), pn.aes("trees", "error"))
    + pn.geom_vline(xintercept=final_trees, color=AXIS)
    + pn.annotate(
        "text",
        x=final_trees,
        y=np.inf,
        label=f"Final model, {final_trees} trees",
        ha="left",
        va="top",
        nudge_x=5,
        color=SECONDARY_INK,
        size=9,
    )
    + pn.geom_line(color=BLUE)
    + pn.geom_point(color=BLUE, size=2.5)
    + pn.labs(title="Out-of-bag error of the final hyperparameters", x="Trees", y="1 - out-of-bag R2")
)
save_figure(plot, args.out_dir, "oob_error", 6, 4)
