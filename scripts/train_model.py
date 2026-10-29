"""
Train the NMD efficiency model nmd_efficiency_rf.onnx, a random forest on the NMD-Scanner features of the TCGA benchmark
(see nmd_scanner.benchmark.datasets.tcga), and predict each training row out of fold. A grid search on folds of whole
chromosomes tunes the random forest (see chromosome_folds). scripts/validate_model.py.py draws the figures.

Run it from the repo root with:

\b
    uv run --extra scripts python scripts/train_model.py

--benchmark-dir is the --out-dir of scripts/build_benchmarks.py, with the same GFF3 and FASTA.

\b
Outputs in --out-dir:
- nmd_efficiency_rf.onnx (see save_onnx), and the same random forest as the joblib pickle nmd_efficiency_rf.pkl
- cv_results.json: the R2 of each fold and the hyperparameters
- oof_predictions.parquet: the out-of-fold prediction of each training row
- feature_importance.parquet, learning_curve.parquet and oob_error.parquet: the numbers of the training figures
"""

import json
import logging
from pathlib import Path

import click
import joblib
import numpy as np
import onnx
import onnxruntime
import polars as pl
from onnx import TensorProto, helper, numpy_helper
from pyfaidx import Fasta
from skl2onnx.common.tree_ensemble import add_tree_to_attribute_pairs, get_default_tree_regressor_attribute_pairs
from sklearn.base import clone
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, learning_curve

from nmd_scanner import read_annotation
from nmd_scanner.benchmark import KEY, join_features, run_nmd_scanner
from nmd_scanner.schema import MODEL_INPUTS

logger = logging.getLogger("train_model")

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
    """
    Split the rows into N_FOLDS folds of whole chromosomes, stratified by N_TARGET_BINS quantile bins of y. The folds
    then have similar sizes and target distributions. Genes can overlap (antisense, nested), so a grouping by gene could
    still leak.
    """
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


@click.command(help=__doc__)
@click.option(
    "--gff3",
    default="data/gencode.v42.annotation.gff3.gz",
    show_default=True,
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    help="GENCODE GFF3 of GRCh38, optionally gzip-compressed",
)
@click.option(
    "--fasta",
    default="data/GRCh38.fa",
    show_default=True,
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    help="GRCh38 FASTA with chr names and a .fai index",
)
@click.option(
    "--benchmark-dir",
    default="data/benchmarks",
    show_default=True,
    type=click.Path(path_type=Path, exists=True, file_okay=False),
    help="--out-dir of build_benchmarks.py, with tcga.parquet",
)
@click.option(
    "--out-dir",
    default="data/train",
    show_default=True,
    type=click.Path(path_type=Path, file_okay=False),
    help="directory for the outputs",
)
def main(gff3: Path, fasta: Path, benchmark_dir: Path, out_dir: Path) -> None:
    logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    for name in ("train_model", "nmd_scanner.benchmark"):
        logging.getLogger(name).setLevel(logging.INFO)
    out_dir.mkdir(parents=True, exist_ok=True)

    tcga_df = pl.read_parquet(benchmark_dir / "tcga.parquet")
    genome = Fasta(str(fasta))
    features_df = run_nmd_scanner(tcga_df, read_annotation(gff3, genome), genome, prediction=False)
    rows_df = join_features(tcga_df, features_df, "TCGA").sort("source_row")
    logger.info("%d TCGA rows, %d of them with nmd_efficiency_prediction_status ok", len(tcga_df), len(rows_df))
    # A DataFrame, so that the random forest keeps the input names as feature_names_in_
    x_df = rows_df.select(pl.col(MODEL_INPUTS).cast(pl.Float64))
    y = rows_df["NMD_efficiency"].to_numpy()
    chrom = rows_df["chrom"].to_numpy()

    # The rows of each fold in a block, in fold order, as validate_model.py.py expects them
    predictions = []
    folds = []
    for fold, (train, test) in enumerate(chromosome_folds(y, chrom)):
        search = tune(x_df[train], y[train], chrom[train])
        prediction = search.predict(x_df[test])
        columns = ["source_row", *KEY, "gene", "nmd_escape", "NMD_efficiency"]
        predictions.append(
            rows_df[test]
            .select(columns)
            .with_columns(fold=pl.lit(fold, dtype=pl.Int64), prediction=pl.Series(prediction))
        )
        chromosomes = sorted(set(chrom[test].tolist()))
        r2 = r2_score(y[test], prediction)
        folds.append({"chromosomes": chromosomes, "rows": len(test), "r2": r2, "params": search.best_params_})
        logger.info("fold %d (%s): R2 %.3f, %s", fold, " ".join(chromosomes), r2, search.best_params_)
    pl.concat(predictions).write_parquet(out_dir / "oof_predictions.parquet")

    final = tune(x_df, y, chrom)
    onnx_path = out_dir / "nmd_efficiency_rf.onnx"
    save_onnx(final.best_estimator_, onnx_path)
    session = onnxruntime.InferenceSession(onnx_path, providers=["CPUExecutionProvider"])
    onnx_prediction = session.run(None, {"input": x_df.to_numpy()})[0].ravel()
    if (onnx_prediction != final.predict(x_df).astype(np.float32)).any():
        raise ValueError(f"{onnx_path} predicts other values than the random forest")

    results = {"rows": len(rows_df), "folds": folds, "inner_cv_r2": final.best_score_, "params": final.best_params_}
    (out_dir / "cv_results.json").write_text(json.dumps(results, indent=2) + "\n")
    joblib.dump(final.best_estimator_, out_dir / "nmd_efficiency_rf.pkl")
    logger.info("final model: inner CV R2 %.3f, %s", final.best_score_, final.best_params_)

    # The impurity importance and the permutation importance of each input in the final random forest. The random
    # forest checks the input names. permutation_importance keeps them only in a pandas DataFrame.
    permutation = permutation_importance(
        final.best_estimator_, x_df.to_pandas(), y, scoring="r2", n_repeats=N_PERMUTATIONS, random_state=SEED, n_jobs=-1
    )
    pl.concat(
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
    ).write_parquet(out_dir / "feature_importance.parquet")

    # The learning curve of the final hyperparameters, in the chromosome folds
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
    # One row per training size and fold
    pl.DataFrame(
        {
            "training_rows": np.repeat(train_sizes, N_FOLDS),
            "fold": np.tile(np.arange(N_FOLDS), len(train_sizes)),
            "training_r2": train_scores.ravel(),
            "held_out_r2": test_scores.ravel(),
        }
    ).write_parquet(out_dir / "learning_curve.parquet")

    # The out-of-bag R2 of the final hyperparameters, fitted on all rows with a growing number of trees
    oob_r2 = [
        clone(final.best_estimator_).set_params(n_estimators=trees, oob_score=True, n_jobs=-1).fit(x_df, y).oob_score_
        for trees in OOB_TREES
    ]
    pl.DataFrame({"trees": OOB_TREES, "oob_r2": oob_r2}).write_parquet(out_dir / "oob_error.parquet")


if __name__ == "__main__":
    main()
