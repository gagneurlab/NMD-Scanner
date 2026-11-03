"""
Train the NMD efficiency model nmd_efficiency_rf.onnx, a random forest on the NMD-Scanner features of the TCGA benchmark
(see nmd_scanner.benchmark.datasets.tcga), and predict each TCGA variant-sample pair out of fold. A grid search on folds
of whole chromosomes tunes the random forest (see chromosome_folds). scripts/validate_model.py.py draws the figures.

Run it from the repo root with:

\b
    uv run --extra scripts python scripts/train_model.py

--benchmark-dir is the --out-dir of scripts/build_benchmarks.py, with the same GFF3 and FASTA.

\b
Outputs in --out-dir:
- nmd_efficiency_rf.onnx (see save_onnx), and the same random forest as the joblib pickle nmd_efficiency_rf.pkl
- cv_results.json: the R2 of each fold, the hyperparameters of the grid search, N_ESTIMATORS and N_PERMUTATIONS
- oof_predictions.parquet: the out-of-fold prediction of each TCGA variant-sample pair
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
from onnx import helper
from pyfaidx import Fasta
from skl2onnx import convert_sklearn
from skl2onnx.common.data_types import FloatTensorType
from sklearn.base import clone
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, learning_curve

from nmd_scanner import annotate_variants, read_annotation
from nmd_scanner.benchmark import DATASETS, KEY, VARIANT_COLUMNS, join_features
from nmd_scanner.schema import MODEL_INPUTS

logger = logging.getLogger("train_model")

SEED = 42
N_FOLDS = 5
# The folds are stratified by this many quantile bins of the target
N_TARGET_BINS = 5
# Each random forest has this many trees. The out-of-bag R2 levels off at 200 trees, as the oob_error figure shows.
N_ESTIMATORS = 200
# The grid of the random forest search: 108 candidates
RF_GRID = {
    "max_depth": [5, 7, 9, 11],
    "max_features": ["sqrt", 0.2, 0.3],
    "min_samples_split": [5, 10, 15],
    "min_samples_leaf": [2, 3, 5],
}
# The permutation importance permutes each input this many times
N_PERMUTATIONS = 10
# In every fold, the learning curve fits the final hyperparameters on these shares. A share is of the training
# variant-sample pairs of the fold with the fewest.
LEARNING_CURVE_SIZES = np.linspace(0.1, 1.0, 8)
# The out-of-bag error curve fits the final hyperparameters with these numbers of trees
OOB_TREES = [50, 100, 200, 300, 400, 500]
# The ONNX model holds and sums the leaf values in float32, so its predictions may differ from the random forest by up
# to this much
ONNX_TOLERANCE = 1e-6


def chromosome_folds(y: np.ndarray, chrom: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    """
    Split the variant-sample pairs into N_FOLDS folds of whole chromosomes, stratified by N_TARGET_BINS quantile bins of
    y. The folds then have similar sizes and target distributions. Genes can overlap (antisense, nested), so a grouping
    by gene could still leak.
    """
    target_bins = np.digitize(y, np.quantile(y, np.linspace(0, 1, N_TARGET_BINS + 1)[1:-1]), right=True)
    splitter = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    return list(splitter.split(np.zeros(len(y)), target_bins, chrom))


def tune(x_df: pl.DataFrame, y: np.ndarray, chrom: np.ndarray) -> GridSearchCV:
    """Run a grid search over RF_GRID on chromosome_folds, and refit the best random forest on all of x_df and y."""
    forest = RandomForestRegressor(n_estimators=N_ESTIMATORS, random_state=SEED, n_jobs=1)
    search = GridSearchCV(forest, RF_GRID, scoring="r2", cv=chromosome_folds(y, chrom), n_jobs=-1)
    return search.fit(x_df, y)


def save_onnx(estimator: RandomForestRegressor, path: Path) -> None:
    """
    Write the random forest as one ONNX TreeEnsembleRegressor node, with the skl2onnx converter. The input "input" has
    type float and shape [N, len(MODEL_INPUTS)], with the columns of MODEL_INPUTS in that order. The output "variable"
    has type float and shape [N, 1]. The metadata key feature_names holds MODEL_INPUTS as a JSON list.
    """
    # Without a name, skl2onnx names the graph at random, so each export of the same random forest gives another file
    model = convert_sklearn(
        estimator,
        name="nmd_efficiency_rf",
        initial_types=[("input", FloatTensorType([None, len(MODEL_INPUTS)]))],
        target_opset=15,
    )
    helper.set_model_props(model, {"feature_names": json.dumps(list(MODEL_INPUTS))})
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
    results_df = annotate_variants(
        tcga_df.select(VARIANT_COLUMNS).unique().to_pandas(),
        read_annotation(gff3, genome),
        genome,
        sequences=False,
        prediction=False,
    )
    training_df = join_features(tcga_df, results_df, DATASETS["tcga"].NAME).sort("source_row")
    logger.info(
        "%d TCGA variant-sample pairs, %d of them with nmd_efficiency_prediction_status ok",
        len(tcga_df),
        len(training_df),
    )
    # A DataFrame, so that the random forest keeps the input names as feature_names_in_
    x_df = training_df.select(pl.col(MODEL_INPUTS).cast(pl.Float64))
    y = training_df["NMD_efficiency"].to_numpy()
    chrom = training_df["chrom"].to_numpy()
    folds = chromosome_folds(y, chrom)

    # The rows of each fold in a block, in fold order
    predictions = []
    fold_results = []
    for fold, (train, test) in enumerate(folds):
        search = tune(x_df[train], y[train], chrom[train])
        prediction = search.predict(x_df[test])
        columns = ["source_row", *KEY, "gene", "nmd_escape", "NMD_efficiency"]
        predictions.append(
            training_df[test]
            .select(columns)
            .with_columns(fold=pl.lit(fold, dtype=pl.Int64), prediction=pl.Series(prediction))
        )
        chromosomes = sorted(set(chrom[test].tolist()))
        r2 = r2_score(y[test], prediction)
        fold_results.append(
            {"chromosomes": chromosomes, "variant_sample_pairs": len(test), "r2": r2, "params": search.best_params_}
        )
        logger.info("fold %d (%s): R2 %.3f, %s", fold, " ".join(chromosomes), r2, search.best_params_)
    pl.concat(predictions).write_parquet(out_dir / "oof_predictions.parquet")

    final = tune(x_df, y, chrom)
    onnx_path = out_dir / "nmd_efficiency_rf.onnx"
    save_onnx(final.best_estimator_, onnx_path)
    session = onnxruntime.InferenceSession(onnx_path, providers=["CPUExecutionProvider"])
    onnx_prediction = session.run(None, {"input": x_df.cast(pl.Float32).to_numpy()})[0].ravel()
    difference = np.abs(onnx_prediction - final.predict(x_df)).max()
    if difference > ONNX_TOLERANCE:
        raise ValueError(f"{onnx_path} predicts values up to {difference:.2g} away from the random forest")

    results = {
        "variant_sample_pairs": len(training_df),
        "folds": fold_results,
        "inner_cv_r2": final.best_score_,
        "params": final.best_params_,
        "n_estimators": N_ESTIMATORS,
        "n_permutations": N_PERMUTATIONS,
    }
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
                    "kind": "impurity",
                    "input": MODEL_INPUTS,
                    "value": final.best_estimator_.feature_importances_,
                }
            ),
            pl.DataFrame(
                {
                    "kind": "permutation",
                    "input": MODEL_INPUTS,
                    "value": permutation.importances_mean,
                    "sd": permutation.importances_std,
                }
            ),
        ],
        how="diagonal",
    ).write_parquet(out_dir / "feature_importance.parquet")

    # The learning curve of the final hyperparameters, in the chromosome folds. learning_curve would scale shares by the
    # training variant-sample pairs of the first fold only.
    learning_sizes = (LEARNING_CURVE_SIZES * min(len(train) for train, _ in folds)).astype(int)
    # A copy of the final random forest with its hyperparameters, unfitted. The final model stays as it is.
    train_sizes, train_scores, test_scores = learning_curve(
        clone(final.best_estimator_),
        x_df,
        y,
        cv=folds,
        train_sizes=learning_sizes,
        scoring="r2",
        shuffle=True,
        random_state=SEED,
        n_jobs=-1,
    )
    # One row per training size and fold
    pl.DataFrame(
        {
            "training_variant_sample_pairs": np.repeat(train_sizes, N_FOLDS),
            "fold": np.tile(np.arange(N_FOLDS), len(train_sizes)),
            "training_r2": train_scores.ravel(),
            "held_out_r2": test_scores.ravel(),
        }
    ).write_parquet(out_dir / "learning_curve.parquet")

    # The out-of-bag R2 of the final hyperparameters, fitted on all variant-sample pairs with a growing number of trees
    oob_r2 = [
        clone(final.best_estimator_).set_params(n_estimators=trees, oob_score=True, n_jobs=-1).fit(x_df, y).oob_score_
        for trees in OOB_TREES
    ]
    pl.DataFrame({"trees": OOB_TREES, "oob_r2": oob_r2}).write_parquet(out_dir / "oob_error.parquet")


if __name__ == "__main__":
    main()
