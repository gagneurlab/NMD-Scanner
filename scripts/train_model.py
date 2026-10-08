# /// script
# requires-python = "==3.12.*"
# dependencies = [
#     "nmd-scanner==0.4.0",
#     "scikit-learn==1.9.1",
#     "pandas==3.0.6",
#     "numpy==2.5.3",
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
Train the NMD efficiency model nmd_efficiency_rf.onnx on NMD-Scanner features of the NMDEff TCGA benchmark, and
predict each training row out of fold.

Run it with:

    uv run scripts/train_model.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa --out-dir out/

Steps:
1. Download tcga_dataset.csv of NMDEff and add the genomic alleles, with the steps of make_benchmark_vcfs.py. Its
   column NMD_efficiency is the target.
2. Run nmd_scanner.annotate on the unique variants. Join each TCGA row to the NMD-Scanner row of its transcript and
   variant, and keep the rows with nmd_model_status "ok". The inputs are nmd_scanner.schema.MODEL_INPUTS as float64.
3. Split the rows into 5 folds with StratifiedGroupKFold over the chromosomes, stratified by 5 quantile bins of the
   target. The folds then have similar sizes and target distributions. Genes can overlap (antisense, nested), so a
   grouping by gene could still leak.
4. For each fold, tune a random forest on the other folds, and predict the rows of the fold. The grid search of the
   tuning splits its rows the same way. validate_model.py reads these out-of-fold predictions.
5. Tune the random forest the same way on all rows, and save it as ONNX (see save_onnx). Fail if the ONNX model
   predicts other values than the random forest, rounded to float32.

Outputs in --out-dir: variants.vcf (the input of nmd_scanner.annotate), oof_predictions.parquet,
nmd_efficiency_rf.onnx and cv_results.json (the R2 of each fold and the hyperparameters).
"""

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import onnx
import onnxruntime
import pandas as pd
from make_benchmark_vcfs import KEY, add_alleles, add_key, read_nmdeff, read_strands, write_vcf
from onnx import TensorProto, helper, numpy_helper
from pyfaidx import Fasta
from skl2onnx.common.tree_ensemble import add_tree_to_attribute_pairs, get_default_tree_regressor_attribute_pairs
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold

import nmd_scanner
from nmd_scanner.schema import MODEL_INPUTS

logger = logging.getLogger("train_model")

SEED = 42
N_FOLDS = 5
# The folds are stratified by this many quantile bins of the target
N_TARGET_BINS = 5
# refined_grid5 of scripts/train_new.ipynb, the grid that found the hyperparameters of best_model.pkl: 324 candidates
RF_GRID = {
    "n_estimators": [200, 300, 400],
    "max_depth": [5, 7, 9, 11],
    "max_features": ["sqrt", 0.2, 0.3],
    "min_samples_split": [5, 10, 15],
    "min_samples_leaf": [2, 3, 5],
}


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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, type=Path, help="GENCODE GFF3 of GRCh38, optionally gzip-compressed")
    parser.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA with chr names and a .fai index")
    parser.add_argument("--out-dir", required=True, type=Path, help="directory for the outputs")
    args = parser.parse_args()
    logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    logger.setLevel(logging.INFO)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tcga = add_alleles(read_nmdeff("tcga_dataset.csv"), read_strands(args.gff3), Fasta(str(args.fasta)))
    vcf = args.out_dir / "variants.vcf"
    write_vcf(tcga.drop_duplicates(["chrom", "pos", "ref", "alt"]), vcf)
    features = add_key(nmd_scanner.annotate(vcf, args.gff3, args.fasta, sequences=False))
    # tcga_row is the index of the row in tcga_dataset.csv
    tcga = tcga[[*KEY, "NMD_efficiency"]].assign(tcga_row=np.arange(len(tcga)))
    rows = tcga.merge(features, on=KEY, how="inner", validate="many_to_one")
    rows = rows[rows["nmd_model_status"] == "ok"].reset_index(drop=True)
    logger.info("%d TCGA rows, %d of them with nmd_model_status ok", len(tcga), len(rows))
    x = rows[MODEL_INPUTS].astype("float64")
    y = rows["NMD_efficiency"].to_numpy(dtype="float64")
    chrom = rows["chrom"].to_numpy(dtype=str)

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
    pd.concat(predictions).to_parquet(args.out_dir / "oof_predictions.parquet", index=False)

    final = tune(x, y, chrom)
    logger.info("final model: inner CV R2 %.3f, %s", final.best_score_, final.best_params_)
    onnx_path = args.out_dir / "nmd_efficiency_rf.onnx"
    save_onnx(final.best_estimator_, onnx_path)
    session = onnxruntime.InferenceSession(onnx_path, providers=["CPUExecutionProvider"])
    onnx_prediction = session.run(None, {"input": x.to_numpy(dtype=np.float64)})[0].ravel()
    if (onnx_prediction != final.predict(x).astype(np.float32)).any():
        raise ValueError(f"{onnx_path} predicts other values than the random forest")

    results = {"rows": len(rows), "folds": folds, "inner_cv_r2": final.best_score_, "params": final.best_params_}
    (args.out_dir / "cv_results.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
