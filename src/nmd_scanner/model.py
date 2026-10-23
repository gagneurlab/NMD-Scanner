"""
The NMD efficiency model: a random forest that predicts the NMD efficiency of a PTC from the columns of MODEL_INPUTS.

scripts/train_model.py trains it and saves it as ONNX. The section "NMD efficiency score" of "Technical Notes.md"
describes the ONNX file.
"""

import importlib.resources
import json
import logging

import numpy as np
import onnxruntime
import pandas as pd

from nmd_scanner.schema import MODEL_INPUTS, MODEL_STATUSES

# The ONNX file of the model, inside the package. In a package on the file system, it is a pathlib.Path.
MODEL_PATH = importlib.resources.files("nmd_scanner") / "nmd_efficiency_rf.onnx"

logger = logging.getLogger(__name__)


def add_nmd_pred_score(results: pd.DataFrame) -> pd.DataFrame:
    """
    Return a copy of ``results`` with the column nmd_pred_score after nmd_model_status. It holds the NMD efficiency
    that the model predicts for the row.

    The model scores only the rows with nmd_model_status "ok". In the other rows, nmd_pred_score is NaN. The column
    has dtype float64, also for zero rows. ``results`` is not changed. ``annotate`` calls this function unless
    ``score=False``.

    :param results: DataFrame with nmd_model_status and the columns of MODEL_INPUTS, e.g. from
        ``add_features_and_rules``. Other columns stay as they are.
    :return: DataFrame with the columns of ``results`` and nmd_pred_score after nmd_model_status
    :raises ValueError: if the metadata key feature_names of the model differs from MODEL_INPUTS
    """

    session = onnxruntime.InferenceSession(MODEL_PATH.read_bytes(), providers=["CPUExecutionProvider"])
    feature_names = json.loads(session.get_modelmeta().custom_metadata_map["feature_names"])
    if feature_names != MODEL_INPUTS:
        raise ValueError(f"The model takes the inputs {feature_names}, which differ from MODEL_INPUTS {MODEL_INPUTS}.")

    status = results["nmd_model_status"]
    counts = status.value_counts()
    for value in MODEL_STATUSES:
        if counts.get(value, 0) > 0:
            logger.info("%d rows with nmd_model_status %s", counts[value], value)

    ok = (status == "ok").to_numpy(dtype=bool, na_value=False)
    score = np.full(len(results), np.nan)
    # Without a row to score, the session does not run
    if ok.any():
        # The bool inputs become 0.0 and 1.0. The output is float32, with shape [N, 1].
        inputs = results.loc[ok, MODEL_INPUTS].to_numpy(dtype=np.float64)
        score[ok] = session.run(None, {session.get_inputs()[0].name: inputs})[0].ravel()

    results = results.copy()
    results.insert(list(results.columns).index("nmd_model_status") + 1, "nmd_pred_score", score)
    return results
