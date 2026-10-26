import hashlib
import logging

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest

import nmd_scanner
import nmd_scanner.model as model_module
from nmd_scanner.cli import annotate
from nmd_scanner.model import MODEL_PATH, add_nmd_efficiency_prediction
from nmd_scanner.schema import MODEL_INPUTS, empty_table, output_column_kinds, to_arrow

ARGS = ("resources/test_files/test_variants.vcf", "resources/chr18.gff3.gz", "resources/chr18.fa.gz")


@pytest.fixture(scope="module")
def results():
    """The annotate result of the chr18 test variants, without the prediction. Some rows have nmd_efficiency_prediction_status ok."""
    return annotate(*ARGS, prediction=False)


def test_the_packaged_model_is_the_trained_model():
    sha256 = hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest()
    assert sha256 == "e19f3c2b6c55c70b6452f657fea7e06a98ac33295b6128e0de0f2ef8fbd700b8"


def test_add_nmd_efficiency_prediction_is_exported_from_the_package():
    assert nmd_scanner.add_nmd_efficiency_prediction is add_nmd_efficiency_prediction


def test_the_prediction_is_finite_exactly_in_the_rows_with_nmd_efficiency_prediction_status_ok(results):
    predicted = add_nmd_efficiency_prediction(results)

    ok = (results["nmd_efficiency_prediction_status"] == "ok").to_numpy()
    assert ok.any() and not ok.all()
    prediction = predicted["nmd_efficiency_prediction"]
    assert prediction.dtype == np.float64
    assert np.isfinite(prediction[ok]).all()
    assert prediction[~ok].isna().all()


def test_the_prediction_follows_nmd_efficiency_prediction_status_and_the_other_columns_stay_unchanged(results):
    original = results.copy()
    with_key = results.assign(my_key="sample_1")

    predicted = add_nmd_efficiency_prediction(with_key)

    assert list(predicted.columns) == [*output_column_kinds(), "my_key"]
    pd.testing.assert_frame_equal(predicted.drop(columns="nmd_efficiency_prediction"), with_key)
    pd.testing.assert_frame_equal(results, original)


@pytest.mark.parametrize("has_rows", [False, True])
def test_without_a_row_to_prediction_the_prediction_is_float64_and_the_model_does_not_run(
    monkeypatch, results, has_rows
):
    import onnxruntime

    def run(*args, **kwargs):
        raise AssertionError("the model ran")

    monkeypatch.setattr(onnxruntime.InferenceSession, "run", run)
    table = results[results["nmd_efficiency_prediction_status"] != "ok"] if has_rows else empty_table()

    predicted = add_nmd_efficiency_prediction(table)

    assert (len(predicted) > 0) == has_rows
    assert list(predicted.columns) == list(output_column_kinds())
    assert predicted["nmd_efficiency_prediction"].dtype == np.float64
    assert predicted["nmd_efficiency_prediction"].isna().all()
    assert to_arrow(predicted).schema.field("nmd_efficiency_prediction").type == pa.float64()


def test_a_model_whose_feature_names_differ_from_model_inputs_raises(monkeypatch, results):
    monkeypatch.setattr(model_module, "MODEL_INPUTS", list(reversed(MODEL_INPUTS)))

    with pytest.raises(ValueError, match="differ from MODEL_INPUTS"):
        add_nmd_efficiency_prediction(results)


def test_annotate_adds_the_prediction_unless_prediction_is_false(results):
    pd.testing.assert_frame_equal(annotate(*ARGS), add_nmd_efficiency_prediction(results))


def test_add_nmd_efficiency_prediction_logs_the_rows_per_nmd_efficiency_prediction_status(results, caplog):
    with caplog.at_level(logging.INFO, logger="nmd_scanner.model"):
        add_nmd_efficiency_prediction(results)

    counts = results["nmd_efficiency_prediction_status"].value_counts()
    for status, count in counts[counts > 0].items():
        assert f"{count} rows with nmd_efficiency_prediction_status {status}" in caplog.text
