import hashlib
import logging

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest

import nmd_scanner
import nmd_scanner.model as model_module
from nmd_scanner.cli import annotate
from nmd_scanner.model import MODEL_PATH, add_nmd_pred_score
from nmd_scanner.schema import MODEL_INPUTS, empty_table, output_column_kinds, to_arrow

ARGS = ("resources/test_files/test_variants.vcf", "resources/chr18.gff3.gz", "resources/chr18.fa.gz")


@pytest.fixture(scope="module")
def results():
    """The annotate result of the chr18 test variants, without the score. Some rows have nmd_model_status ok."""
    return annotate(*ARGS, score=False)


def test_the_packaged_model_is_the_trained_model():
    sha256 = hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest()
    assert sha256 == "e19f3c2b6c55c70b6452f657fea7e06a98ac33295b6128e0de0f2ef8fbd700b8"


def test_add_nmd_pred_score_is_exported_from_the_package():
    assert nmd_scanner.add_nmd_pred_score is add_nmd_pred_score


def test_the_score_is_finite_exactly_in_the_rows_with_nmd_model_status_ok(results):
    scored = add_nmd_pred_score(results)

    ok = (results["nmd_model_status"] == "ok").to_numpy()
    assert ok.any() and not ok.all()
    score = scored["nmd_pred_score"]
    assert score.dtype == np.float64
    assert np.isfinite(score[ok]).all()
    assert score[~ok].isna().all()


def test_the_score_follows_nmd_model_status_and_the_other_columns_stay_unchanged(results):
    original = results.copy()
    with_key = results.assign(my_key="sample_1")

    scored = add_nmd_pred_score(with_key)

    assert list(scored.columns) == [*output_column_kinds(), "my_key"]
    pd.testing.assert_frame_equal(scored.drop(columns="nmd_pred_score"), with_key)
    pd.testing.assert_frame_equal(results, original)


@pytest.mark.parametrize("has_rows", [False, True])
def test_without_a_row_to_score_the_score_is_float64_and_the_model_does_not_run(monkeypatch, results, has_rows):
    import onnxruntime

    def run(*args, **kwargs):
        raise AssertionError("the model ran")

    monkeypatch.setattr(onnxruntime.InferenceSession, "run", run)
    table = results[results["nmd_model_status"] != "ok"] if has_rows else empty_table()

    scored = add_nmd_pred_score(table)

    assert (len(scored) > 0) == has_rows
    assert list(scored.columns) == list(output_column_kinds())
    assert scored["nmd_pred_score"].dtype == np.float64
    assert scored["nmd_pred_score"].isna().all()
    assert to_arrow(scored).schema.field("nmd_pred_score").type == pa.float64()


def test_a_model_whose_feature_names_differ_from_model_inputs_raises(monkeypatch, results):
    monkeypatch.setattr(model_module, "MODEL_INPUTS", list(reversed(MODEL_INPUTS)))

    with pytest.raises(ValueError, match="differ from MODEL_INPUTS"):
        add_nmd_pred_score(results)


def test_annotate_adds_the_score_unless_score_is_false(results):
    pd.testing.assert_frame_equal(annotate(*ARGS), add_nmd_pred_score(results))


def test_add_nmd_pred_score_logs_the_rows_per_nmd_model_status(results, caplog):
    with caplog.at_level(logging.INFO, logger="nmd_scanner.model"):
        add_nmd_pred_score(results)

    counts = results["nmd_model_status"].value_counts()
    for status, count in counts[counts > 0].items():
        assert f"{count} rows with nmd_model_status {status}" in caplog.text
