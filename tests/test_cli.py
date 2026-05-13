import pandas as pd
import pytest

from nmd_scanner.cli import is_valid_output_path, write_results


def test_is_valid_output_path_accepts_csv_in_existing_dir(tmp_path):
    assert is_valid_output_path(str(tmp_path / "results.csv")) is True


def test_is_valid_output_path_accepts_parquet_extensions(tmp_path):
    assert is_valid_output_path(str(tmp_path / "results.parquet")) is True
    assert is_valid_output_path(str(tmp_path / "results.pq")) is True


def test_is_valid_output_path_accepts_bare_filename_in_cwd():
    assert is_valid_output_path("results.csv") is True


def test_is_valid_output_path_rejects_existing_directory(tmp_path):
    assert is_valid_output_path(str(tmp_path)) is False


def test_is_valid_output_path_rejects_missing_parent_directory(tmp_path):
    assert is_valid_output_path(str(tmp_path / "missing" / "results.csv")) is False


def test_is_valid_output_path_rejects_unsupported_extension(tmp_path):
    assert is_valid_output_path(str(tmp_path / "results.tsv")) is False
    assert is_valid_output_path(str(tmp_path / "results")) is False


def test_write_results_csv_roundtrip(tmp_path):
    df = pd.DataFrame({"transcript_id": ["t1", "t2"], "nmd_escape": [True, False]})
    out = tmp_path / "results.csv"

    write_results(df, str(out))

    assert out.exists()
    loaded = pd.read_csv(out)
    assert list(loaded.columns) == ["transcript_id", "nmd_escape"]
    assert list(loaded["transcript_id"]) == ["t1", "t2"]


def test_write_results_rejects_unsupported_extension(tmp_path):
    df = pd.DataFrame({"x": [1]})
    with pytest.raises(ValueError, match="Unsupported output extension"):
        write_results(df, str(tmp_path / "results.tsv"))


def test_write_results_parquet_roundtrip(tmp_path):
    pytest.importorskip("pyarrow")
    df = pd.DataFrame({"transcript_id": ["t1", "t2"], "nmd_escape": [True, False]})
    out = tmp_path / "results.parquet"

    write_results(df, str(out))

    assert out.exists()
    loaded = pd.read_parquet(out)
    assert list(loaded.columns) == ["transcript_id", "nmd_escape"]
    assert list(loaded["transcript_id"]) == ["t1", "t2"]
