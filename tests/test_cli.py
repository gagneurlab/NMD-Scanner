import pandas as pd
import pytest

from nmd_scanner.cli import is_valid_output_path, main, write_results


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


def test_main_end_to_end_smoke(tmp_path):
    """Run the full pipeline on the bundled chr18 test data and check the output file."""

    out = tmp_path / "smoke_results.csv"
    results = main(
        vcf_path="resources/test_files/test_variants.vcf",
        gtf_path="resources/chr18.gtf.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(out),
    )

    # Pipeline returns a non-empty DataFrame
    assert isinstance(results, pd.DataFrame)
    assert not results.empty

    # Output file written and reloadable
    assert out.exists()
    loaded = pd.read_csv(out)
    assert len(loaded) == len(results)

    # Columns from each pipeline stage are present
    for col in [
        "transcript_id",
        "variant_id",
        "ref_cds_seq",
        "alt_cds_seq",
        "utr3_length",
        "total_exon_count",
        "nmd_escape",
    ]:
        assert col in results.columns, f"missing column: {col}"
