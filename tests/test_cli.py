import pandas as pd
import pytest

from nmd_scanner.cli import is_valid_output_path, main, to_parquet_safe, write_results


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


def test_main_end_to_end_parquet_typed_columns(tmp_path):
    """
    Run the full pipeline on the bundled chr18 test data, write Parquet, and read it back.

    ``ref_all_stop_codons`` and ``alt_all_stop_codons`` hold (position, codon) tuples,
    e.g. (5442, "TGA"); pyarrow cannot infer a single type for a tuple mixing int and str,
    so they need a typed struct schema instead. ``transcript_exon_info`` holds
    (exon_number, exon_length) tuples; exon_number used to come from the GTF as a string in
    this column but as an int everywhere else, which pyarrow also rejects.
    """

    pytest.importorskip("pyarrow")
    import pyarrow as pa
    import pyarrow.parquet as pq

    out = tmp_path / "typed_results.parquet"
    results = main(
        vcf_path="resources/test_files/test_variants.vcf",
        gtf_path="resources/chr18.gtf.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(out),
    )

    assert out.exists()
    schema = pq.read_schema(str(out))

    stop_codon_type = pa.list_(pa.struct([pa.field("position", pa.int64()), pa.field("codon", pa.string())]))
    for column in ["ref_all_stop_codons", "alt_all_stop_codons"]:
        assert schema.field(column).type.equals(stop_codon_type), (
            f"{column} has unexpected parquet type: {schema.field(column).type}"
        )

    loaded = pd.read_parquet(out)
    assert len(loaded) == len(results)

    # Values are preserved, just reshaped from (position, codon) tuples to records
    has_stop_codons = results["ref_all_stop_codons"].apply(lambda v: isinstance(v, list) and len(v) > 0)
    sample_pos = results.index[has_stop_codons][0]
    expected = [{"position": pos, "codon": codon} for pos, codon in results.loc[sample_pos, "ref_all_stop_codons"]]
    assert list(loaded.loc[sample_pos, "ref_all_stop_codons"]) == expected

    # Exon numbers in transcript_exon_info are one type (int), not a mix of int and str
    exon_info_samples = loaded["transcript_exon_info"].dropna()
    exon_info_samples = exon_info_samples[exon_info_samples.apply(len) > 0]
    assert not exon_info_samples.empty
    for exon_number, exon_length in exon_info_samples.iloc[0]:
        assert not isinstance(exon_number, str)
        assert not isinstance(exon_length, str)


def test_write_results_parquet_types_stop_codon_columns(tmp_path):
    """
    ``to_parquet_safe`` (used by ``write_results``) turns (position, codon) tuples into
    {"position": ..., "codon": ...} records for every stop-codon column, including
    ``transcript_all_stop_codons`` and rows holding None, so Parquet gets a typed struct
    schema instead of raising ArrowInvalid.
    """

    pytest.importorskip("pyarrow")
    import pyarrow as pa
    import pyarrow.parquet as pq

    df = pd.DataFrame(
        {
            "transcript_id": ["t1", "t2"],
            "ref_all_stop_codons": [[(5442, "TGA"), (10, "TAA")], []],
            "alt_all_stop_codons": [[(3, "TGA")], None],
            "transcript_all_stop_codons": [None, [(7, "TAG")]],
        }
    )
    out = tmp_path / "stop_codons.parquet"

    write_results(df, str(out))

    stop_codon_type = pa.list_(pa.struct([pa.field("position", pa.int64()), pa.field("codon", pa.string())]))
    schema = pq.read_schema(str(out))
    for column in ["ref_all_stop_codons", "alt_all_stop_codons", "transcript_all_stop_codons"]:
        assert schema.field(column).type.equals(stop_codon_type)

    table = pq.read_table(out)
    assert table.column("ref_all_stop_codons").to_pylist() == [
        [{"position": 5442, "codon": "TGA"}, {"position": 10, "codon": "TAA"}],
        [],
    ]
    assert table.column("alt_all_stop_codons").to_pylist() == [[{"position": 3, "codon": "TGA"}], None]
    assert table.column("transcript_all_stop_codons").to_pylist() == [None, [{"position": 7, "codon": "TAG"}]]

    # write_results (the CSV path) and the in-memory df passed in are untouched
    assert df["ref_all_stop_codons"].iloc[0] == [(5442, "TGA"), (10, "TAA")]


def test_to_parquet_safe_leaves_other_columns_untouched():
    df = pd.DataFrame({"transcript_id": ["t1"], "nmd_escape": [True]})
    safe = to_parquet_safe(df)
    assert safe is df
