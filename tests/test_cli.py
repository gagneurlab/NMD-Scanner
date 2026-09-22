import logging
import sys
from pathlib import Path

import pandas as pd
import pytest

import nmd_scanner.cli as cli_module
from nmd_scanner.cli import OUTPUT_COLUMN_KINDS, is_valid_output_path, main, main_cli, to_parquet_safe, write_results


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
    assert to_parquet_safe(df) is df

    df = pd.DataFrame(
        {"transcript_id": ["t1"], "transcript_exon_info": [[(1, 36)]], "ref_all_stop_codons": [[(5442, "TGA")]]}
    )
    original = df.copy()
    safe = to_parquet_safe(df)
    assert safe is not df
    pd.testing.assert_frame_equal(safe.drop(columns="ref_all_stop_codons"), df.drop(columns="ref_all_stop_codons"))
    assert safe["ref_all_stop_codons"].tolist() == [[{"position": 5442, "codon": "TGA"}]]
    pd.testing.assert_frame_equal(df, original)


@pytest.mark.parametrize("missing", [None, float("nan"), pd.NA])
def test_write_results_parquet_keeps_any_missing_stop_codon_value_null(tmp_path, missing):
    """None, np.nan and pd.NA in a stop-codon column are all written as null."""

    pytest.importorskip("pyarrow")
    import pyarrow.parquet as pq

    df = pd.DataFrame({"transcript_id": ["t1", "t2"], "alt_all_stop_codons": [[(3, "TGA")], missing]})
    out = tmp_path / "missing.parquet"

    write_results(df, str(out))

    assert pq.read_table(out).column("alt_all_stop_codons").to_pylist() == [[{"position": 3, "codon": "TGA"}], None]


def _results_schema(tmp_path, vcf_path, name):
    import pyarrow.parquet as pq

    out = tmp_path / name
    main(
        vcf_path=vcf_path,
        gtf_path="resources/chr18.gtf.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(out),
    )
    return pq.read_schema(str(out)), pd.read_parquet(out)


def test_parquet_schema_is_the_same_for_every_run(tmp_path):
    """
    Without an explicit schema, columns that are only None in a run (e.g. the transcript_*
    stop-codon columns when no variant has a start or stop loss) are written as ``null``.
    A run without start or stop loss, a run with them and an empty table must agree.
    """

    pytest.importorskip("pyarrow")
    import pyarrow as pa

    schema_without, _ = _results_schema(tmp_path, "resources/test_files/test_variants.vcf", "without.parquet")
    schema_with, loaded = _results_schema(tmp_path, "resources/test_files/variants.vcf", "with.parquet")
    assert loaded["start_loss"].any() or loaded["stop_loss"].any()

    empty = pd.DataFrame(columns=schema_with.names)
    empty_out = tmp_path / "empty.parquet"
    write_results(empty, str(empty_out))
    import pyarrow.parquet as pq

    schema_empty = pq.read_schema(str(empty_out))

    assert schema_without.equals(schema_with)
    assert schema_empty.equals(schema_with)
    for field in schema_with:
        assert not pa.types.is_null(field.type), field.name
    assert schema_with.field("transcript_all_stop_codons").type.equals(
        pa.list_(pa.struct([pa.field("position", pa.int64()), pa.field("codon", pa.string())]))
    )
    assert schema_with.field("transcript_stop_codon_exons").type.equals(pa.list_(pa.int64()))
    assert schema_with.field("transcript_start_codon_exon").type.equals(pa.int64())
    assert schema_with.field("transcript_valid_stop").type.equals(pa.bool_())


def test_parquet_values_roundtrip_unchanged_and_none_stays_null(tmp_path):
    pytest.importorskip("pyarrow")
    import pyarrow.parquet as pq

    out = tmp_path / "roundtrip.parquet"
    results = main(
        vcf_path="resources/test_files/variants.vcf",
        gtf_path="resources/chr18.gtf.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(out),
    )
    table = pq.read_table(out)

    assert table.column_names == list(results.columns)
    for column in results.columns:
        expected = results[column].tolist()
        actual = table.column(column).to_pylist()
        assert len(actual) == len(expected)
        for exp, act in zip(expected, actual):
            if pd.api.types.is_scalar(exp) and pd.isna(exp):
                assert act is None, column
            elif column in ("ref_all_stop_codons", "alt_all_stop_codons", "transcript_all_stop_codons"):
                assert act == [{"position": p, "codon": c} for p, c in exp], column
            elif isinstance(exp, list):
                assert [list(x) if isinstance(x, tuple) else x for x in exp] == act, column
            else:
                assert exp == act, column


def test_main_without_cds_overlap_writes_empty_csv(tmp_path, intergenic_vcf, caplog):
    out = tmp_path / "empty.csv"
    with caplog.at_level(logging.INFO):
        results = main(
            vcf_path=intergenic_vcf,
            gtf_path="resources/chr18.gtf.gz",
            fasta_path="resources/chr18.fa.gz",
            output=str(out),
        )

    assert results.empty
    assert list(results.columns) == list(OUTPUT_COLUMN_KINDS)
    loaded = pd.read_csv(out)
    assert list(loaded.columns) == list(OUTPUT_COLUMN_KINDS)
    assert len(loaded) == 0
    assert "No variant overlapped a CDS" in caplog.text


def test_main_without_cds_overlap_writes_empty_parquet_with_the_usual_schema(tmp_path, intergenic_vcf):
    pytest.importorskip("pyarrow")
    import pyarrow.parquet as pq

    empty_out = tmp_path / "empty.parquet"
    main(
        vcf_path=intergenic_vcf,
        gtf_path="resources/chr18.gtf.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(empty_out),
    )
    full_out = tmp_path / "full.parquet"
    main(
        vcf_path="resources/test_files/variants.vcf",
        gtf_path="resources/chr18.gtf.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(full_out),
    )

    assert pq.read_table(empty_out).num_rows == 0
    assert pq.read_schema(str(empty_out)).equals(pq.read_schema(str(full_out)))


def test_main_with_only_reference_mismatches_writes_empty_csv(tmp_path, reference_mismatch_vcf, caplog):
    out = tmp_path / "mismatch.csv"
    with caplog.at_level(logging.INFO):
        results = main(
            vcf_path=reference_mismatch_vcf,
            gtf_path="resources/chr18.gtf.gz",
            fasta_path="resources/chr18.fa.gz",
            output=str(out),
        )

    assert results.empty
    loaded = pd.read_csv(out)
    assert list(loaded.columns) == list(OUTPUT_COLUMN_KINDS)
    assert len(loaded) == 0
    assert "No variant left after the reference check" in caplog.text


def test_main_with_only_reference_mismatches_writes_empty_parquet_with_the_usual_schema(
    tmp_path, reference_mismatch_vcf
):
    pytest.importorskip("pyarrow")
    import pyarrow.parquet as pq

    empty_out = tmp_path / "mismatch.parquet"
    main(
        vcf_path=reference_mismatch_vcf,
        gtf_path="resources/chr18.gtf.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(empty_out),
    )
    full_out = tmp_path / "full.parquet"
    main(
        vcf_path="resources/test_files/variants.vcf",
        gtf_path="resources/chr18.gtf.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(full_out),
    )

    assert pq.read_table(empty_out).num_rows == 0
    assert pq.read_schema(str(empty_out)).equals(pq.read_schema(str(full_out)))


def test_main_end_to_end_reassign_exons(tmp_path):
    """``--reassign_exons`` runs on the bundled chr18 data and yields int exon numbers."""

    out = tmp_path / "reassigned.csv"
    results = main(
        vcf_path="resources/test_files/test_variants.vcf",
        gtf_path="resources/chr18.gtf.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(out),
        reassign_exons=True,
    )

    assert not results.empty
    assert out.exists()
    for exon_info in results["transcript_exon_info"]:
        for exon_number, exon_length in exon_info:
            assert not isinstance(exon_number, str)
            assert not isinstance(exon_length, str)

    # chr18.gtf.gz is hg38, where the annotated exon numbers already follow transcript order.
    # extract_ptc casts the annotated ones to int, too.
    annotated = main(
        vcf_path="resources/test_files/test_variants.vcf",
        gtf_path="resources/chr18.gtf.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(tmp_path / "annotated.csv"),
    )
    pd.testing.assert_frame_equal(results, annotated)


# --annotation / --gtf CLI option tests


def _patch_main(monkeypatch):
    calls = {}

    def fake_main(vcf_path, gtf_path, fasta_path, output, reassign_exons=False, annotation_path=None):
        calls["gtf_path"] = gtf_path
        calls["annotation_path"] = annotation_path
        return pd.DataFrame()

    monkeypatch.setattr(cli_module, "main", fake_main)
    return calls


def test_main_cli_requires_annotation_or_gtf(monkeypatch, tmp_path):
    out = tmp_path / "out.csv"
    monkeypatch.setattr(sys, "argv", ["nmd-scanner", "--vcf", "in.vcf", "--fasta", "ref.fa", "--output", str(out)])
    with pytest.raises(SystemExit):
        main_cli()


def test_main_cli_rejects_both_annotation_and_gtf(monkeypatch, tmp_path):
    out = tmp_path / "out.csv"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "nmd-scanner",
            "--vcf",
            "in.vcf",
            "--annotation",
            "a.gtf",
            "--gtf",
            "b.gtf",
            "--fasta",
            "ref.fa",
            "--output",
            str(out),
        ],
    )
    with pytest.raises(SystemExit):
        main_cli()


def test_main_cli_annotation_option_reaches_main(monkeypatch, tmp_path):
    out = tmp_path / "out.csv"
    calls = _patch_main(monkeypatch)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "nmd-scanner",
            "--vcf",
            "in.vcf",
            "--annotation",
            "annotation.gff3",
            "--fasta",
            "ref.fa",
            "--output",
            str(out),
        ],
    )
    main_cli()
    assert calls["annotation_path"] == "annotation.gff3"
    assert calls["gtf_path"] is None


def test_main_cli_gtf_alias_reaches_main(monkeypatch, tmp_path):
    """--gtf is a deprecated alias for --annotation; it must keep working unchanged."""
    out = tmp_path / "out.csv"
    calls = _patch_main(monkeypatch)
    monkeypatch.setattr(
        sys,
        "argv",
        ["nmd-scanner", "--vcf", "in.vcf", "--gtf", "annotation.gtf", "--fasta", "ref.fa", "--output", str(out)],
    )
    main_cli()
    assert calls["gtf_path"] == "annotation.gtf"
    assert calls["annotation_path"] is None


def test_main_cli_rejects_empty_annotation(monkeypatch, tmp_path, capsys):
    out = tmp_path / "out.csv"
    monkeypatch.setattr(
        sys, "argv", ["nmd-scanner", "--vcf", "in.vcf", "--annotation=", "--fasta", "ref.fa", "--output", str(out)]
    )
    with pytest.raises(SystemExit):
        main_cli()
    assert "--annotation" in capsys.readouterr().err


class _GtfReaderUsed(Exception):
    pass


def test_main_gtf_path_always_uses_the_gtf_reader(monkeypatch):
    """A GTF under any file name works through gtf_path (--gtf); a GFF3 name does not switch the reader."""

    def fake_read_gtf(path):
        raise _GtfReaderUsed(path)

    monkeypatch.setattr(cli_module, "read_gtf", fake_read_gtf)
    for name in ("plain.txt", "annotation.gff3"):
        with pytest.raises(_GtfReaderUsed, match=name):
            main(
                "resources/part-00241-61a0abbf-fbf9-444f-8287-4e46ad4b9b7b-c000.vcf",
                name,
                "resources/chr18.fa.gz",
                "unused.csv",
            )


def test_main_accepts_a_path_object_for_the_annotation(tmp_path):
    out = str(tmp_path / "out.csv")
    vcf = "resources/part-00241-61a0abbf-fbf9-444f-8287-4e46ad4b9b7b-c000.vcf"
    via_gtf_path = main(vcf, Path("resources/chr18.gtf.gz"), "resources/chr18.fa.gz", out)
    via_annotation = main(vcf, None, "resources/chr18.fa.gz", out, annotation_path=Path("resources/chr18.gtf.gz"))
    pd.testing.assert_frame_equal(via_gtf_path, via_annotation)


def test_main_needs_exactly_one_of_gtf_path_and_annotation_path():
    with pytest.raises(ValueError, match="exactly one"):
        main("in.vcf", None, "ref.fa", "out.csv")
    with pytest.raises(ValueError, match="exactly one"):
        main("in.vcf", "a.gtf", "ref.fa", "out.csv", annotation_path="b.gtf")
