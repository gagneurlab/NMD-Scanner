import logging
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd
import pytest

import nmd_scanner.cli as cli_module
from nmd_scanner.cli import (
    annotate,
    is_valid_output_path,
    main,
    main_cli,
    write_results,
)
from nmd_scanner.schema import (
    MODEL_INPUTS,
    MODEL_STATUSES,
    NMD_RULE_COLUMN_KINDS,
    OUTPUT_COLUMN_KINDS,
    SEQUENCE_COLUMNS,
    STOP_CODON_COLUMNS,
    output_column_kinds,
    to_arrow,
)
from nmd_scanner.variant_placement import EXON_BOUNDARY_AMBIGUOUS, SPLICE_SITE_DESTROYED

RESOURCES = Path(__file__).resolve().parent.parent / "resources"


@pytest.fixture(autouse=True)
def _no_process_setup(monkeypatch):
    """main_cli sets up the process it runs in, but here that is the process of pytest."""
    monkeypatch.setattr(cli_module, "_set_up_process", lambda: None)


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
        annotation_path="resources/chr18.gff3.gz",
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

    ``ref_stop_codons`` and ``alt_stop_codons`` hold (position, codon) tuples,
    e.g. (5442, "TGA"); pyarrow cannot infer a single type for a tuple mixing int and str,
    so they need a typed struct schema instead. ``transcript_exons`` holds
    (exon_number, exon_length) tuples, both ints.
    """

    pytest.importorskip("pyarrow")
    import pyarrow as pa
    import pyarrow.parquet as pq

    out = tmp_path / "typed_results.parquet"
    results = main(
        vcf_path="resources/test_files/test_variants.vcf",
        annotation_path="resources/chr18.gff3.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(out),
    )

    assert out.exists()
    schema = pq.read_schema(str(out))

    stop_codon_type = pa.list_(pa.struct([pa.field("position", pa.int64()), pa.field("codon", pa.string())]))
    for column in ["ref_stop_codons", "alt_stop_codons"]:
        assert schema.field(column).type.equals(stop_codon_type), (
            f"{column} has unexpected parquet type: {schema.field(column).type}"
        )

    loaded = pd.read_parquet(out)
    assert len(loaded) == len(results)

    # Values are preserved, just reshaped from (position, codon) tuples to records
    has_stop_codons = results["ref_stop_codons"].apply(lambda v: isinstance(v, list) and len(v) > 0)
    sample_pos = results.index[has_stop_codons][0]
    expected = [{"position": pos, "codon": codon} for pos, codon in results.loc[sample_pos, "ref_stop_codons"]]
    assert list(loaded.loc[sample_pos, "ref_stop_codons"]) == expected

    # Exon numbers in transcript_exons are one type (int), not a mix of int and str
    exon_info_samples = loaded["transcript_exons"].dropna()
    exon_info_samples = exon_info_samples[exon_info_samples.apply(len) > 0]
    assert not exon_info_samples.empty
    for exon_number, exon_length in exon_info_samples.iloc[0]:
        assert not isinstance(exon_number, str)
        assert not isinstance(exon_length, str)


def test_write_results_parquet_types_stop_codon_columns(tmp_path):
    """
    ``to_arrow`` (used by ``write_results``) turns (position, codon) tuples into
    {"position": ..., "codon": ...} records for every stop-codon column, including
    ``alt_scan_stop_codons`` and rows holding None, so Parquet gets a typed struct
    schema instead of raising ArrowInvalid.
    """

    pytest.importorskip("pyarrow")
    import pyarrow as pa
    import pyarrow.parquet as pq

    df = pd.DataFrame(
        {
            "transcript_id": ["t1", "t2"],
            "ref_stop_codons": [[(5442, "TGA"), (10, "TAA")], []],
            "alt_stop_codons": [[(3, "TGA")], None],
            "alt_scan_stop_codons": [None, [(7, "TAG")]],
        }
    )
    out = tmp_path / "stop_codons.parquet"

    write_results(df, str(out))

    stop_codon_type = pa.list_(pa.struct([pa.field("position", pa.int64()), pa.field("codon", pa.string())]))
    schema = pq.read_schema(str(out))
    for column in ["ref_stop_codons", "alt_stop_codons", "alt_scan_stop_codons"]:
        assert schema.field(column).type.equals(stop_codon_type)

    table = pq.read_table(out)
    assert table.column("ref_stop_codons").to_pylist() == [
        [{"position": 5442, "codon": "TGA"}, {"position": 10, "codon": "TAA"}],
        [],
    ]
    assert table.column("alt_stop_codons").to_pylist() == [[{"position": 3, "codon": "TGA"}], None]
    assert table.column("alt_scan_stop_codons").to_pylist() == [None, [{"position": 7, "codon": "TAG"}]]

    # write_results (the CSV path) and the in-memory df passed in are untouched
    assert df["ref_stop_codons"].iloc[0] == [(5442, "TGA"), (10, "TAA")]


def test_the_deprecated_aliases_parquet_schema_and_to_parquet_safe_still_work():
    no_stop_codons = pd.DataFrame({"transcript_id": ["t1"], "nmd_escape": [True]})
    with pytest.warns(DeprecationWarning, match="to_arrow"):
        assert cli_module.to_parquet_safe(no_stop_codons) is no_stop_codons

    df = pd.DataFrame({"transcript_id": ["t1"], "transcript_exons": [[(1, 36)]], "ref_stop_codons": [[(5442, "TGA")]]})
    original = df.copy()
    with pytest.warns(DeprecationWarning, match="to_arrow"):
        safe = cli_module.to_parquet_safe(df)
    with pytest.warns(DeprecationWarning, match="to_arrow"):
        schema = cli_module.parquet_schema(df)

    assert safe is not df
    pd.testing.assert_frame_equal(safe.drop(columns="ref_stop_codons"), df.drop(columns="ref_stop_codons"))
    assert safe["ref_stop_codons"].tolist() == [[{"position": 5442, "codon": "TGA"}]]
    pd.testing.assert_frame_equal(df, original)
    assert schema.equals(to_arrow(df).schema)
    assert cli_module.STOP_CODON_COLUMNS == STOP_CODON_COLUMNS


@pytest.mark.parametrize("missing", [None, float("nan"), pd.NA])
def test_write_results_parquet_keeps_any_missing_stop_codon_value_null(tmp_path, missing):
    """None, np.nan and pd.NA in a stop-codon column are all written as null."""

    pytest.importorskip("pyarrow")
    import pyarrow.parquet as pq

    df = pd.DataFrame({"transcript_id": ["t1", "t2"], "alt_stop_codons": [[(3, "TGA")], missing]})
    out = tmp_path / "missing.parquet"

    write_results(df, str(out))

    assert pq.read_table(out).column("alt_stop_codons").to_pylist() == [[{"position": 3, "codon": "TGA"}], None]


def _results_schema(tmp_path, vcf_path, name):
    import pyarrow.parquet as pq

    out = tmp_path / name
    main(
        vcf_path=vcf_path,
        annotation_path="resources/chr18.gff3.gz",
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
    assert schema_with.field("alt_scan_stop_codons").type.equals(
        pa.list_(pa.struct([pa.field("position", pa.int64()), pa.field("codon", pa.string())]))
    )
    assert schema_with.field("alt_scan_stop_codon_exons").type.equals(pa.list_(pa.int64()))
    assert schema_with.field("alt_scan_start_codon_exon").type.equals(pa.int64())
    assert schema_with.field("alt_scan_first_stop_codon").type.equals(pa.string())


def test_parquet_values_roundtrip_unchanged_and_none_stays_null(tmp_path):
    pytest.importorskip("pyarrow")
    import pyarrow.parquet as pq

    out = tmp_path / "roundtrip.parquet"
    results = main(
        vcf_path="resources/test_files/variants.vcf",
        annotation_path="resources/chr18.gff3.gz",
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
            elif column in ("ref_stop_codons", "alt_stop_codons", "alt_scan_stop_codons"):
                assert act == [{"position": p, "codon": c} for p, c in exp], column
            elif isinstance(exp, list):
                assert [list(x) if isinstance(x, tuple) else x for x in exp] == act, column
            else:
                assert exp == act, column


@pytest.mark.parametrize("ids", [("12345", "NA"), ("007", "0123")])
@pytest.mark.parametrize("output_name", ["numeric_ids.parquet", "numeric_ids.csv"])
def test_main_keeps_numeric_and_NA_variant_ids_as_text(tmp_path, output_name, ids):
    """
    ClinVar-style numeric IDs, also with leading zeros, and an ID of "NA" must reach variant_id as written, in either
    output format.
    """

    vcf = tmp_path / "numeric_ids.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr18\t21383521\t{ids[0]}\tG\tGT\t.\t.\t.\n"
        f"chr18\t21383521\t{ids[1]}\tG\tGTTT\t.\t.\t.\n"
    )
    out = tmp_path / output_name
    results = main(
        vcf_path=str(vcf),
        annotation_path="resources/chr18.gff3.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(out),
    )

    assert set(results["variant_id"]) == set(ids)
    if output_name.endswith(".parquet"):
        loaded = pd.read_parquet(out)
    else:
        loaded = pd.read_csv(out, dtype={"variant_id": str}, keep_default_na=False)
    assert set(loaded["variant_id"]) == set(ids)


def test_main_without_cds_overlap_writes_empty_csv(tmp_path, intergenic_vcf, caplog):
    out = tmp_path / "empty.csv"
    with caplog.at_level(logging.INFO):
        results = main(
            vcf_path=intergenic_vcf,
            annotation_path="resources/chr18.gff3.gz",
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
    """
    The full run has rows with an unknown alt transcript. Their start and stop loss and NMD rules are null, and
    those columns stay bool in Parquet, as in the empty run.
    """

    pytest.importorskip("pyarrow")
    import pyarrow as pa
    import pyarrow.parquet as pq

    empty_out = tmp_path / "empty.parquet"
    main(
        vcf_path=intergenic_vcf,
        annotation_path="resources/chr18.gff3.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(empty_out),
    )
    full_out = tmp_path / "full.parquet"
    main(
        vcf_path="resources/test_files/variants.vcf",
        annotation_path="resources/chr18.gff3.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(full_out),
    )

    assert pq.read_table(empty_out).num_rows == 0
    schema = pq.read_schema(str(full_out))
    assert pq.read_schema(str(empty_out)).equals(schema)

    loaded = pd.read_parquet(full_out)
    unknown = loaded["unknown_reason"].notna()
    assert unknown.any() and not unknown.all()
    assert set(loaded.loc[unknown, "unknown_reason"]) <= {SPLICE_SITE_DESTROYED, EXON_BOUNDARY_AMBIGUOUS}
    assert schema.field("unknown_reason").type.equals(pa.string())
    for column in ["start_loss", "stop_loss", *NMD_RULE_COLUMN_KINDS]:
        assert schema.field(column).type.equals(pa.bool_()), column
        assert loaded.loc[unknown, column].isna().all(), column
        assert loaded.loc[~unknown, column].notna().any(), column
    assert (loaded.loc[unknown, "nmd_model_status"] == "unknown_effect").all()
    assert (loaded.loc[~unknown, "nmd_model_status"] != "unknown_effect").all()


def test_main_without_reference_mismatches_does_not_warn_about_them(tmp_path, caplog):
    vcf = tmp_path / "matching.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr18\t21383521\tv1\tG\tGT\t.\t.\t.\n"
    )
    with caplog.at_level(logging.INFO):
        results = main(
            vcf_path=str(vcf),
            annotation_path="resources/chr18.gff3.gz",
            fasta_path="resources/chr18.fa.gz",
            output=str(tmp_path / "out.csv"),
        )
    assert not results.empty
    assert "reference mismatches" not in caplog.text


def test_main_with_a_reference_mismatch_warns_about_it(tmp_path, caplog):
    # v1 matches the FASTA. "mismatch" has REF C, but chr18:21383519 is A, the first base of the GREB1L CDS.
    vcf = tmp_path / "one_mismatch.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr18\t21383519\tmismatch\tC\tT\t.\t.\t.\n"
        "chr18\t21383521\tv1\tG\tGT\t.\t.\t.\n"
    )
    with caplog.at_level(logging.INFO):
        results = main(
            vcf_path=str(vcf),
            annotation_path="resources/chr18.gff3.gz",
            fasta_path="resources/chr18.fa.gz",
            output=str(tmp_path / "out.csv"),
        )

    assert set(results["variant_id"]) == {"v1"}
    warnings = [record.getMessage() for record in caplog.records if record.levelno == logging.WARNING]
    # The count is not checked.
    assert any(message.endswith(" variant-transcript pairs due to reference mismatches.") for message in warnings)
    listing = next(message for message in warnings if message.startswith("Reference-mismatched variants:"))
    # one row per transcript: transcript_id, Chromosome, Start_variant, End_variant, Ref, Alt
    rows = [line.split() for line in listing.splitlines()[2:]]
    assert rows
    assert {tuple(row[1:]) for row in rows} == {("chr18", "21383518", "21383519", "C", "T")}


def test_main_with_only_reference_mismatches_writes_empty_csv(tmp_path, reference_mismatch_vcf, caplog):
    out = tmp_path / "mismatch.csv"
    with caplog.at_level(logging.INFO):
        results = main(
            vcf_path=reference_mismatch_vcf,
            annotation_path="resources/chr18.gff3.gz",
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
        annotation_path="resources/chr18.gff3.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(empty_out),
    )
    full_out = tmp_path / "full.parquet"
    main(
        vcf_path="resources/test_files/variants.vcf",
        annotation_path="resources/chr18.gff3.gz",
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
        annotation_path="resources/chr18.gff3.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(out),
        reassign_exons=True,
    )

    assert not results.empty
    assert out.exists()
    for exon_info in results["transcript_exons"]:
        for exon_number, exon_length in exon_info:
            assert not isinstance(exon_number, str)
            assert not isinstance(exon_length, str)

    # chr18.gff3.gz is hg38, where the annotated exon numbers already follow transcript order.
    # extract_ptc casts the annotated ones to int, too.
    annotated = main(
        vcf_path="resources/test_files/test_variants.vcf",
        annotation_path="resources/chr18.gff3.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(tmp_path / "annotated.csv"),
    )
    pd.testing.assert_frame_equal(results, annotated)


def test_main_keeps_variants_with_unknown_alt_transcript(tmp_path):
    """Variants over a splice site get a row without prediction, and unknown_reason says why."""

    results = main(
        vcf_path="resources/test_files/variants.vcf",
        annotation_path="resources/chr18.gff3.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(tmp_path / "results.csv"),
    )

    assert list(results.columns) == list(OUTPUT_COLUMN_KINDS)
    unknown = results[results["unknown_reason"].notna()]
    assert not unknown.empty
    assert set(unknown["unknown_reason"]) <= {SPLICE_SITE_DESTROYED, EXON_BOUNDARY_AMBIGUOUS}
    assert unknown["ref_cds_seq"].notna().all()
    for column in ["alt_cds_seq", "alt_has_ptc", "start_loss", "stop_loss", "nmd_escape"]:
        assert unknown[column].isna().all(), column
    assert results.loc[results["unknown_reason"].isna(), "nmd_escape"].notna().all()


@pytest.fixture(scope="module")
def cli_stderr(tmp_path_factory):
    """stderr of the CLI in its own process, on the bundled chr18 test data"""

    output = tmp_path_factory.mktemp("cli") / "results.csv"
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            "from nmd_scanner.cli import main_cli; main_cli()",
            "--vcf",
            str(RESOURCES / "test_files" / "test_variants.vcf"),
            "--annotation",
            str(RESOURCES / "chr18.gff3.gz"),
            "--fasta",
            str(RESOURCES / "chr18.fa.gz"),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        check=True,
        env={name: value for name, value in os.environ.items() if name != "TQDM_DISABLE"},
    )
    return completed.stderr


def test_main_cli_logs_the_info_messages_of_nmd_scanner_in_its_format(cli_stderr):
    assert "INFO nmd_scanner.cli: Reading VCF file" in cli_stderr
    assert "WARNING nmd_scanner.rules: Skipping 88 variant-transcript pairs" in cli_stderr
    # polars-bio's Rust code logs at INFO too
    info = [line for line in cli_stderr.splitlines() if " INFO " in line]
    assert all(" INFO nmd_scanner." in line for line in info)


def test_main_cli_shows_no_progress_bars(cli_stderr):
    # polars-bio shows a tqdm bar for every read
    assert "rows/s" not in cli_stderr


# --annotation CLI option tests


def _patch_main(monkeypatch):
    calls = {}

    def fake_main(vcf_path, annotation_path, fasta_path, output, reassign_exons=False, sequences=True):
        calls["annotation_path"] = annotation_path
        calls["sequences"] = sequences
        return pd.DataFrame()

    monkeypatch.setattr(cli_module, "main", fake_main)
    return calls


def test_main_cli_requires_annotation(monkeypatch, tmp_path, capsys):
    out = tmp_path / "out.csv"
    monkeypatch.setattr(sys, "argv", ["nmd-scanner", "--vcf", "in.vcf", "--fasta", "ref.fa", "--output", str(out)])
    with pytest.raises(SystemExit):
        main_cli()
    assert "--annotation" in capsys.readouterr().err


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


@pytest.mark.parametrize("flags, sequences", [([], True), (["--no-sequences"], False)])
def test_main_cli_no_sequences_flag_reaches_main(monkeypatch, tmp_path, flags, sequences):
    out = tmp_path / "out.csv"
    calls = _patch_main(monkeypatch)
    argv = ["nmd-scanner", "--vcf", "in.vcf", "--annotation", "annotation.gff3", "--fasta", "ref.fa", "--output"]
    monkeypatch.setattr(sys, "argv", [*argv, str(out), *flags])
    main_cli()
    assert calls["sequences"] is sequences


def test_main_cli_rejects_empty_annotation(monkeypatch, tmp_path, capsys):
    out = tmp_path / "out.csv"
    monkeypatch.setattr(
        sys, "argv", ["nmd-scanner", "--vcf", "in.vcf", "--annotation=", "--fasta", "ref.fa", "--output", str(out)]
    )
    with pytest.raises(SystemExit):
        main_cli()
    assert "--annotation" in capsys.readouterr().err


def test_main_cli_rejects_an_annotation_of_unknown_format(monkeypatch, tmp_path, capsys):
    out = tmp_path / "out.csv"
    calls = _patch_main(monkeypatch)
    monkeypatch.setattr(
        sys,
        "argv",
        ["nmd-scanner", "--vcf", "in.vcf", "--annotation", "annotation.txt", "--fasta", "ref.fa", "--output", str(out)],
    )
    with pytest.raises(SystemExit):
        main_cli()
    assert "Cannot detect annotation format" in capsys.readouterr().err
    assert calls == {}


def test_main_accepts_a_path_object_for_the_annotation(tmp_path):
    out = str(tmp_path / "out.csv")
    vcf = "resources/part-00241-61a0abbf-fbf9-444f-8287-4e46ad4b9b7b-c000.vcf"
    via_path = main(vcf, Path("resources/chr18.gff3.gz"), "resources/chr18.fa.gz", out)
    via_str = main(vcf, "resources/chr18.gff3.gz", "resources/chr18.fa.gz", out)
    assert not via_path.empty
    pd.testing.assert_frame_equal(via_path, via_str)


_TUBB8B_GFF3 = """\
##gff-version 3
chr18\tensembl\tgene\t47221\t49615\t.\t-\t.\tID=gene:G1;biotype=protein_coding
chr18\tensembl\tmRNA\t47221\t49615\t.\t-\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding
chr18\tensembl\texon\t49501\t49615\t.\t-\t.\tParent=transcript:T1;rank=1
chr18\tensembl\texon\t49129\t49237\t.\t-\t.\tParent=transcript:T1;rank=2
chr18\tensembl\texon\t48940\t49050\t.\t-\t.\tParent=transcript:T1;rank=3
chr18\tensembl\texon\t47221\t48447\t.\t-\t.\tParent=transcript:T1;rank=4
chr18\tensembl\tCDS\t49501\t49557\t.\t-\t0\tID=CDS:P1;Parent=transcript:T1
chr18\tensembl\tCDS\t49129\t49237\t.\t-\t0\tID=CDS:P1;Parent=transcript:T1
chr18\tensembl\tCDS\t48940\t49050\t.\t-\t2\tID=CDS:P1;Parent=transcript:T1
chr18\tensembl\tCDS\t47390\t48447\t.\t-\t2\tID=CDS:P1;Parent=transcript:T1
"""


def test_annotate_returns_what_main_writes(tmp_path):
    out = tmp_path / "main.csv"
    expected = main(
        vcf_path="resources/test_files/test_variants.vcf",
        annotation_path="resources/chr18.gff3.gz",
        fasta_path="resources/chr18.fa.gz",
        output=str(out),
    )

    results = annotate(
        "resources/test_files/test_variants.vcf",
        "resources/chr18.gff3.gz",
        "resources/chr18.fa.gz",
    )

    assert isinstance(results, pd.DataFrame)
    assert not results.empty
    assert list(results.columns) == list(expected.columns)
    pd.testing.assert_frame_equal(results, expected)
    written = tmp_path / "annotate.csv"
    write_results(results, str(written))
    assert written.read_bytes() == out.read_bytes()


def test_annotate_does_not_write_files_or_configure_logging(tmp_path):
    # in its own process, because pytest has imported polars-bio already, which configures the root logger
    code = textwrap.dedent(
        f"""
        import logging

        root = logging.getLogger()
        root.setLevel(logging.INFO)
        import nmd_scanner

        nmd_scanner.annotate(
            {str(RESOURCES / "test_files" / "test_variants.vcf")!r},
            {str(RESOURCES / "chr18.gff3.gz")!r},
            {str(RESOURCES / "chr18.fa.gz")!r},
        )
        print(len(root.handlers), root.level)
        """
    )

    completed = subprocess.run([sys.executable, "-c", code], cwd=tmp_path, capture_output=True, text=True, check=True)

    assert list(tmp_path.iterdir()) == []
    assert completed.stdout.split() == ["0", str(logging.INFO)]


# Each symbolic allele and breakend sits at the position of v1, inside the GREB1L CDS. Before they were skipped, a
# symbolic allele with a padding base went into the alt CDS as text, e.g. "<DEL>".
SYMBOLIC_ALTS = ["<DEL>", "<DUP>", "<INS>", "<INV>", "<CNV>", "<DUP:TANDEM>", "G]chr2:100]", "[chr2:100[G", "G.", ".G"]


def _symbolic_vcf(tmp_path, with_snv):
    records = [f"chr18\t21383521\tsv{i}\tG\t{alt}\t.\t.\t.\n" for i, alt in enumerate(SYMBOLIC_ALTS)]
    if with_snv:
        records.append("chr18\t21383521\tv1\tG\tGT\t.\t.\t.\n")
    vcf = tmp_path / "symbolic.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "".join(records))
    return str(vcf)


def test_annotate_skips_symbolic_alleles_and_breakends_with_a_warning(tmp_path, caplog):
    with caplog.at_level(logging.INFO):
        results = annotate(_symbolic_vcf(tmp_path, True), "resources/chr18.gff3.gz", "resources/chr18.fa.gz")

    assert not results.empty
    assert set(results["variant_id"]) == {"v1"}
    assert not results["alt_cds_seq"].str.contains("<|>|\\[|\\]|\\.").any()
    warnings = [record.getMessage() for record in caplog.records if record.levelno == logging.WARNING]
    assert any(
        message.startswith(f"Skipping {len(SYMBOLIC_ALTS)} variant(s) with a symbolic ALT") for message in warnings
    )


def test_annotate_without_symbolic_alleles_does_not_warn_about_them(caplog):
    with caplog.at_level(logging.INFO):
        annotate("resources/test_files/test_variants.vcf", "resources/chr18.gff3.gz", "resources/chr18.fa.gz")

    assert "symbolic ALT" not in caplog.text


def test_annotate_raises_if_the_fasta_lacks_a_chromosome_with_variants_and_cds_rows(tmp_path):
    fasta = tmp_path / "chr1.fa"
    fasta.write_text(">chr1\nACGT\n")

    with pytest.raises(ValueError) as error:
        annotate("resources/test_files/test_variants.vcf", "resources/chr18.gff3.gz", str(fasta))

    assert str(error.value) == (
        "The FASTA has no sequence for 1 chromosome(s) with variants and CDS rows: chr18. Each chromosome with "
        "variants and CDS rows needs a sequence of the same name in the FASTA."
    )


def test_annotate_ignores_a_chromosome_without_cds_rows_that_the_fasta_lacks(tmp_path):
    vcf = tmp_path / "chrZ.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrZ\t1000\tz1\tA\tT\t.\t.\t.\n"
    )

    results = annotate(str(vcf), "resources/chr18.gff3.gz", "resources/chr18.fa.gz")

    assert results.empty
    assert list(results.columns) == list(OUTPUT_COLUMN_KINDS)


def test_annotate_sets_nmd_model_status_ok_exactly_for_the_rows_the_model_can_score():
    results = annotate("resources/test_files/test_variants.vcf", "resources/chr18.gff3.gz", "resources/chr18.fa.gz")

    new_ptc = results["alt_has_ptc"].fillna(False) & ~results["ref_has_ptc"].fillna(False)
    scorable = new_ptc & results[MODEL_INPUTS].notna().all(axis=1)
    assert (results["nmd_model_status"] == "ok").tolist() == scorable.tolist()
    assert scorable.any()
    assert set(results["nmd_model_status"]) <= set(MODEL_STATUSES)


def test_annotate_reassign_exons_matches_main(tmp_path):
    args = ("resources/test_files/test_variants.vcf", "resources/chr18.gff3.gz", "resources/chr18.fa.gz")
    expected = main(*args, str(tmp_path / "main.csv"), reassign_exons=True)

    results = annotate(*args, reassign_exons=True)

    pd.testing.assert_frame_equal(results, expected)


def test_annotate_reads_gff3_with_the_fasta(tmp_path):
    from pyfaidx import Fasta

    fasta_path = "resources/chr18.fa.gz"
    ref = str(Fasta(fasta_path)["chr18"][48000:48001]).upper()
    alt = "A" if ref != "A" else "C"
    vcf = tmp_path / "in_cds.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr18\t48001\tin_cds\t{ref}\t{alt}\t.\t.\t.\n"
    )
    gff3 = tmp_path / "tubb8b.gff3"
    gff3.write_text(_TUBB8B_GFF3)

    results = annotate(str(vcf), str(gff3), fasta_path)
    expected = main(str(vcf), str(gff3), fasta_path, str(tmp_path / "main.csv"))

    assert list(results.columns) == list(OUTPUT_COLUMN_KINDS)
    assert list(results["transcript_id"]) == ["T1"]
    # The FASTA has TAG (reverse complement of chr18:47390-47392) at the end of the last CDS row, so the
    # coding region ends in a stop codon: 57 + 109 + 111 + 1058 nt.
    assert results["strand"].tolist() == ["-"]
    assert results["has_stop_codon"].tolist() == [True]
    assert results["ref_cds_length"].tolist() == [1335]
    pd.testing.assert_frame_equal(results, expected)


@pytest.mark.parametrize("extension", [".csv", ".parquet"])
@pytest.mark.parametrize("vcf", ["resources/test_files/variants.vcf", "intergenic_vcf"])
def test_main_without_sequences_writes_the_full_output_without_the_4_sequence_columns(
    request, tmp_path, vcf, extension
):
    import pyarrow.parquet as pq

    vcf_path = request.getfixturevalue(vcf) if vcf.endswith("_vcf") else vcf
    args = (vcf_path, "resources/chr18.gff3.gz", "resources/chr18.fa.gz")
    full_out = tmp_path / f"full{extension}"
    reduced_out = tmp_path / f"reduced{extension}"
    write_results(annotate(*args), str(full_out))

    reduced = main(*args, str(reduced_out), sequences=False)

    assert list(reduced.columns) == list(output_column_kinds(sequences=False))
    if extension == ".csv":
        loaded = pd.read_csv(reduced_out)
        assert list(loaded.columns) == list(reduced.columns)
        pd.testing.assert_frame_equal(loaded, pd.read_csv(full_out).drop(columns=list(SEQUENCE_COLUMNS)))
    else:
        table = pq.read_table(reduced_out)
        assert table.column_names == list(reduced.columns)
        assert table.equals(pq.read_table(full_out).drop_columns(list(SEQUENCE_COLUMNS)))
