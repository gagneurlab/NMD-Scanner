"""Tests of nmd_scanner.benchmark on tiny tables, without network."""

import hashlib
import io
import logging
import pkgutil
import zipfile

import fsspec
import pandas as pd
import polars as pl
import pytest
from polars.testing import assert_frame_equal
from pyfaidx import Fasta

from nmd_scanner.benchmark import DATASETS, datasets, join_features
from nmd_scanner.benchmark.datasets import geuvadis, gtex, mmrf_target, tcga
from nmd_scanner.benchmark.sources import check_sha256, fetch, hgvsc_alleles

# ENST3 is missing from the GFF3
TRANSCRIPTS = pl.DataFrame({"transcript": ["ENST1", "ENST2"], "strand": ["+", "-"], "gene": ["ENSG1", "ENSG2"]})
SCHEMA = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "ref": pl.String,
    "alt": pl.String,
    "transcript": pl.String,
    "gene": pl.String,
    "NMD_efficiency": pl.Float64,
}


@pytest.fixture
def fasta(tmp_path):
    """chr1 of 8 bases, with A at position 1 and T at position 8."""
    path = tmp_path / "genome.fa"
    path.write_text(">chr1\nACGTACGT\n")
    return Fasta(str(path))


def kim_table(**columns) -> pl.DataFrame:
    """
    A table of Kim et al. 2024 with 4 rows: a plus-strand transcript, a minus-strand transcript, and two
    rows on a transcript missing from the GFF3, whose REF matches the FASTA on the minus and on the plus strand.
    """
    return pl.DataFrame(
        {
            "chromosome": ["chr1"] * 4,
            "start": [2, 3, 4, 5],
            "Transcript_ID": ["ENST1", "ENST2", "ENST3", "ENST3"],
            "HGVSc": ["c.10C>T", "ENST2:c.7C>A", "c.1A>G", "c.2A>G"],
            "NMD_efficiency": [1.5, 0.5, 2.0, -0.25],
            **columns,
        }
    )


def benchmark(rows: dict, extra_schema: dict | None = None) -> pl.DataFrame:
    return pl.DataFrame(rows, schema={**SCHEMA, **(extra_schema or {}), "source_row": pl.UInt32})


def test_hgvsc_alleles_on_both_strands(fasta):
    table_df = hgvsc_alleles(kim_table(), TRANSCRIPTS, fasta)
    assert table_df.select("chrom", "pos", "transcript", "ref", "alt").rows() == [
        ("chr1", 2, "ENST1", "C", "T"),
        ("chr1", 3, "ENST2", "G", "T"),
        ("chr1", 4, "ENST3", "T", "C"),
        ("chr1", 5, "ENST3", "A", "G"),
    ]


def test_hgvsc_alleles_fail_on_other_variants(fasta):
    with pytest.raises(ValueError, match="c.10del .ENST1. is not a single-base substitution"):
        hgvsc_alleles(kim_table().with_columns(HGVSc=pl.lit("c.10del")), TRANSCRIPTS, fasta)


def test_build_fails_on_ref_mismatch(fasta):
    with pytest.raises(ValueError, match="chr1:2 .ENST1.: REF G but the FASTA has C"):
        tcga.build(kim_table().with_columns(HGVSc=pl.lit("c.10G>T")), TRANSCRIPTS, fasta)


@pytest.mark.parametrize("dataset", [tcga, mmrf_target])
def test_build_nmdeff_tables(dataset, fasta):
    expected_df = benchmark(
        {
            "chrom": ["chr1"] * 4,
            "start": [1, 2, 3, 4],
            "end": [2, 3, 4, 5],
            "ref": ["C", "G", "T", "A"],
            "alt": ["T", "T", "C", "G"],
            "transcript": ["ENST1", "ENST2", "ENST3", "ENST3"],
            "gene": ["ENSG1", "ENSG2", None, None],
            "NMD_efficiency": [1.5, 0.5, 2.0, -0.25],
            "source_row": [0, 1, 2, 3],
        }
    )
    assert_frame_equal(dataset.build(kim_table(), TRANSCRIPTS, fasta), expected_df)


def test_build_gtex(fasta, caplog):
    caplog.set_level(logging.INFO)
    raw_df = kim_table(Tissue_type=["LCL", "WHLBLD", "LCL", "WHLBLD"], ASE_altered=[0.25, 0.0, None, 0.5])
    expected_df = benchmark(
        {
            "chrom": ["chr1"] * 2,
            "start": [1, 4],
            "end": [2, 5],
            "ref": ["C", "A"],
            "alt": ["T", "G"],
            "transcript": ["ENST1", "ENST3"],
            "gene": ["ENSG1", None],
            "NMD_efficiency": [1.5, -0.25],
            "tissue": ["Cells - EBV-transformed lymphocytes", "Whole Blood"],
            "source_row": [0, 3],
        },
        {"tissue": pl.String},
    )
    assert_frame_equal(gtex.build(raw_df, TRANSCRIPTS, fasta), expected_df)
    assert "GTEx: 2 of 4 rows lack ASE_altered or have 0" in caplog.text


def test_build_gtex_rejects_an_unknown_tissue(fasta):
    raw_df = kim_table(Tissue_type=["LCL", "WHLBLD", "LCL", "BLOOD"], ASE_altered=[0.25, 0.5, 0.5, 0.5])
    with pytest.raises(pl.exceptions.InvalidOperationError, match="replace_strict"):
        gtex.build(raw_df, TRANSCRIPTS, fasta)


def test_build_geuvadis(fasta, caplog):
    caplog.set_level(logging.INFO)
    nonsense_df = pl.DataFrame(
        {
            "Variant_ID": [f"v{i}" for i in range(5)],
            "CHROM": ["chr1"] * 5,
            "POS": [1, 2, 3, 4, 8],
            "REF": ["A", "C", "G", "T", "T"],
            "ALT": ["T"] * 4 + ["A"],
            "TranscriptID": ["ENST1", "ENST1", "ENST2", "ENST2", "ENST3"],
            "VAF_mean": [0.25, None, 0.5, 0.125, 0.0],
        }
    )
    # v0 lacks flags, as a variant that the study did not evaluate
    isoforms_df = pl.DataFrame(
        {
            "Variant_ID": ["v2", "v3", "v4"],
            "resMNV": [False, True, False],
            "nonTranslating": [False, False, False],
            "Multiisoforms": [False, None, False],
        }
    )
    expected_df = benchmark(
        {
            "chrom": ["chr1"] * 2,
            "start": [0, 2],
            "end": [1, 3],
            "ref": ["A", "G"],
            "alt": ["T", "T"],
            "transcript": ["ENST1", "ENST2"],
            "gene": ["ENSG1", "ENSG2"],
            "NMD_efficiency": [1.0, 0.0],
            "source_row": [0, 2],
        }
    )
    assert_frame_equal(geuvadis.build((nonsense_df, isoforms_df), TRANSCRIPTS, fasta), expected_df)
    assert "Geuvadis: 2 of 5 rows lack VAF_mean or have 0" in caplog.text
    assert "Geuvadis: 1 of 3 rows are outside the study's accurate annotation" in caplog.text


def test_check_sha256():
    data = b"NMD"
    assert check_sha256(data, hashlib.sha256(data).hexdigest(), "data") == data
    with pytest.raises(ValueError, match="data has sha256"):
        check_sha256(data, hashlib.sha256(b"other").hexdigest(), "data")


def test_fetch_reads_downloads(tmp_path):
    """fetch reads a file or a zip member that downloads has, and checks its sha256."""
    (tmp_path / "table.csv").write_bytes(b"a,b\n")
    (tmp_path / "member.txt").write_bytes(b"c\n")
    sha256 = hashlib.sha256(b"a,b\n").hexdigest()
    assert fetch(tmp_path, "https://example.org/files/table.csv", sha256) == b"a,b\n"
    member = fetch(
        tmp_path, "https://example.org/files/archive.zip", hashlib.sha256(b"c\n").hexdigest(), "dir/member.txt"
    )
    assert member == b"c\n"
    with pytest.raises(ValueError, match="table.csv has sha256"):
        fetch(tmp_path, "https://example.org/files/table.csv", hashlib.sha256(b"other").hexdigest())


def test_fetch_downloads_a_file_and_a_zip_member(tmp_path):
    """fetch downloads with fsspec and saves in downloads. The memory filesystem stands in for HTTP."""
    memory = fsspec.filesystem("memory")
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("dir/other.txt", b"x" * 1000)
        archive.writestr("dir/member.txt", b"c\n")
    memory.pipe("/files/table.csv", b"a,b\n")
    memory.pipe("/files/archive.zip", buffer.getvalue())
    sha256 = hashlib.sha256(b"a,b\n").hexdigest()
    assert fetch(tmp_path / "new", "memory://files/table.csv", sha256) == b"a,b\n"
    assert (tmp_path / "new" / "table.csv").read_bytes() == b"a,b\n"
    member = fetch(tmp_path / "new", "memory://files/archive.zip", hashlib.sha256(b"c\n").hexdigest(), "dir/member.txt")
    assert member == b"c\n"
    assert (tmp_path / "new" / "member.txt").read_bytes() == b"c\n"
    with pytest.raises(ValueError, match="memory://files/table.csv has sha256"):
        fetch(tmp_path / "other", "memory://files/table.csv", hashlib.sha256(b"other").hexdigest())


def test_datasets_cover_every_module():
    modules = {module.name for module in pkgutil.iter_modules(datasets.__path__)}
    expected = {name: f"{datasets.__name__}.{name}" for name in modules}
    assert {name: dataset.__name__ for name, dataset in DATASETS.items()} == expected
    for name, dataset in DATASETS.items():
        assert isinstance(dataset.NAME, str), name
        assert callable(dataset.download), name
        assert callable(dataset.build), name


def test_join_features(caplog):
    caplog.set_level(logging.INFO)
    table_df = pl.DataFrame(
        {
            "transcript": ["ENST1", "ENST1", "ENST2", "ENST3"],
            "chrom": ["chr1"] * 4,
            "start": [0, 1, 2, 3],
            "end": [1, 2, 3, 4],
            "ref": ["A"] * 4,
            "alt": ["T"] * 4,
            "gene": ["ENSG1", "ENSG1", "ENSG2", None],
        }
    )
    # The result of annotate_variants has transcript_id with version
    results_df = pd.DataFrame(
        {
            "transcript_id": ["ENST1.1", "ENST1.1"],
            "chrom": ["chr1"] * 2,
            "start": [0, 1],
            "end": [1, 2],
            "ref": ["A"] * 2,
            "alt": ["T"] * 2,
            "nmd_efficiency_prediction_status": ["ok", "ref_ptc"],
        }
    )
    kept_df = join_features(table_df, results_df, "Test")
    assert kept_df.select("transcript", "start", "nmd_efficiency_prediction_status").rows() == [("ENST1", 0, "ok")]
    assert (
        "Test: 2 of 4 rows have no NMD-Scanner row of their transcript and variant: 1 on transcripts missing from the"
        " GFF3, 1 with a variant that touches no coding region of its transcript (ENST2)"
    ) in caplog.text
    assert (
        "Test: 1 of the 2 joined rows have an nmd_efficiency_prediction_status other than ok: 1 ref_ptc" in caplog.text
    )
