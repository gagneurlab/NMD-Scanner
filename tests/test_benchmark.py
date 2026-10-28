"""Tests of nmd_scanner.benchmark on tiny tables, without network."""

import hashlib
import io
import logging
import pkgutil
import zipfile

import polars as pl
import pytest
from polars.testing import assert_frame_equal
from pyfaidx import Fasta

from nmd_scanner.benchmark import DATASETS, datasets, join_features
from nmd_scanner.benchmark.datasets import geuvadis, gtex, mmrf_target, tcga
from nmd_scanner.benchmark.sources import RemoteFile, check_sha256, fetch, hgvsc_alleles, read_transcripts

# ENST3 is missing from the GFF3
TRANSCRIPTS = pl.DataFrame({"transcript": ["ENST1", "ENST2"], "strand": ["+", "-"], "gene": ["ENSG1", "ENSG2"]})
SCHEMA = {
    "chrom": pl.String,
    "pos": pl.Int64,
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
            "pos": [2, 3, 4, 5],
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
    raw_df = kim_table(Tissue_type=["LIVER", "LUNG", "LIVER", "LUNG"], ASE_altered=[0.25, 0.0, None, 0.5])
    expected_df = benchmark(
        {
            "chrom": ["chr1"] * 2,
            "pos": [2, 5],
            "ref": ["C", "A"],
            "alt": ["T", "G"],
            "transcript": ["ENST1", "ENST3"],
            "gene": ["ENSG1", None],
            "NMD_efficiency": [1.5, -0.25],
            "tissue": ["LIVER", "LUNG"],
            "source_row": [0, 3],
        },
        {"tissue": pl.String},
    )
    assert_frame_equal(gtex.build(raw_df, TRANSCRIPTS, fasta), expected_df)
    assert "GTEx: 2 of 4 rows lack ASE_altered or have 0" in caplog.text


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
            "pos": [1, 3],
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


def test_read_transcripts(tmp_path):
    gff3 = tmp_path / "annotation.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "chr1\tHAVANA\tgene\t1\t8\t.\t-\t.\tID=ENSG2.1;gene_id=ENSG2.1\n"
        "chr1\tHAVANA\ttranscript\t1\t8\t.\t-\t.\tID=ENST2.4;Parent=ENSG2.1\n"
        "chr1\tHAVANA\texon\t1\t8\t.\t-\t.\tID=exon:ENST2.4:1;Parent=ENST2.4\n"
    )
    assert read_transcripts(gff3).rows() == [("ENST2", "-", "ENSG2")]


def test_read_transcripts_ensembl(tmp_path):
    gff3 = tmp_path / "annotation.gff3"
    gff3.write_text(
        "##gff-version 3\n"
        "chr1\tensembl\tgene\t1\t8\t.\t-\t.\tID=gene:ENSG2;biotype=protein_coding\n"
        "chr1\tensembl\tmRNA\t1\t8\t.\t-\t.\tID=transcript:ENST2;Parent=gene:ENSG2\n"
        "chr1\tensembl\texon\t1\t8\t.\t-\t.\tID=ENSE2;Parent=transcript:ENST2\n"
        "chr1\tensembl\tCDS\t1\t8\t.\t-\t0\tID=CDS:ENSP2;Parent=transcript:ENST2\n"
    )
    assert read_transcripts(gff3).rows() == [("ENST2", "-", "ENSG2")]


def test_read_transcripts_fails_without_enst_ids(tmp_path):
    gff3 = tmp_path / "annotation.gff3"
    gff3.write_text("chr1\tensembl\tgene\t1\t8\t.\t-\t.\tID=gene:ENSG2.1\n")
    with pytest.raises(ValueError, match="has no transcript line with an ID that starts with ENST"):
        read_transcripts(gff3)


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
            "pos": [1, 2, 3, 4],
            "ref": ["A"] * 4,
            "alt": ["T"] * 4,
            "gene": ["ENSG1", "ENSG1", "ENSG2", None],
        }
    )
    features_df = pl.DataFrame(
        {
            "transcript": ["ENST1", "ENST1"],
            "chrom": ["chr1"] * 2,
            "pos": [1, 2],
            "ref": ["A"] * 2,
            "alt": ["T"] * 2,
            "nmd_efficiency_prediction_status": ["ok", "ref_ptc"],
        }
    )
    kept_df = join_features(table_df, features_df, "Test")
    assert kept_df.select("transcript", "pos", "nmd_efficiency_prediction_status").rows() == [("ENST1", 1, "ok")]
    assert (
        "Test: 2 of 4 rows have no NMD-Scanner row of their transcript and variant: 1 on transcripts missing from the"
        " GFF3, 1 with a variant that touches no coding region of its transcript (ENST2)"
    ) in caplog.text
    assert (
        "Test: 1 of the 2 joined rows have an nmd_efficiency_prediction_status other than ok: 1 ref_ptc" in caplog.text
    )


class FakeResponse(io.BytesIO):
    def __init__(self, data: bytes, status: int) -> None:
        super().__init__(data)
        self.status = status


def serve(monkeypatch, data: bytes, status: int = 206) -> list[str]:
    """Make urlopen answer with the byte range of data that a request asks for, and return the ranges asked."""
    ranges = []

    def urlopen(request, timeout):
        ranges.append(request.headers["Range"])
        first, last = map(int, request.headers["Range"].removeprefix("bytes=").split("-"))
        return FakeResponse(data[first : last + 1], status)

    monkeypatch.setattr("urllib.request.urlopen", urlopen)
    return ranges


def test_fetch_downloads_a_zip_member(tmp_path, monkeypatch):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("dir/other.txt", b"x" * 1000)
        archive.writestr("dir/member.txt", b"c\n")
    data = buffer.getvalue()
    ranges = serve(monkeypatch, data)
    sha256 = hashlib.sha256(b"c\n").hexdigest()
    url = "https://example.org/files/archive.zip"
    assert fetch(tmp_path, url, sha256, "dir/member.txt", len(data)) == b"c\n"
    assert (tmp_path / "member.txt").read_bytes() == b"c\n"
    assert len(ranges) > 0
    assert all(int(r.removeprefix("bytes=").split("-")[1]) < len(data) for r in ranges)
    with pytest.raises(ValueError, match="needs the size of the zip"):
        fetch(tmp_path / "other", url, sha256, "dir/member.txt")


def test_remote_file_fails_if_the_server_ignores_the_range(monkeypatch):
    serve(monkeypatch, b"abcdef", status=200)
    with pytest.raises(ValueError, match="answered a range request with HTTP 200"):
        RemoteFile("https://example.org/files/archive.zip", 6).read(3)


def test_fetch_downloads_a_file(tmp_path, monkeypatch):
    monkeypatch.setattr("urllib.request.urlopen", lambda url, timeout: io.BytesIO(b"a,b\n"))
    sha256 = hashlib.sha256(b"a,b\n").hexdigest()
    assert fetch(tmp_path / "new", "https://example.org/files/table.csv", sha256) == b"a,b\n"
    assert (tmp_path / "new" / "table.csv").read_bytes() == b"a,b\n"
