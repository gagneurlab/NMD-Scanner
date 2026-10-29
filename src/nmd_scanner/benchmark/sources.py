"""
The sources of the benchmarks, and the steps that the dataset modules share to turn a study table into a benchmark
table (see nmd_scanner.benchmark).
"""

import hashlib
import logging
import urllib.parse
from pathlib import Path

import polars as pl
from pyfaidx import Fasta

logger = logging.getLogger(__name__)

# Both tables of NMDEff (https://github.com/hjkng/nmdeff) were last changed in this commit, which is also its tag v1.0
NMDEFF_URL = "https://raw.githubusercontent.com/hjkng/nmdeff/08c92768fcb689236a833db6a2f2d9bcbe919f12/"
# The columns of every benchmark table, in order. tissue, if a table has it, and source_row follow them.
COLUMNS = ["chrom", "start", "end", "ref", "alt", "transcript", "gene", "NMD_efficiency"]
HGVSC_SUBSTITUTION = r"^(?:[^:]+:)?c\.\d+(?P<ref>[ACGT])>(?P<alt>[ACGT])$"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def check_sha256(data: bytes, expected: str, name: str) -> bytes:
    """Return data if its sha256 is expected, and fail otherwise. name names data in the error."""
    digest = hashlib.sha256(data).hexdigest()
    if digest != expected:
        raise ValueError(f"{name} has sha256 {digest}, expected {expected}")
    return data


def fetch(downloads: Path, url: str, sha256: str, member: str | None = None) -> bytes:
    """
    Return the file at url, or with member, the member of the zip at url. Check its sha256.

    The file is read from the directory downloads if it is there, and otherwise downloaded with fsspec (extra scripts)
    and saved there. Its name in downloads is the last part of url, or of member. fsspec reads a member through HTTP
    range requests, so only its bytes are downloaded.
    """
    path = downloads / Path(member or urllib.parse.urlparse(url).path).name
    if path.exists():
        return check_sha256(path.read_bytes(), sha256, str(path))
    import fsspec

    with fsspec.open(f"zip://{member}::{url}" if member else url) as file:
        data = file.read()
    check_sha256(data, sha256, url)
    downloads.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return data


def hgvsc_alleles(table_df: pl.DataFrame, transcripts_df: pl.DataFrame, fasta: Fasta) -> pl.DataFrame:
    """
    Rename the columns chromosome, start and Transcript_ID of a table of Kim et al. 2024 to chrom, pos and transcript,
    and add ref and alt on the forward strand. transcripts_df is the transcript table of build (see
    nmd_scanner.benchmark).

    HGVSc (c.123G>T, optionally prefixed by the transcript ID) has the alleles on the transcript strand. On a
    minus-strand transcript of the GFF3, complement them. For a transcript missing from the GFF3, use the orientation
    whose REF matches the FASTA. Fail if HGVSc is not a single-base substitution. The transcript IDs of the table must
    lack the version, like those of transcripts_df.
    """
    table_df = table_df.rename({"chromosome": "chrom", "start": "pos", "Transcript_ID": "transcript"})
    table_df = table_df.with_columns(
        pl.col("HGVSc").str.extract_groups(HGVSC_SUBSTITUTION).struct.unnest(),
        genomic_ref=pl.Series(
            [fasta[chrom][pos - 1 : pos].seq.upper() for chrom, pos in zip(table_df["chrom"], table_df["pos"])],
            dtype=pl.String,
        ),
        strand=pl.col("transcript").replace_strict(
            transcripts_df["transcript"], transcripts_df["strand"], default=None
        ),
    )
    if (invalid_df := table_df.filter(pl.col("ref").is_null())).height > 0:
        row = invalid_df.row(0, named=True)
        raise ValueError(f"{row['HGVSc']} ({row['transcript']}) is not a single-base substitution")
    strand = pl.col("strand")
    complement = (strand == "-") | (strand.is_null() & (pl.col("ref") != pl.col("genomic_ref")))
    return table_df.with_columns(
        pl.when(complement).then(pl.col("ref", "alt").replace(COMPLEMENT)).otherwise(pl.col("ref", "alt"))
    ).drop("genomic_ref", "strand")


def finish_table(table_df: pl.DataFrame, transcripts_df: pl.DataFrame, fasta: Fasta) -> pl.DataFrame:
    """
    Add the 0-based half-open start and end of the 1-based pos of table_df, and the gene of each transcript. Check each
    REF against the FASTA, and keep the benchmark columns: COLUMNS, tissue if table_df has it, and source_row. gene
    comes from transcripts_df (see nmd_scanner.benchmark), and it is null for a transcript that the GFF3 lacks. Raise a
    ValueError at the first REF that differs from the FASTA.
    """
    table_df = table_df.with_columns(start=pl.col("pos").cast(pl.Int64) - 1).with_columns(
        end=pl.col("start") + pl.col("ref").str.len_chars().cast(pl.Int64)
    )
    table_df = table_df.with_columns(
        fasta_ref=pl.Series(
            [
                fasta[chrom][start:end].seq.upper()
                for chrom, start, end in table_df.select("chrom", "start", "end").iter_rows()
            ],
            dtype=pl.String,
        ),
        gene=pl.col("transcript").replace_strict(transcripts_df["transcript"], transcripts_df["gene"], default=None),
    )
    if (mismatches_df := table_df.filter(pl.col("ref") != pl.col("fasta_ref"))).height > 0:
        row = mismatches_df.row(0, named=True)
        raise ValueError(
            f"{row['chrom']}:{row['pos']} ({row['transcript']}): REF {row['ref']} but the FASTA has {row['fasta_ref']}"
        )
    return table_df.select(*COLUMNS, *(["tissue"] if "tissue" in table_df.columns else []), "source_row")
