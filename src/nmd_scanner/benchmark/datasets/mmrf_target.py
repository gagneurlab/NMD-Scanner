"""
MMRF-TARGET: somatic nonsense variants in tumors of the projects MMRF and TARGET, the table MMRF_TARGET_dataset.csv
of NMDEff (Kim et al. 2024).

NMD_efficiency is the column of the table. Every row has VAF_RNA above 0.
"""

from pathlib import Path

import polars as pl
from pyfaidx import Fasta

from nmd_scanner.benchmark.sources import NMDEFF_URL, fetch, hgvsc_alleles, to_schema

NAME = "MMRF-TARGET"
SHA256 = "c063914da054292c6c6b867eb58650d5378c74d5c852f54c874eab8e51cb4a19"


def download(downloads: Path) -> pl.DataFrame:
    """Read MMRF_TARGET_dataset.csv of NMDEff with fetch."""
    return pl.read_csv(fetch(downloads, NMDEFF_URL + "MMRF_TARGET_dataset.csv", SHA256), infer_schema_length=None)


def build(raw: pl.DataFrame, transcripts: pl.DataFrame, fasta: Fasta) -> pl.DataFrame:
    """The benchmark table of raw, the table of download, with the alleles of hgvsc_alleles."""
    return to_schema(hgvsc_alleles(raw.with_row_index("source_row"), transcripts, fasta), transcripts, fasta)
