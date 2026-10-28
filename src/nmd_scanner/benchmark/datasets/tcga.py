"""
TCGA: somatic nonsense variants in tumors of TCGA, the table tcga_dataset.csv of NMDEff (Kim et al. 2024).

NMD_efficiency is the column of the table. In 29 rows, it differs from -log2(VAF_DNA_RNA_ratio) by up to 0.36. Every
row has VAF_RNA above 0. scripts/train_model.py trains the NMD efficiency model on this table.
"""

from pathlib import Path

import polars as pl
from pyfaidx import Fasta

from nmd_scanner.benchmark.sources import NMDEFF_URL, fetch, finish_table, hgvsc_alleles

NAME = "TCGA"
SHA256 = "7458af15a699a643de3638114af99a17d7eb95b3f80bce187e45347b9c1dc3ca"


def download(downloads: Path) -> pl.DataFrame:
    """Read tcga_dataset.csv of NMDEff with fetch."""
    return pl.read_csv(fetch(downloads, NMDEFF_URL + "tcga_dataset.csv", SHA256), infer_schema_length=None)


def build(raw_df: pl.DataFrame, transcripts_df: pl.DataFrame, fasta: Fasta) -> pl.DataFrame:
    """The benchmark table of raw_df, the table of download, with the alleles of hgvsc_alleles."""
    return finish_table(
        hgvsc_alleles(raw_df.with_row_index("source_row"), transcripts_df, fasta), transcripts_df, fasta
    )
