"""
GTEx: germline nonsense variants in the tissues of GTEx v8, Supplementary Data 2 of Kim et al. 2024 (Commun Biol
7:1461).

The table has one row per variant and tissue. NMD_efficiency is its column NMD_efficiency, and ASE_altered is the
variant allele fraction in the RNA. For ASE_altered 0, the table caps NMD_efficiency at 10, so build keeps the rows with
ASE_altered above 0. In 4,102 of the 17,816 rows with ASE_altered above 0, NMD_efficiency differs from
-log2(ASE_altered / 0.5). The difference moves the Spearman correlations of scripts/validate_model.py by only about
0.001.
"""

import logging
from pathlib import Path

import polars as pl
from pyfaidx import Fasta

from nmd_scanner.benchmark.sources import fetch, finish_table, hgvsc_alleles

logger = logging.getLogger(__name__)

NAME = "GTEx"
# Supplementary Data 1 to 12 of Kim et al. 2024, one sheet each
URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-024-07136-y/MediaObjects/"
    "42003_2024_7136_MOESM3_ESM.xlsx"
)
SHA256 = "637c92b17095ef6147b74ae0ab6549f674adc6a67ee04c06d483fa7331ff34c6"
SHEET = "Supplementary Data 2"


def download(downloads: Path) -> pl.DataFrame:
    """Read the sheet of Supplementary Data 2 with fetch. polars reads it with fastexcel, of the extra scripts."""
    # The first line of the sheet is its title
    return pl.read_excel(fetch(downloads, URL, SHA256), sheet_name=SHEET, read_options={"header_row": 1})


def build(raw_df: pl.DataFrame, transcripts_df: pl.DataFrame, fasta: Fasta) -> pl.DataFrame:
    """
    The benchmark table of raw_df, the table of download, with the alleles of hgvsc_alleles and the column tissue. It has
    the rows with ASE_altered above 0.
    """
    table_df = raw_df.with_row_index("source_row").filter(pl.col("ASE_altered") > 0)
    logger.info("%s: %d of %d rows lack ASE_altered or have 0", NAME, len(raw_df) - len(table_df), len(raw_df))
    table_df = table_df.rename({"Tissue_type": "tissue"})
    return finish_table(hgvsc_alleles(table_df, transcripts_df, fasta), transcripts_df, fasta)
