"""
Geuvadis: germline nonsense SNVs in lymphoblastoid cell lines (Iha et al. 2025), Zenodo record 16666299.

Nonsense_list.txt has one row per variant and MANE Select transcript, with REF and ALT on the forward strand. VAF_mean
is the variant allele fraction in the RNA. NMD_efficiency is -log2(VAF_mean / 0.5), because a heterozygous germline
variant has a variant allele fraction of 0.5 in the DNA. build keeps the rows with VAF_mean above 0 in the study's
accurate annotation: no rescued MNV, no nontranslating gene, no mixed NMD isoforms.
"""

import logging
from pathlib import Path

import polars as pl
from pyfaidx import Fasta

from nmd_scanner.benchmark.sources import fetch, finish_table

logger = logging.getLogger(__name__)

NAME = "Geuvadis"
# Zenodo record 16666299 has one file, NMD-rules.zip, of ZENODO_SIZE bytes
ZENODO_URL = "https://zenodo.org/api/records/16666299/files/NMD-rules.zip/content"
ZENODO_SIZE = 3_189_701_188
NONSENSE_LIST = "NMD-rules/Datasets/Nonsense_list.txt"
ISOFORMS = "NMD-rules/Datasets/IsoformLevel.NMDannotation.txt"
SHA256 = {
    NONSENSE_LIST: "5a248be88a6b2a287c2baf12e542d638257a06daadc21a4c440af812518bce26",
    ISOFORMS: "bc16f2b869149f07e0a5d7ce4173df4fbca9ea9c4fab2dd9f3c07a82c4011a7d",
}
# A variant with one of these flags in IsoformLevel.NMDannotation.txt is not in the accurate annotation
ACCURATE_ANNOTATION_FLAGS = ["resMNV", "nonTranslating", "Multiisoforms"]
# The variant allele fraction in the DNA of a heterozygous germline variant
GERMLINE_VAF_DNA = 0.5


def download(downloads: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Read Nonsense_list.txt and IsoformLevel.NMDannotation.txt of NMD-rules.zip with fetch. The record is one zip of
    3.2 GB, so fetch reads just these 2 members.
    """
    nonsense, isoforms = (
        fetch(downloads, ZENODO_URL, SHA256[member], member, ZENODO_SIZE) for member in (NONSENSE_LIST, ISOFORMS)
    )
    return (
        pl.read_csv(nonsense, separator="\t", infer_schema_length=None),
        pl.read_csv(isoforms, separator="\t", columns=["Variant_ID", *ACCURATE_ANNOTATION_FLAGS]),
    )


def build(raw: tuple[pl.DataFrame, pl.DataFrame], transcripts_df: pl.DataFrame, fasta: Fasta) -> pl.DataFrame:
    """
    The benchmark table of raw, the tables Nonsense_list.txt and IsoformLevel.NMDannotation.txt of download. It has
    the rows with VAF_mean above 0 in the accurate annotation. IsoformLevel.NMDannotation.txt has the evaluated
    variants, and a variant that it lacks has no flag.
    """
    nonsense_df, isoforms_df = raw
    table_df = nonsense_df.with_row_index("source_row").filter(pl.col("VAF_mean") > 0)
    logger.info("%s: %d of %d rows lack VAF_mean or have 0", NAME, len(nonsense_df) - len(table_df), len(nonsense_df))
    flagged = isoforms_df.filter(pl.any_horizontal(ACCURATE_ANNOTATION_FLAGS))["Variant_ID"]
    accurate_df = table_df.filter(~pl.col("Variant_ID").is_in(flagged.implode()))
    logger.info(
        "%s: %d of %d rows are outside the study's accurate annotation (flags %s)",
        NAME,
        len(table_df) - len(accurate_df),
        len(table_df),
        " ".join(ACCURATE_ANNOTATION_FLAGS),
    )
    table_df = accurate_df.rename(
        {"CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt", "TranscriptID": "transcript"}
    )
    return finish_table(
        table_df.with_columns(NMD_efficiency=-(pl.col("VAF_mean") / GERMLINE_VAF_DNA).log(2)), transcripts_df, fasta
    )
