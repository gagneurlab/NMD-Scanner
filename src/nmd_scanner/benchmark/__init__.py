"""
The NMD efficiency benchmarks: tables of the NMD efficiency that public studies measured for nonsense variants.
scripts/build_benchmarks.py builds them. scripts/train_model.py trains the NMD efficiency model on TCGA, and
scripts/validate_model.py.py tests it on all of them.

DATASETS maps the name of each benchmark to its module in nmd_scanner.benchmark.datasets. A module has NAME, the name
of the benchmark in logs and figures, and two functions:
- download(downloads) returns the raw tables of the study. It downloads them with sources.fetch, which keeps the files
  in the directory downloads.
- build(raw, transcripts_df, fasta) returns the benchmark table of the raw tables, without network. fasta is a
  pyfaidx.Fasta of GRCh38. transcripts_df has one row per transcript of nmd_scanner.read_annotation, with the columns
  transcript (transcript_id without version), strand (Strand) and gene (gene_id without version).
  scripts/build_benchmarks.py builds it. build applies the study's own filters, and it logs the rows that they drop.

Every benchmark table has these columns, with one row per measurement:
- chrom, pos, ref and alt: the variant on GRCh38. pos is 1-based, and ref and alt are on the forward strand.
- transcript: the Ensembl ID of the measured transcript, as the study gives it.
- gene: the Ensembl ID without version of the gene of the transcript in the GFF3. It is null for a transcript missing
  from the GFF3.
- NMD_efficiency: the NMD efficiency that the study measured. 0 means no NMD, and 1 means that the variant allele
  fraction in the RNA is half of that in the DNA.
- tissue: the tissue of the measurement, only in a table with one row per variant and tissue.
- source_row: the index of the row in the study table, from 0. Sorting on it gives the rows a fixed order.

run_nmd_scanner runs NMD-Scanner on the variants of benchmark tables, and join_features keeps the rows that the NMD
efficiency model can predict. Read the annotation_df of run_nmd_scanner from the GFF3 of the benchmark tables, so
that a null gene marks a transcript missing from it.
"""

import logging

import pandas as pd
import polars as pl
from pyfaidx import Fasta

import nmd_scanner
from nmd_scanner.benchmark.datasets import geuvadis, gtex, mmrf_target, tcga

logger = logging.getLogger(__name__)

DATASETS = {"tcga": tcga, "geuvadis": geuvadis, "gtex": gtex, "mmrf_target": mmrf_target}
# The columns that join a row of a benchmark table to the NMD-Scanner row of its transcript and variant
KEY = ["transcript", "chrom", "pos", "ref", "alt"]


def run_nmd_scanner(
    variants_df: pl.DataFrame, annotation_df: pd.DataFrame, fasta: Fasta, prediction: bool = True
) -> pl.DataFrame:
    """
    Run nmd_scanner.annotate_variants on the unique variants of the columns chrom, pos, ref and alt of variants_df,
    without the sequence columns. annotation_df is the table of nmd_scanner.read_annotation, and fasta is the
    pyfaidx.Fasta that it took. Return the result as a polars table with the columns of KEY: transcript is
    transcript_id without version, and pos is the 1-based start.
    """
    # vcf_df is the table of nmd_scanner.read_vcf, whose Start and End are 0-based half-open
    start = pl.col("pos").cast(pl.Int64) - 1
    vcf_df = (
        variants_df.select("chrom", "pos", "ref", "alt")
        .unique()
        .sort("chrom", "pos", "ref", "alt")
        .select(
            Chromosome="chrom",
            Start=start,
            End=start + pl.col("ref").str.len_chars().cast(pl.Int64),
            ID=pl.lit("."),
            Ref="ref",
            Alt="alt",
        )
        .to_pandas()
    )
    result = nmd_scanner.annotate_variants(vcf_df, annotation_df, fasta, sequences=False, prediction=prediction)
    # NMD-Scanner's start is 0-based
    return pl.DataFrame(nmd_scanner.to_arrow(result)).with_columns(
        transcript=pl.col("transcript_id").str.split(".").list.first(), pos=pl.col("start") + 1
    )


def join_features(table_df: pl.DataFrame, features_df: pl.DataFrame, dataset: str) -> pl.DataFrame:
    """
    Join each row of a benchmark table to the row of features_df, the result of run_nmd_scanner, of its transcript and
    variant.
    Keep the rows with nmd_efficiency_prediction_status "ok". Log the rows that this drops, and why, with dataset as
    the name of the table. The rows have no fixed order.

    A row has no NMD-Scanner row if the GFF3 lacks its transcript, or if its variant touches no coding region of the
    transcript (see "Technical Notes.md").
    """
    joined_df = table_df.join(features_df, on=KEY, validate="m:1")
    unmatched_df = table_df.join(features_df, on=KEY, how="anti")
    missing = unmatched_df["gene"].is_null()
    if len(unmatched_df) > 0:
        logger.info(
            "%s: %d of %d rows have no NMD-Scanner row of their transcript and variant: %d on transcripts missing from"
            " the GFF3, %d with a variant that touches no coding region of its transcript (%s)",
            dataset,
            len(unmatched_df),
            len(table_df),
            missing.sum(),
            (~missing).sum(),
            " ".join(unmatched_df.filter(~missing)["transcript"].unique().sort()) or "none",
        )
    status = joined_df["nmd_efficiency_prediction_status"]
    dropped_df = status.filter(status != "ok").value_counts(sort=True)
    if len(dropped_df) > 0:
        logger.info(
            "%s: %d of the %d joined rows have an nmd_efficiency_prediction_status other than ok: %s",
            dataset,
            dropped_df["count"].sum(),
            len(joined_df),
            ", ".join(f"{count} {value}" for value, count in dropped_df.iter_rows()),
        )
    return joined_df.filter(status == "ok")
