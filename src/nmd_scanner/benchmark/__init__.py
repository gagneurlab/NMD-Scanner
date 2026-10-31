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
  scripts/build_benchmarks.py builds it. build applies the study's own filters, and it logs the measurements that they
  drop.

Every benchmark table has these columns, with one row per measurement:
- chrom, start, end, ref and alt: the variant on GRCh38, with the names and the coordinates of the NMD-Scanner
  output. start and end are 0-based half-open, and ref and alt are on the forward strand.
- transcript: the Ensembl ID of the measured transcript, as the study gives it.
- gene: the Ensembl ID without version of the gene of the transcript in the GFF3. It is null for a transcript missing
  from the GFF3.
- NMD_efficiency: the NMD efficiency that the study measured. 0 means no NMD, and 1 means that the variant allele
  fraction in the RNA is half of that in the DNA.
- tissue: the GTEx v8 name (SMTSD) of the tissue of the measurement, only in a table with one row per variant and
  tissue.
- source_row: the index of the row in the study table, from 0. Sorting on it gives the rows a fixed order.

nmd_scanner.annotate_variants takes the unique VARIANT_COLUMNS of a benchmark table as its vcf_df, in pandas.
join_features joins the result of annotate_variants to the benchmark table, and it keeps the measurements that the NMD
efficiency model can predict. Read the annotation_df of annotate_variants from the GFF3 of the benchmark tables, so
that a null gene marks a transcript missing from it.
"""

import logging

import pandas as pd
import polars as pl

import nmd_scanner
from nmd_scanner.benchmark.datasets import geuvadis, gtex, mmrf_target, tcga

logger = logging.getLogger(__name__)

DATASETS = {"tcga": tcga, "geuvadis": geuvadis, "gtex": gtex, "mmrf_target": mmrf_target}
# The columns of the variant of a measurement
VARIANT_COLUMNS = ["chrom", "start", "end", "ref", "alt"]
# The columns that join a measurement to the NMD-Scanner result of its variant-transcript pair
KEY = ["transcript", *VARIANT_COLUMNS]


def join_features(table_df: pl.DataFrame, results_df: pd.DataFrame, dataset: str) -> pl.DataFrame:
    """
    Join each measurement of a benchmark table to its variant-transcript pair in results_df, the result of
    nmd_scanner.annotate_variants.
    Keep the measurements with nmd_efficiency_prediction_status "ok". Log the measurements that this drops, and why,
    with dataset as the name of the table. The rows of the result have no fixed order.

    A measurement has no variant-transcript pair in results_df if the GFF3 lacks its transcript, or if its variant
    touches no coding region of the transcript (see "Technical Notes.md").
    """
    # transcript is transcript_id without version
    features_df = pl.DataFrame(nmd_scanner.to_arrow(results_df)).with_columns(
        transcript=pl.col("transcript_id").str.split(".").list.first()
    )
    joined_df = table_df.join(features_df, on=KEY, validate="m:1")
    unmatched_df = table_df.join(features_df, on=KEY, how="anti")
    missing = unmatched_df["gene"].is_null()
    if len(unmatched_df) > 0:
        logger.info(
            "%s: %d of %d measurements have no NMD-Scanner result for their variant-transcript pair: %d on transcripts"
            " missing from the GFF3, %d with a variant that touches no coding region of its transcript (%s)",
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
            "%s: %d of the %d joined measurements have an nmd_efficiency_prediction_status other than ok: %s",
            dataset,
            dropped_df["count"].sum(),
            len(joined_df),
            ", ".join(f"{count} {value}" for value, count in dropped_df.iter_rows()),
        )
    return joined_df.filter(status == "ok")
