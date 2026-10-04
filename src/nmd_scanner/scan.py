# Import dependencies
import logging
import os

import pandas as pd
import pyranges as pr
from pyfaidx import Fasta

# Create the functions used for reading in the files (VCF, GTF, FASTA)

logger = logging.getLogger(__name__)


def read_vcf(vcf_path):
    # TODO: adjust this function to also get structural variants and then also adjust downstream analysis
    #  (especially the inclusion of the Variants into the reference CDS sequence to create the alternative CDS)

    """
    Read a single VCF file into a PyRanges object with adjusted coordinates.

    The VCF must be left-normalized and single-allelic (one ALT allele per
    record), e.g. produced by ``bcftools norm -m- -f reference.fa``.
    Multi-allelic records (comma-separated ALT) are rejected because the
    downstream variant application assumes exactly one ALT allele per row.
    """
    df = pd.read_csv(
        vcf_path,
        comment="#",
        sep="\t",
        header=None,
        names=["Chromosome", "Start", "ID", "Ref", "Alt", "Qual", "Filter", "Info"],
    )

    # Reject multi-allelic records: the VCF spec allows comma-separated ALT,
    # but the rest of the pipeline assumes one ALT allele per row.
    multiallelic = df["Alt"].astype(str).str.contains(",")
    if multiallelic.any():
        n_multiallelic = int(multiallelic.sum())
        raise ValueError(
            f"{n_multiallelic} multi-allelic record(s) found in {vcf_path}. "
            "NMD-Scanner requires a left-normalized, single-allelic VCF. "
            "Split and normalize first, e.g. `bcftools norm -m- -f reference.fa`."
        )

    # Adjust coordinates to 0-based
    df["Start"] = df["Start"] - 1
    df["End"] = df["Start"] + df["Ref"].str.len()

    # Keep only relevant columns
    gr = pr.PyRanges(df[["Chromosome", "Start", "End", "ID", "Ref", "Alt", "Qual", "Filter", "Info"]])
    return gr


def read_gtf(gtf_path):
    """
    Reads a GTF file into a PyRanges object, with its rows as they are in the file. A GTF CDS excludes
    the stop codon, which has its own stop_codon rows. ``merge_stop_codons_into_cds`` builds the coding
    regions that ``rules.extract_ptc`` takes.
    """
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")
    return pr.read_gtf(gtf_path)


def merge_stop_codons_into_cds(df, transcript_col="transcript_id"):
    """
    Returns the coding regions of a GTF, in the form ``rules.extract_ptc`` takes: one CDS row per
    transcript and exon, the union of its CDS and stop_codon rows.

    A GTF CDS excludes the stop codon, which has its own stop_codon rows. A stop codon split across an
    intron has two of them, and the second one can lie in an exon without CDS. A transcript without
    stop_codon rows, e.g. one tagged cds_end_NF, keeps its CDS as it is: its coding region ends without
    a stop codon.

    The rows are merged per transcript and exon_number. To recompute the exon numbers with
    ``compute_exon_numbers``, do so before the merge.

    :param df: GTF rows (DataFrame) with Feature, Start, End, exon_number and ``transcript_col``. Only
        the CDS and stop_codon rows are used.
    :param transcript_col: The name of the column that indicates the transcript ID
    :return: DataFrame with the CDS rows, each extended by the stop codon bases of its exon, plus one
        CDS row for each exon that holds only stop codon bases. exon_number is int. The column
        has_stop_codon says whether the coding region of the transcript ends in an annotated stop codon,
        i.e. whether the transcript has stop_codon rows.
    :raises ValueError: if stop codon bases do not touch or overlap the CDS of their exon.
    """
    df = df[df["Feature"].isin(["CDS", "stop_codon"])].copy()
    df["exon_number"] = df["exon_number"].astype(int)
    keys = [transcript_col, "exon_number"]

    is_stop = df["Feature"] == "stop_codon"
    stop_rows = df[is_stop]
    if stop_rows.empty and not df.empty:
        logger.warning(
            "No stop_codon rows found next to the CDS rows: every transcript is treated as having no annotated "
            "stop codon (no 3'UTR length, no stop codon distance, every in-frame stop is premature). "
            "A GTF CDS excludes the stop codon, so the merge needs the stop_codon rows of the GTF."
        )
    stops = stop_rows.groupby(keys, observed=True).agg(stop_start=("Start", "min"), stop_end=("End", "max"))

    cds = df[~is_stop].merge(stops, left_on=keys, right_index=True, how="left")
    has_stop = cds["stop_start"].notna()
    # coordinates are half-open, so touching intervals share one coordinate
    gap = has_stop & ((cds["stop_start"] > cds["End"]) | (cds["stop_end"] < cds["Start"]))
    if gap.any():
        raise ValueError(
            "stop_codon rows do not touch the CDS of their exon in transcripts: "
            + ", ".join(sorted(cds.loc[gap, transcript_col].astype(str).unique()))
        )
    start_dtype, end_dtype = cds["Start"].dtype, cds["End"].dtype
    cds["Start"] = cds[["Start", "stop_start"]].min(axis=1).astype(start_dtype)
    cds["End"] = cds[["End", "stop_end"]].max(axis=1).astype(end_dtype)
    cds = cds.drop(columns=["stop_start", "stop_end"])

    # exons with stop codon bases but no CDS row, e.g. the second part of a split stop codon
    cds_keys = pd.MultiIndex.from_frame(cds[keys])
    stop_only = stops[~stops.index.isin(cds_keys)]
    extra = stop_rows.drop_duplicates(keys).set_index(keys).loc[stop_only.index].reset_index()
    extra["Start"] = stop_only["stop_start"].to_numpy().astype(start_dtype)
    extra["End"] = stop_only["stop_end"].to_numpy().astype(end_dtype)
    extra["Feature"] = "CDS"

    logger.info(
        "Stop codons from stop_codon rows: %d transcripts, %d of them with stop codon bases in an exon without CDS.",
        stops.index.get_level_values(transcript_col).nunique(),
        extra[transcript_col].nunique(),
    )
    coding = pd.concat([cds, extra[cds.columns]], ignore_index=True)
    coding["has_stop_codon"] = coding[transcript_col].isin(stop_rows[transcript_col])
    return coding


def read_fasta(fasta_path):
    """
    Reads a genome FASTA file using pyfaidx.Fasta and returns a pyfaidx.Fasta object.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    return Fasta(fasta_path)


def compute_exon_numbers(gtf):
    """
    Compute exon numbers for the Features exon, CDS and stop_codon in a GTF PyRanges object.
    Exon numbers are assigned based on genomic order per transcript and strand.
    CDS and stop_codon features inherit the exon number of the exon they overlap.

    On + Strand: Smallest exon number is the Start, Largest exon number is the End.
    On - Strand: Smallest exon number is the Start, Largest exon number is the End.
    (was different for hg19: the smallest exon number was the end, that is why we need to adjust it here.)

    :param gtf: PyRanges object of the GTF
    :return: PyRanges object with new column 'exon_number_computed'
    """
    gtf_df = gtf.df.copy()

    # A GTF read from file has exon_number as str (pandas 3) with missing values on features
    # without one. The computed numbers are ints, so hold the column as nullable integer.
    if "exon_number" in gtf_df.columns:
        gtf_df["exon_number"] = gtf_df["exon_number"].astype("Int64")
    else:
        gtf_df["exon_number"] = pd.Series(pd.NA, index=gtf_df.index, dtype="Int64")

    # Step 1: Compute exon numbers for exon features
    exons = gtf_df[gtf_df.Feature == "exon"].copy()
    for tx, group in exons.groupby("transcript_id"):
        strand = group["Strand"].iloc[0]
        if strand == "+":
            sorted_group = group.sort_values("Start")
        else:
            sorted_group = group.sort_values("Start", ascending=False)
        exons.loc[sorted_group.index, "exon_number"] = range(1, len(sorted_group) + 1)

    # Step 2: Assign exon numbers to CDS and stop_codon features
    cds = gtf_df[gtf_df.Feature.isin(["CDS", "stop_codon"])].copy()
    for tx, exon_group in exons.groupby("transcript_id"):
        cds_group = cds[cds.transcript_id == tx]
        for idx, cds_row in cds_group.iterrows():
            overlaps = exon_group[(exon_group["Start"] <= cds_row["End"]) & (exon_group["End"] >= cds_row["Start"])]
            if not overlaps.empty:
                # choose exon with maximum overlap
                overlap_idx = overlaps.apply(
                    lambda row: min(row["End"], cds_row["End"]) - max(row["Start"], cds_row["Start"]), axis=1
                ).idxmax()
                gtf_df.loc[idx, "exon_number"] = exon_group.loc[overlap_idx, "exon_number"]

    # Step 3: Update exon features
    gtf_df.loc[exons.index, "exon_number"] = exons["exon_number"]

    return pr.PyRanges(gtf_df)
