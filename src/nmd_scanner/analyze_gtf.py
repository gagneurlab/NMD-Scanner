# analyze_gtf.py is not needed anymore

import pandas as pd
import pyranges as pr
from nmd_scanner.scan import read_gtf

def extract_features(gr: pr.PyRanges, feature: str, chromosomes=None) -> pd.DataFrame:
    """
    Filter a PyRanges by feature type (e.g., 'exon', 'CDS', 'transcript'),
    and optional list of chromosomes, then return the underlying DataFrame.
    """
    df = gr.df
    df = df[df.Feature == feature]
    if chromosomes is not None:
        df = df[df.Chromosome.isin(chromosomes)]
    return df.copy()

def annotate_exon_lengths(exons_df: pd.DataFrame) -> pd.DataFrame:
    """Add an 'exon_length' column = End - Start."""
    exons_df = exons_df.copy()
    exons_df['exon_length'] = exons_df['End'] - exons_df['Start']
    return exons_df

def count_exons_per_transcript(exons_df: pd.DataFrame) -> pd.DataFrame:
    """
    Given an exon DataFrame (with a 'transcript_id' column),
    returns a DF with columns ['transcript_id', 'num_exons_per_transcript'].
    """
    counts = (
        exons_df
        .groupby('transcript_id')
        .size()
        .reset_index(name='num_exons_per_transcript')
    )
    return counts

def merge_exon_counts(exons_df: pd.DataFrame, exon_counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge exon counts back onto the exons_df, filling missing with 0.
    """
    merged = exons_df.merge(
        exon_counts_df,
        on='transcript_id',
        how='left'
    )
    merged['num_exons_per_transcript'] = (
        merged['num_exons_per_transcript']
        .fillna(0)
        .astype(int)
    )
    return merged

def run_gtf_analysis(gtf_path, chromosomes=None):
    """Full workflow from GTF → exons with lengths and exon counts per transcript."""
    gtf = read_gtf(gtf_path)
    exons_df = extract_features(gtf, 'exon', chromosomes)
    exons_with_len = annotate_exon_lengths(exons_df)
    exon_counts = count_exons_per_transcript(exons_with_len)
    merged = merge_exon_counts(exons_with_len, exon_counts)
    return {
        'exons': exons_df,
        'exon_lengths': exons_with_len,
        'exon_counts': exon_counts,
        'merged': merged
    }