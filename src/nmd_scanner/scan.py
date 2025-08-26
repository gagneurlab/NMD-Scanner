# Import dependencies
import os
import pandas as pd
import pyranges as pr
from pyfaidx import Fasta

# Create the functions used for reading in the files (VCF, GTF, FASTA)

def read_vcf(vcf_path):
    """
    Read a single VCF file into a PyRanges object with adjusted coordinates.
    """
    df = pd.read_csv(
        vcf_path,
        comment='#',
        sep='\t',
        header=None,
        names=['Chromosome', 'Start', 'ID', 'Ref', 'Alt', 'Qual', 'Filter', 'Info']
    )

    # Adjust coordinates to 0-based
    df['Start'] = df['Start'] - 1
    df['End'] = df['Start'] + df['Ref'].str.len()

    # Keep only relevant columns
    gr = pr.PyRanges(df[['Chromosome', 'Start', 'End', 'ID', 'Ref', 'Alt', 'Qual', 'Filter', 'Info']])
    return gr

def read_gtf(gtf_path):
    """
    Reads a GTF file into a PyRanges object.
    """
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")
    return pr.read_gtf(gtf_path)

def read_fasta(fasta_path):
    """
    Reads a genome FASTA file using pyfaidx.Fasta and returns a pyfaidx.Fasta object.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    return Fasta(fasta_path)

def compute_exon_numbers(gtf):
    """
    Compute exon numbers for the Features exon and CDS in a GTF PyRanges object.
    Exon numbers are assigned based on genomic order per transcript and strand.
    CDS features inherit the exon number of the exon they overlap.

    On + Strand: Smallest exon number is the Start, Largest exon number is the End.
    On - Strand: Smallest exon number is the Start, Largest exon number is the End.
    (was different for hg19: the smallest exon number was the end, that is why we need to adjust it here.)

    :param gtf: PyRanges object of the GTF
    :return: PyRanges object with new column 'exon_number_computed'
    """
    gtf_df = gtf.df.copy()

    # Step 1: Compute exon numbers for exon features
    exons = gtf_df[gtf_df.Feature == "exon"].copy()
    for tx, group in exons.groupby("transcript_id"):
        strand = group["Strand"].iloc[0]
        if strand == "+":
            sorted_group = group.sort_values("Start")
        else:
            sorted_group = group.sort_values("Start", ascending=False)
        exons.loc[sorted_group.index, "exon_number"] = range(1, len(sorted_group) + 1)

    # Step 2: Assign exon numbers to CDS features
    cds = gtf_df[gtf_df.Feature == "CDS"].copy()
    for tx, exon_group in exons.groupby("transcript_id"):
        cds_group = cds[cds.transcript_id == tx]
        for idx, cds_row in cds_group.iterrows():
            overlaps = exon_group[(exon_group["Start"] <= cds_row["End"]) & (exon_group["End"] >= cds_row["Start"])]
            if not overlaps.empty:
                # choose exon with maximum overlap
                overlap_idx = overlaps.apply(lambda row: min(row["End"], cds_row["End"]) - max(row["Start"], cds_row["Start"]), axis=1).idxmax()
                gtf_df.loc[idx, "exon_number"] = exon_group.loc[overlap_idx, "exon_number"]

    # Step 3: Update exon features
    gtf_df.loc[exons.index, "exon_number"] = exons["exon_number"]

    return pr.PyRanges(gtf_df)