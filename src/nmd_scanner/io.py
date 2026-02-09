# Import dependencies
import logging
import os
import pandas as pd
import pyranges as pr
from pyfaidx import Fasta
from functools import lru_cache

logger = logging.getLogger(__name__)

__all__ = [
    "read_vcf", "read_gtf", "read_fasta"
    ]

# Global cache to map id --> actual fasta object
_fasta_cache = {}

def read_vcf(vcf_path):

    # TODO: adjust this function to also get structural variants and then also adjust downstream analysis
    #  (especially the inclusion of the Variants into the reference CDS sequence to create the alternative CDS)
    """
    Read a single VCF file into a PyRanges object with adjusted coordinates.
    """
    df = pd.read_csv(
        vcf_path,
        comment='#',
        sep='\t',
        header=None,
        names=['Chromosome', 'Start', 'ID', 'Ref',
               'Alt', 'Qual', 'Filter', 'Info']
    )

    # Adjust coordinates to 0-based
    df['Start'] = df['Start'] - 1
    df['End'] = df['Start'] + df['Ref'].str.len()

    # Keep only relevant columns
    gr = pr.PyRanges(df[['Chromosome', 'Start', 'End', 'ID',
                     'Ref', 'Alt', 'Qual', 'Filter', 'Info']])
    return gr


def read_gtf(gtf_path, reassign_exons=False, canonical_only=False):
    """
    Reads a GTF file into a PyRanges object.
    """
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")
    gtf = pr.read_gtf(gtf_path)

    if reassign_exons:
        # Adjust exon number in GTF (need this for the (old) hg19 version)
        logger.info("Adjusting exon numbers")
        gtf = compute_exon_numbers(gtf)
        logger.info("Exon numbers adjusted")
    
    if canonical_only:
        logger.info("Filtering to Ensembl canonical transcripts only")
        # Filter to canonical transcripts only based on Ensembl_canonical tag
        # The tag is expected to be in the format: "Ensembl_canonical,other
        # check if Ensembl_canonical is in the set of tags
        gtf = gtf[gtf.tag.apply(lambda x: False if pd.isna(x) else ('Ensembl_canonical' in x.split(',')))]
        logger.debug(f"After filtering to canonical transcripts, GTF shape: {gtf.df.shape}")

    return gtf


def read_fasta(fasta_path):
    """
    Reads a genome FASTA file using pyfaidx.Fasta and returns a pyfaidx.Fasta object.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    return Fasta(fasta_path)


@lru_cache(maxsize=None)
def _fetch_seq_cached(chrom, start, end, fasta_id):
    fasta = _fasta_cache[fasta_id]  # lookup the actual fasta object
    return fasta[chrom][start:end].seq.upper()

def extract_fasta_sequence(chrom, start, end, fasta):
    """
    For faster access of sequences from a FASTA reference using caching
    Returns original DataFrame with an additional column containing the extracted sequence per row
    """
    _fasta_cache[id(fasta)] = fasta  # store for later access
    return _fetch_seq_cached(chrom, start, end, id(fasta))


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
        exons.loc[sorted_group.index, "exon_number"] = range(
            1, len(sorted_group) + 1)

    # Step 2: Assign exon numbers to CDS features
    cds = gtf_df[gtf_df.Feature == "CDS"].copy()
    for tx, exon_group in exons.groupby("transcript_id"):
        cds_group = cds[cds.transcript_id == tx]
        for idx, cds_row in cds_group.iterrows():
            overlaps = exon_group[(exon_group["Start"] <= cds_row["End"]) & (
                exon_group["End"] >= cds_row["Start"])]
            if not overlaps.empty:
                # choose exon with maximum overlap
                overlap_idx = overlaps.apply(lambda row: min(
                    row["End"], cds_row["End"]) - max(row["Start"], cds_row["Start"]), axis=1).idxmax()
                gtf_df.loc[idx, "exon_number"] = exon_group.loc[overlap_idx, "exon_number"]

    # Step 3: Update exon features
    gtf_df.loc[exons.index, "exon_number"] = exons["exon_number"]

    return pr.PyRanges(gtf_df)