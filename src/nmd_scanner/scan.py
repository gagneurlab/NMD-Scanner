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

    # Adjust coordinates to 0-based, half-open (BED-like)
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

