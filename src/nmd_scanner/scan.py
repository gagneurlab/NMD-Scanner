# Import dependencies
import os
import pandas as pd
import pyranges as pr
from pyfaidx import Fasta

# Create the functions used for reading in the files (VCF, GTF, FASTA)

def list_vcf_files(vcf_path):
    """
    Accepts either a path to a VCF file or a directory containing VCF files.

    Parameters:
        vcf_path (str): Path to a VCF file or directory.

    Returns:
        Tuple[List[str], bool]: (VCF file paths, is_single_file)

    Raises:
        FileNotFoundError: If no valid VCF file(s) are found.
        ValueError: If the input path is invalid.
    """
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF path does not exist: {vcf_path}")

    if os.path.isfile(vcf_path):
        if vcf_path.endswith(".vcf") or vcf_path.endswith(".vcf.gz"):
            return [vcf_path], True
        else:
            raise ValueError(f"Provided file is not a VCF file: {vcf_path}")

    elif os.path.isdir(vcf_path):
        all_files = os.listdir(vcf_path)
        vcf_files = [os.path.join(vcf_path, f) for f in all_files if f.endswith(".vcf") or f.endswith(".vcf.gz")]
        if not vcf_files:
            raise FileNotFoundError(f"No VCF files found in the directory: {vcf_path}")
        return vcf_files, False

    else:
        raise ValueError(f"Provided path is neither a file nor a directory: {vcf_path}")


def read_vcf(vcf_path):
    """
    Read a single VCF file into a PyRanges object with adjusted coordinates.

    Parameters:
        vcf_path (str): Path to the VCF file.

    Returns:
        PyRanges: Genomic intervals from VCF.
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
    """Reads a GTF file into a PyRanges object."""
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

