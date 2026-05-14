# Import dependencies
import argparse
import os

import pandas as pd
from pyfaidx import Fasta

from nmd_scanner.extra_features import add_nmd_features, evaluate_nmd_escape_rules
from nmd_scanner.rules import extract_ptc
from nmd_scanner.scan import compute_exon_numbers, read_gtf, read_vcf

SUPPORTED_OUTPUT_EXTENSIONS = (".csv", ".parquet", ".pq")


def main(vcf_path, gtf_path, fasta_path, output, reassign_exons=False):
    """
    Main function for NMD scanner

    Steps:
    1. Read input files (VCF, GTF, FASTA)
    2. Assign exon numbers (optional, recommended for hg19)
    3. Parse and preprocess gene annotations (CDS, exons)
    4. Extract premature termination codons (PTCs) & Evaluate NMD escape rules
    5. Add extra features to output (e.g. 3' & 5'UTR length, downstream & upstream exon counts, etc.)
    6. Return and Save output results

    :param vcf_path: path to the input VCF file
    :param gtf_path: path to the input GTF annotation file
    :param fasta_path: path to the reference FASTA file
    :param output: path to the output file (.csv, .parquet, or .pq)
    :return: DataFrame summarizing all annotated variants
    """

    # read VCF file (variants)
    print(f"Reading VCF file: {vcf_path}")
    vcf = read_vcf(vcf_path)
    print(f"VCF shape: {vcf.df.shape}")

    # read GTF file (gene annotation)
    print(f"Reading GTF file: {gtf_path}")
    gtf = read_gtf(gtf_path)
    print(f"GTF File shape: {gtf.df.shape}")

    # read FASTA file (genome sequence)
    print(f"Reading FASTA file: {fasta_path}")
    fasta = Fasta(fasta_path)

    # Adjust exon number in GTF (need this for the (old) hg19 version)
    if reassign_exons:
        print("Adjust exon numbers")
        gtf = compute_exon_numbers(gtf)
        print("Exon numbers adjusted.")

    # extract CDS regions from the GTF file
    cds = gtf[gtf.Feature == "CDS"]
    cds_df = cds.df

    # extract exon regions from the GTF file and compute exon related metrics:
    # exon length & number of exons contained in each transcript
    exons = gtf[gtf.Feature == "exon"]
    exons_df = exons.df
    exons_df["exon_length"] = exons_df["End"] - exons_df["Start"]

    # Create reference and alternative CDS and transcript sequences (+ metadata) and analyze for start and stop codons & -loss
    print("Creating sequences and analyzing...")
    results = extract_ptc(cds_df, vcf, fasta, exons_df)

    # Add additional features (inspired by NMD efficiency benchmark dataset)
    extra_features = results.apply(add_nmd_features, axis=1, result_type="expand")
    results = pd.concat([results, extra_features], axis=1)

    # Compute NMD-rules as last step
    nmd_results = results.apply(evaluate_nmd_escape_rules, axis=1, result_type="expand")
    results = pd.concat([results, nmd_results], axis=1)

    # Write output
    print(f"Writing results to {output}")
    write_results(results, output)

    return results


def write_results(results, output):
    """
    Write the results DataFrame to a CSV or Parquet file based on the output extension.
    """

    ext = os.path.splitext(output)[1].lower()
    if ext == ".csv":
        results.to_csv(output, index=False)
    elif ext in (".parquet", ".pq"):
        try:
            results.to_parquet(output, index=False)
        except ImportError as e:
            raise SystemExit(
                f"Writing parquet requires pyarrow. Install it via: pip install pyarrow\nOriginal error: {e}"
            ) from e
    else:
        raise ValueError(f"Unsupported output extension: {ext!r}. Supported: {', '.join(SUPPORTED_OUTPUT_EXTENSIONS)}")


def is_valid_output_path(path):
    """
    Validate that ``path`` is a writable output file path.

    Rules:
    - Must not point at an existing directory (single-file output only).
    - Parent directory must already exist (we do not create directories).
    - Extension must be one of the supported formats.
    """

    if os.path.isdir(path):
        return False
    parent = os.path.dirname(path) or "."
    if not os.path.isdir(parent):
        return False
    ext = os.path.splitext(path)[1].lower()
    return ext in SUPPORTED_OUTPUT_EXTENSIONS


def main_cli():
    """Console-script entry point: parse arguments and run the pipeline."""

    parser = argparse.ArgumentParser(description="Run NMD pipeline")
    parser.add_argument("--vcf", required=True, help="Path to VCF file")
    parser.add_argument("--gtf", required=True, help="Path to GTF file")
    parser.add_argument("--fasta", required=True, help="Path to FASTA file")
    parser.add_argument(
        "--output",
        required=True,
        help=(
            "Path to the output file. Extension determines format: "
            ".csv for CSV, .parquet or .pq for Parquet (requires pyarrow). "
            "Parent directory must exist; the file is overwritten if present."
        ),
    )

    # If user adds flag, reassign exon numbers
    parser.add_argument(
        "--reassign_exons", action="store_true", help="Recompute exon numbers (recommended for hg19; may be slow)"
    )

    args = parser.parse_args()

    # Check that the output path is valid
    if not is_valid_output_path(args.output):
        raise SystemExit(
            f"Invalid output path: {args.output!r}. "
            f"Must be a non-existing or overwritable file with one of these extensions: "
            f"{', '.join(SUPPORTED_OUTPUT_EXTENSIONS)}, "
            f"and its parent directory must exist."
        )

    # Run the main pipeline
    main(args.vcf, args.gtf, args.fasta, args.output, reassign_exons=args.reassign_exons)


if __name__ == "__main__":
    main_cli()
