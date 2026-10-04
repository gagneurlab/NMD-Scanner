# Import dependencies
import argparse
import logging
import os

import pandas as pd
from pyfaidx import Fasta

from nmd_scanner.extra_features import add_features_and_rules
from nmd_scanner.rules import extract_ptc
from nmd_scanner.scan import compute_exon_numbers, read_gtf, read_vcf
from nmd_scanner.schema import OUTPUT_COLUMN_KINDS

SUPPORTED_OUTPUT_EXTENSIONS = (".csv", ".parquet", ".pq")

# Columns that hold lists of (position, codon) tuples, e.g. (5442, "TGA"). pyarrow's
# pandas conversion treats each tuple as a flat, homogeneously-typed sub-list rather
# than a struct: it infers the element type from the tuple's first field (an int) and
# then fails on the second, string field. Parquet output needs these turned into
# {"position": ..., "codon": ...} records instead, so pyarrow can infer
# list<struct<position: int64, codon: string>>. CSV output and the in-memory results
# table are unaffected; only the parquet copy is rewritten.
STOP_CODON_COLUMNS = ("ref_all_stop_codons", "alt_all_stop_codons", "transcript_all_stop_codons")

logger = logging.getLogger(__name__)


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
    :return: DataFrame summarizing all annotated variants, with the columns and dtypes of OUTPUT_COLUMN_KINDS
             (see nmd_scanner.schema). It has zero rows if no variant gives a result.
    """

    # read VCF file (variants)
    logger.info("Reading VCF file: %s", vcf_path)
    vcf = read_vcf(vcf_path)
    logger.info("VCF shape: %s", vcf.df.shape)

    # read GTF file (gene annotation)
    logger.info("Reading GTF file: %s", gtf_path)
    gtf = read_gtf(gtf_path)
    logger.info("GTF File shape: %s", gtf.df.shape)

    # read FASTA file (genome sequence)
    logger.info("Reading FASTA file: %s", fasta_path)
    fasta = Fasta(fasta_path)

    # Adjust exon number in GTF (need this for the (old) hg19 version)
    if reassign_exons:
        logger.info("Adjust exon numbers")
        gtf = compute_exon_numbers(gtf)
        logger.info("Exon numbers adjusted.")

    # extract CDS regions from the GTF file
    cds = gtf[gtf.Feature == "CDS"]
    cds_df = cds.df

    # extract exon regions from the GTF file and compute exon related metrics:
    # exon length & number of exons contained in each transcript
    exons = gtf[gtf.Feature == "exon"]
    exons_df = exons.df
    exons_df["exon_length"] = exons_df["End"] - exons_df["Start"]

    # Create reference and alternative CDS and transcript sequences (+ metadata) and analyze for start and stop codons & -loss
    logger.info("Creating sequences and analyzing...")
    results = extract_ptc(cds_df, vcf, fasta, exons_df)

    # Add the NMD features (inspired by the NMD efficiency benchmark dataset) and the NMD escape rules
    results = add_features_and_rules(results)

    # Write output
    logger.info("Writing results to %s", output)
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
            import pyarrow as pa
            import pyarrow.parquet as pq

            table = pa.Table.from_pandas(to_parquet_safe(results), schema=parquet_schema(results), preserve_index=False)
            pq.write_table(table, output)
        except ImportError as e:
            raise SystemExit(
                f'Writing parquet requires pyarrow. Install it via: pip install "nmd_scanner[parquet]"\nOriginal error: {e}'
            ) from e
    else:
        raise ValueError(f"Unsupported output extension: {ext!r}. Supported: {', '.join(SUPPORTED_OUTPUT_EXTENSIONS)}")


def parquet_schema(results):
    """
    Return the pyarrow schema for the columns of ``results``, with the types listed in
    ``OUTPUT_COLUMN_KINDS``. A column that is not listed raises a KeyError.
    """

    import pyarrow as pa

    stop_codon = pa.struct([pa.field("position", pa.int64()), pa.field("codon", pa.string())])
    kind_types = {
        "string": pa.string(),
        "int": pa.int64(),
        "bool": pa.bool_(),
        "pair_list": pa.list_(pa.list_(pa.int64())),
        "int_list": pa.list_(pa.int64()),
        "stop_codon_list": pa.list_(stop_codon),
    }
    return pa.schema([pa.field(column, kind_types[OUTPUT_COLUMN_KINDS[column]]) for column in results.columns])


def to_parquet_safe(results):
    """
    Return a copy of ``results`` with the stop-codon columns given a parquet-friendly,
    typed representation. See ``STOP_CODON_COLUMNS`` for why this is needed. Every other
    column, and the ``results`` table passed in, is left untouched.
    """

    columns_present = [column for column in STOP_CODON_COLUMNS if column in results.columns]
    if not columns_present:
        return results

    results = results.copy()
    for column in columns_present:
        results[column] = results[column].apply(_stop_codons_to_records)
    return results


def _stop_codons_to_records(stop_codons):
    """
    Turn a list of (position, codon) tuples into {"position": ..., "codon": ...} records.
    A missing value (None, np.nan, pd.NA) stays missing.
    """

    if pd.api.types.is_scalar(stop_codons) and pd.isna(stop_codons):
        return None
    return [{"position": position, "codon": codon} for position, codon in stop_codons]


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

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    parser = argparse.ArgumentParser(description="Run NMD pipeline")
    parser.add_argument("--vcf", required=True, help="Path to VCF file")
    parser.add_argument("--gtf", required=True, help="Path to GTF file")
    parser.add_argument("--fasta", required=True, help="Path to FASTA file")
    parser.add_argument(
        "--output",
        required=True,
        help=(
            "Path to the output file. Extension determines format: "
            ".csv for CSV, .parquet or .pq for Parquet (requires the parquet extra). "
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
