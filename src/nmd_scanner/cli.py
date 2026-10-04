# Import dependencies
import argparse
import logging
import os

import pandas as pd
from pyfaidx import Fasta

from nmd_scanner.extra_features import add_features_and_rules
from nmd_scanner.rules import extract_ptc
from nmd_scanner.scan import read_annotation, read_vcf
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


def annotate(vcf_path, annotation_path, fasta_path, reassign_exons=False, annotation_format=None):
    """
    Annotate the variants of a VCF file with NMD features and return the result table.

    Nothing is written to disk and logging is not configured. Use `write_results` to save the table.

    Steps:
    1. Read input files (VCF, FASTA, annotation)
    2. Assign exon numbers (optional, recommended for hg19)
    3. Parse and preprocess gene annotations (coding regions, i.e. CDS plus stop codon, and exons)
    4. Extract premature termination codons (PTCs) & Evaluate NMD escape rules
    5. Add extra features to output (e.g. 3' & 5'UTR length, downstream & upstream exon counts, etc.)

    :param vcf_path: path to the input VCF file
    :param annotation_path: path to the input gene annotation file (GTF or GFF3, optionally
                            gzip-compressed)
    :param fasta_path: path to the reference FASTA file. For a GFF3 annotation, it also shows whether a CDS ends
                       in a stop codon.
    :param reassign_exons: recompute the exon numbers of the annotation (recommended for hg19; may be slow)
    :param annotation_format: "gtf" or "gff3", or None to detect the format from the file suffix
    :return: DataFrame summarizing all annotated variants, with the columns and dtypes of OUTPUT_COLUMN_KINDS
             (see nmd_scanner.schema). It has zero rows if no variant gives a result.
    """

    # read VCF file (variants)
    logger.info("Reading VCF file: %s", vcf_path)
    vcf = read_vcf(vcf_path)
    logger.info("VCF shape: %s", vcf.df.shape)

    # read FASTA file (genome sequence)
    logger.info("Reading FASTA file: %s", fasta_path)
    fasta = Fasta(fasta_path)

    # read gene annotation file (GTF or GFF3) into exon rows and coding regions (CDS rows with has_stop_codon).
    # reassign_exons recomputes the exon numbers (need this for the (old) hg19 version).
    logger.info("Reading annotation file: %s", annotation_path)
    gtf = read_annotation(annotation_path, fasta, fmt=annotation_format, reassign_exons=reassign_exons)
    logger.info("Annotation file shape: %s", gtf.df.shape)
    cds_df = gtf[gtf.Feature == "CDS"].df

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

    return results


def main(vcf_path, gtf_path, fasta_path, output, reassign_exons=False, annotation_path=None):
    """
    Main function for NMD scanner: annotate the variants and write the results to a file

    :param vcf_path: path to the input VCF file
    :param gtf_path: path to the input GTF file (whatever its file name; optionally gzip-compressed).
                      Give either this or ``annotation_path``, the other one as None.
    :param fasta_path: path to the reference FASTA file
    :param output: path to the output file (.csv, .parquet, or .pq)
    :param reassign_exons: recompute the exon numbers of the annotation
    :param annotation_path: path to the input gene annotation file (GTF or GFF3, optionally
                      gzip-compressed; format is auto-detected from the file suffix), instead of gtf_path
    :return: DataFrame summarizing all annotated variants: the table that annotate() returns
    """

    if (gtf_path is None) == (annotation_path is None):
        raise ValueError("Give exactly one of gtf_path and annotation_path.")

    if gtf_path is None:
        results = annotate(vcf_path, annotation_path, fasta_path, reassign_exons=reassign_exons)
    else:
        results = annotate(vcf_path, gtf_path, fasta_path, reassign_exons=reassign_exons, annotation_format="gtf")

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
    annotation_group = parser.add_mutually_exclusive_group(required=True)
    annotation_group.add_argument(
        "--annotation",
        help=(
            "Path to gene annotation file (GTF or GFF3, optionally gzip-compressed). "
            "Format is auto-detected from the file suffix. GENCODE and Ensembl GFF3 flavors "
            "are supported."
        ),
    )
    annotation_group.add_argument(
        "--gtf",
        help="Path to a GTF file, read as GTF whatever its file name. Deprecated: use --annotation.",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to reference genome FASTA file. For a GFF3 annotation, it also shows whether a CDS ends in a stop codon.",
    )
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

    if args.annotation == "":
        parser.error("argument --annotation: expected a path, got an empty string")

    # Check that the output path is valid
    if not is_valid_output_path(args.output):
        raise SystemExit(
            f"Invalid output path: {args.output!r}. "
            f"Must be a non-existing or overwritable file with one of these extensions: "
            f"{', '.join(SUPPORTED_OUTPUT_EXTENSIONS)}, "
            f"and its parent directory must exist."
        )

    # Run the main pipeline
    main(
        args.vcf,
        args.gtf,
        args.fasta,
        args.output,
        reassign_exons=args.reassign_exons,
        annotation_path=args.annotation,
    )


if __name__ == "__main__":
    main_cli()
