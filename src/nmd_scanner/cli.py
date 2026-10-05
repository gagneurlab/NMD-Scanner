# Import dependencies
import argparse
import functools
import logging
import os

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import tqdm
from pyfaidx import Fasta

from nmd_scanner.extra_features import add_features_and_rules
from nmd_scanner.rules import extract_ptc
from nmd_scanner.scan import detect_annotation_format, read_annotation, read_vcf
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


def annotate(vcf_path, annotation_path, fasta_path, reassign_exons=False):
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
    :param annotation_path: path to the input gene annotation file (GFF3, optionally gzip-compressed)
    :param fasta_path: path to the reference FASTA file. It also shows whether a CDS ends in a stop codon.
    :param reassign_exons: recompute the exon numbers of the annotation (recommended for hg19; may be slow)
    :return: DataFrame summarizing all annotated variants, with the columns and dtypes of OUTPUT_COLUMN_KINDS
             (see nmd_scanner.schema). It has zero rows if no variant gives a result.
    """

    # read VCF file (variants)
    logger.info("Reading VCF file: %s", vcf_path)
    vcf = read_vcf(vcf_path)
    logger.info("VCF shape: %s", vcf.shape)

    # read FASTA file (genome sequence)
    logger.info("Reading FASTA file: %s", fasta_path)
    fasta = Fasta(fasta_path)

    # read gene annotation file (GFF3) into exon rows and coding regions (CDS rows with has_stop_codon).
    # reassign_exons recomputes the exon numbers (need this for the (old) hg19 version).
    logger.info("Reading annotation file: %s", annotation_path)
    annotation = read_annotation(annotation_path, fasta, reassign_exons=reassign_exons)
    logger.info("Annotation file shape: %s", annotation.shape)
    cds_df = annotation[annotation["Feature"] == "CDS"]

    # extract exon regions from the annotation and compute exon related metrics:
    # exon length & number of exons contained in each transcript
    exons_df = annotation[annotation["Feature"] == "exon"].copy()
    exons_df["exon_length"] = exons_df["End"] - exons_df["Start"]

    # Create reference and alternative CDS and transcript sequences (+ metadata) and analyze for start and stop codons & -loss
    logger.info("Creating sequences and analyzing...")
    results = extract_ptc(cds_df, vcf, fasta, exons_df)

    # Add the NMD features (inspired by the NMD efficiency benchmark dataset) and the NMD escape rules
    results = add_features_and_rules(results)

    return results


def main(vcf_path, annotation_path, fasta_path, output, reassign_exons=False):
    """
    Main function for NMD scanner: annotate the variants and write the results to a file

    :param vcf_path: path to the input VCF file
    :param annotation_path: path to the input gene annotation file (GFF3, optionally gzip-compressed)
    :param fasta_path: path to the reference FASTA file
    :param output: path to the output file (.csv, .parquet, or .pq)
    :param reassign_exons: recompute the exon numbers of the annotation
    :return: DataFrame summarizing all annotated variants: the table that annotate() returns
    """

    results = annotate(vcf_path, annotation_path, fasta_path, reassign_exons=reassign_exons)

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
        table = pa.Table.from_pandas(to_parquet_safe(results), schema=parquet_schema(results), preserve_index=False)
        pq.write_table(table, output)
    else:
        raise ValueError(f"Unsupported output extension: {ext!r}. Supported: {', '.join(SUPPORTED_OUTPUT_EXTENSIONS)}")


def parquet_schema(results):
    """
    Return the pyarrow schema for the columns of ``results``, with the types listed in
    ``OUTPUT_COLUMN_KINDS``. A column that is not listed raises a KeyError.
    """

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


def _set_up_process():
    """
    Sets up the process that the CLI owns: NMD-Scanner logs its INFO messages to stderr, and all
    other loggers their warnings. tqdm shows no progress bars.
    """

    logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(levelname)s %(name)s: %(message)s", force=True)
    # only NMD-Scanner logs at INFO: at INFO, the Rust code of polars-bio logs too
    logging.getLogger("nmd_scanner").setLevel(logging.INFO)
    # under `python -m nmd_scanner.cli`, this module logs as __main__
    logger.setLevel(logging.INFO)
    # polars-bio shows a tqdm bar for every read, also in a pipe or a log file
    tqdm.tqdm.__init__ = functools.partialmethod(tqdm.tqdm.__init__, disable=True)


def main_cli():
    """Console-script entry point: parse arguments and run the pipeline."""

    _set_up_process()

    parser = argparse.ArgumentParser(description="Run NMD pipeline")
    parser.add_argument("--vcf", required=True, help="Path to VCF file")
    parser.add_argument(
        "--annotation",
        required=True,
        help=(
            "Path to gene annotation file (GFF3 with a .gff3 or .gff suffix, optionally gzip-compressed). "
            "GENCODE and Ensembl GFF3 flavors are supported. GTF is not supported: use the GFF3 of the same release."
        ),
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to reference genome FASTA file. It also shows whether a CDS ends in a stop codon.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help=(
            "Path to the output file. Extension determines format: "
            ".csv for CSV, .parquet or .pq for Parquet. "
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
    # a GTF or an unknown suffix fails here, before the VCF and FASTA are read
    try:
        detect_annotation_format(args.annotation)
    except ValueError as error:
        parser.error(f"argument --annotation: {error}")

    # Check that the output path is valid
    if not is_valid_output_path(args.output):
        raise SystemExit(
            f"Invalid output path: {args.output!r}. "
            f"Must be a non-existing or overwritable file with one of these extensions: "
            f"{', '.join(SUPPORTED_OUTPUT_EXTENSIONS)}, "
            f"and its parent directory must exist."
        )

    # Run the main pipeline
    main(args.vcf, args.annotation, args.fasta, args.output, reassign_exons=args.reassign_exons)


if __name__ == "__main__":
    main_cli()
