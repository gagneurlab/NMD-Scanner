# Import dependencies
import argparse
import functools
import json
import logging
import os
import warnings

import pandas as pd
import pyarrow.parquet as pq
import tqdm
from pyfaidx import Fasta

from nmd_scanner.extra_features import add_features_and_rules
from nmd_scanner.rules import extract_ptc
from nmd_scanner.scan import detect_annotation_format, read_annotation, read_vcf
from nmd_scanner.schema import OUTPUT_COLUMN_KINDS, SEQUENCE_COLUMNS, _arrow_schema, to_arrow

# Deprecated alias of schema.STOP_CODON_COLUMNS, kept for one release together with parquet_schema and
# to_parquet_safe. It goes in a later release.
from nmd_scanner.schema import STOP_CODON_COLUMNS as STOP_CODON_COLUMNS

SUPPORTED_OUTPUT_EXTENSIONS = (".csv", ".parquet", ".pq")

logger = logging.getLogger(__name__)


def annotate(
    vcf_path: str | os.PathLike,
    annotation_path: str | os.PathLike,
    fasta_path: str | os.PathLike,
    reassign_exons: bool = False,
    sequences: bool = True,
) -> pd.DataFrame:
    """
    Annotate the variants of a VCF file with NMD features and return the result table.

    Nothing is written to disk and logging is not configured. Use `write_results` to save the table, or `to_arrow`
    to convert it to a typed pyarrow Table.

    Steps:
    1. Read input files (VCF, FASTA, annotation)
    2. Assign exon numbers (optional, recommended for hg19)
    3. Parse and preprocess gene annotations (coding regions, i.e. CDS plus stop codon, and exons)
    4. Extract premature termination codons (PTCs) & Evaluate NMD escape rules
    5. Add extra features to output (e.g. 3' & 5'UTR length, downstream & upstream exon counts, etc.)

    :param vcf_path: path to the input VCF file
    :param annotation_path: path to the input gene annotation file (GFF3, optionally gzip-compressed)
    :param fasta_path: path to the reference FASTA file. It also shows whether a CDS ends in a stop codon, and for
                       an Ensembl GFF3 whether it starts with one.
    :param reassign_exons: recompute the exon numbers of the annotation (recommended for hg19; may be slow)
    :param sequences: keep the 4 sequence columns of SEQUENCE_COLUMNS. With False, they are left out, which saves
                      most of the table's memory. The other columns stay the same.
    :return: DataFrame summarizing all annotated variants, with the columns and dtypes of
             output_column_kinds(sequences) (see nmd_scanner.schema). It has zero rows if no variant gives a result.
    """

    # read VCF file (variants)
    logger.info("Reading VCF file: %s", vcf_path)
    vcf = read_vcf(vcf_path)
    logger.info("VCF shape: %s", vcf.shape)

    # read FASTA file (genome sequence)
    logger.info("Reading FASTA file: %s", fasta_path)
    fasta = Fasta(fasta_path)

    # read gene annotation file (GFF3) into exon rows and coding regions (CDS rows with has_start_codon and
    # has_stop_codon).
    # reassign_exons recomputes the exon numbers (recommended for hg19).
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

    # add_features_and_rules applies the schema with the sequences, so they are dropped after it
    if not sequences:
        results = results.drop(columns=list(SEQUENCE_COLUMNS))

    return results


def main(vcf_path, annotation_path, fasta_path, output, reassign_exons=False, sequences=True):
    """
    Main function for NMD scanner: annotate the variants and write the results to a file

    :param vcf_path: path to the input VCF file
    :param annotation_path: path to the input gene annotation file (GFF3, optionally gzip-compressed)
    :param fasta_path: path to the reference FASTA file
    :param output: path to the output file (.csv, .parquet, or .pq)
    :param reassign_exons: recompute the exon numbers of the annotation
    :param sequences: keep the 4 sequence columns of SEQUENCE_COLUMNS (see annotate())
    :return: DataFrame summarizing all annotated variants: the table that annotate() returns
    """

    results = annotate(vcf_path, annotation_path, fasta_path, reassign_exons=reassign_exons, sequences=sequences)

    # Write output
    logger.info("Writing results to %s", output)
    write_results(results, output)

    return results


def write_results(results, output):
    """
    Write the results DataFrame to a CSV or Parquet file based on the output extension.
    A Parquet file gets the table of ``to_arrow``, with the same types for every input. A CSV file gets each list
    column as JSON, e.g. [{"exon_number": 1, "length": 36}], which ``json.loads`` reads back.
    """

    ext = os.path.splitext(output)[1].lower()
    if ext == ".csv":
        _json_lists(results).to_csv(output, index=False)
    elif ext in (".parquet", ".pq"):
        pq.write_table(to_arrow(results), output)
    else:
        raise ValueError(f"Unsupported output extension: {ext!r}. Supported: {', '.join(SUPPORTED_OUTPUT_EXTENSIONS)}")


def _json_lists(results):
    """
    Return a copy of ``results`` in which each list column (kind pair_list, int_list or stop_codon_list) holds JSON
    text. A missing value stays missing, so CSV writes it as an empty field.
    """

    results = results.copy()
    for column in results.columns:
        if OUTPUT_COLUMN_KINDS.get(column) in ("pair_list", "int_list", "stop_codon_list"):
            results[column] = results[column].map(_json_list, na_action="ignore")
    return results


def _json_list(value):
    """Return a list value as JSON text. A numpy int, e.g. of a column read from Parquet, becomes a plain int."""

    return json.dumps(list(value), default=lambda item: item.item())


def parquet_schema(results):
    """
    Deprecated: use ``nmd_scanner.to_arrow``, whose Table has this schema. This alias goes in a later release.

    Return the pyarrow schema for the columns of ``results``, with the types of schema.KIND_ARROW_TYPES.
    A column that OUTPUT_COLUMN_KINDS does not list raises a KeyError.
    """

    warnings.warn(
        "nmd_scanner.cli.parquet_schema is deprecated, use nmd_scanner.to_arrow", DeprecationWarning, stacklevel=2
    )
    return _arrow_schema(results.columns)


def to_parquet_safe(results):
    """
    Deprecated: use ``nmd_scanner.to_arrow``, which does this conversion. This alias goes in a later release.

    Return ``results`` itself. Its stop codon columns (schema.STOP_CODON_COLUMNS) hold {"position": ..., "codon": ...}
    records already, which this function made from (position, codon) tuples before 0.4.0.
    """

    warnings.warn(
        "nmd_scanner.cli.to_parquet_safe is deprecated, use nmd_scanner.to_arrow", DeprecationWarning, stacklevel=2
    )
    return results


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
            "GENCODE and Ensembl GFF3 flavors are supported."
        ),
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help=(
            "Path to reference genome FASTA file. It also shows whether a CDS ends in a stop codon, "
            "and for an Ensembl GFF3 whether it starts with one."
        ),
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
    parser.add_argument(
        "--no-sequences",
        dest="sequences",
        action="store_false",
        help=(
            "Leave out the 4 sequence columns ref_cds_seq, alt_cds_seq, transcript_seq and alt_transcript_seq. "
            "They make up most of the output size."
        ),
    )

    args = parser.parse_args()

    if args.annotation == "":
        parser.error("argument --annotation: expected a path, got an empty string")
    # an unknown suffix fails here, before the VCF and FASTA are read
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
    main(
        args.vcf,
        args.annotation,
        args.fasta,
        args.output,
        reassign_exons=args.reassign_exons,
        sequences=args.sequences,
    )


if __name__ == "__main__":
    main_cli()
