"""
Schema of the result table: its columns, their order and their dtypes.

Each function that returns results gives them this schema, also for a table without rows. The section "Output
columns" of "Technical Notes.md" gives the meaning of each column and says when it is null. ``to_arrow`` converts
a result to a pyarrow Table with the same types for every input.
"""

import pandas as pd
import pyarrow as pa

from nmd_scanner.variant_placement import EXON_BOUNDARY_AMBIGUOUS, SPLICE_SITE_DESTROYED

# The values of nmd_model_status, in the order they are checked. A row gets the first value whose condition holds:
# - unknown_effect: unknown_reason is set. The alt transcript is unknown, so alt_has_ptc and 15 of the 19 model
#   inputs are null. unknown_reason says why.
# - no_ptc: alt_has_ptc is not True, so there is no PTC to score.
# - ref_ptc: ref_has_ptc is True. The reference has a PTC already, so the variant does not create it.
# - no_annotated_stop: has_stop_codon is False, which makes annotated_stop_distance and utr3_length null.
# - no_annotated_start: has_start_codon is False, which makes ptc_to_start_codon null. The true start codon lies
#   upstream of the CDS, at an unknown distance (e.g. cds_start_NF).
# - start_loss: start_loss is True and ptc_to_start_codon is null. After a start loss, the scan of alt_transcript_seq
#   takes the next ATG, and ptc_to_start_codon runs from it to the PTC. A scanned row is a PTC row only if the scan
#   finds an ATG and a stop codon after it, so it always has a ptc_to_start_codon. Only a start-loss row without
#   alt_transcript_seq gets this value, e.g. one of a transcript without exon rows. It is not scanned and keeps the
#   PTC of the alt CDS.
# - missing_input: another model input is null.
# - ok: the variant creates the PTC, and no model input is null.
# no_annotated_stop and no_annotated_start hold for every variant of the transcript. So they go before the values
# start_loss and missing_input, which depend on the variant.
MODEL_STATUSES = (
    "unknown_effect",
    "no_ptc",
    "ref_ptc",
    "no_annotated_stop",
    "no_annotated_start",
    "start_loss",
    "missing_input",
    "ok",
)

# The values of stop_classification: the path of rules.analyze_transcript that set alt_has_ptc, stop_loss and
# annotated_stop_distance. alt_transcript: the first in-frame stop codon of the alt transcript. start_loss_scan: the
# rescued ORF after a start loss. alt_cds: the codon scan of the alt CDS, for a row that keeps the flags from the CDS.
STOP_CLASSIFICATIONS = ("alt_transcript", "start_loss_scan", "alt_cds")

# The categories of each categorical kind: the closed value set of the column whose name the kind has. A categorical
# column holds one of these values or a null, and it has all of them as categories, also if a value does not occur.
CATEGORIES = {
    "strand": ("+", "-"),
    "unknown_reason": (SPLICE_SITE_DESTROYED, EXON_BOUNDARY_AMBIGUOUS),
    "nmd_model_status": MODEL_STATUSES,
    "stop_classification": STOP_CLASSIFICATIONS,
}

# pandas dtype of each column kind. The int, bool and string dtypes are nullable, so a missing value is
# pd.NA, and a column with missing values keeps its dtype. The string dtype has "python" storage: it
# behaves the same on pandas 2 and 3, with or without pyarrow. The default string dtype of pandas 3,
# "str", marks a missing value as NaN instead. The list kinds hold Python lists in object columns:
# "pair_list" holds {"exon_number": ..., "length": ...} records (dicts), "int_list" exon numbers and
# "stop_codon_list" {"position": ..., "codon": ...} records. A record has the field names of its Arrow
# struct (KIND_ARROW_TYPES). So a value keeps the shape of its elements through to_arrow, Parquet and
# pd.read_parquet with default arguments. Only the container changes: pd.read_parquet gives a numpy array of
# the same dicts, or of numpy ints. A categorical kind (CATEGORIES) has a CategoricalDtype with its fixed
# categories, so a missing value is NaN.
KIND_DTYPES = {
    "string": pd.StringDtype("python"),
    "int": pd.Int64Dtype(),
    "bool": pd.BooleanDtype(),
    "pair_list": object,
    "int_list": object,
    "stop_codon_list": object,
    **{kind: pd.CategoricalDtype(list(categories)) for kind, categories in CATEGORIES.items()},
}

# Arrow type of each column kind, which to_arrow and the Parquet output use. A record of a pair_list or a
# stop_codon_list becomes a struct with the same field names. A categorical kind becomes a dictionary of strings
# with int8 indices.
KIND_ARROW_TYPES = {
    "string": pa.string(),
    "int": pa.int64(),
    "bool": pa.bool_(),
    "pair_list": pa.list_(pa.struct([pa.field("exon_number", pa.int64()), pa.field("length", pa.int64())])),
    "int_list": pa.list_(pa.int64()),
    "stop_codon_list": pa.list_(pa.struct([pa.field("position", pa.int64()), pa.field("codon", pa.string())])),
    **dict.fromkeys(CATEGORIES, pa.dictionary(pa.int8(), pa.string())),
}

# Kind of every column that extract_ptc returns, in output order
PTC_COLUMN_KINDS = {
    "transcript_id": "string",
    "variant_id": "string",
    "cds_start": "int",
    "cds_end": "int",
    "ref_cds_seq": "string",
    "ref_cds_length": "int",
    "alt_cds_seq": "string",
    "alt_cds_length": "int",
    "chromosome": "string",
    "gene_id": "string",
    "strand": "strand",
    "has_start_codon": "bool",
    "has_stop_codon": "bool",
    "cds_frame": "int",
    "ref": "string",
    "alt": "string",
    "variant_start": "int",
    "variant_end": "int",
    "ref_cds_exons": "pair_list",
    "alt_cds_exons": "pair_list",
    "cds_in_transcript": "bool",
    "start_codon_exon": "int",
    "ref_last_codon": "string",
    "ref_valid_stop": "bool",
    "ref_first_stop_codon": "string",
    "ref_first_stop_pos": "int",
    "ref_stop_codon_count": "int",
    "ref_stop_codons": "stop_codon_list",
    "ref_stop_codon_exons": "int_list",
    "ref_has_ptc": "bool",
    "alt_last_codon": "string",
    "alt_valid_stop": "bool",
    "alt_first_stop_codon": "string",
    "alt_first_stop_pos": "int",
    "alt_stop_codon_count": "int",
    "alt_stop_codons": "stop_codon_list",
    "alt_stop_codon_exons": "int_list",
    "alt_has_ptc": "bool",
    "start_loss": "bool",
    "stop_loss": "bool",
    "stop_classification": "stop_classification",
    "transcript_start": "int",
    "transcript_end": "int",
    "transcript_seq": "string",
    "transcript_length": "int",
    "cds_start_in_transcript": "int",
    "cds_end_in_transcript": "int",
    "alt_transcript_seq": "string",
    "alt_transcript_length": "int",
    "alt_cds_start_in_transcript": "int",
    "transcript_exons": "pair_list",
    "alt_transcript_exons": "pair_list",
    "alt_scan_start_codon_pos": "int",
    "alt_scan_start_codon_exon": "int",
    "alt_scan_first_stop_codon": "string",
    "alt_scan_first_stop_pos": "int",
    "alt_scan_stop_codon_count": "int",
    "alt_scan_stop_codons": "stop_codon_list",
    "alt_scan_stop_codon_exons": "int_list",
    "unknown_reason": "unknown_reason",
}

# Kind of every column that add_nmd_features returns, in output order
NMD_FEATURE_COLUMN_KINDS = {
    "utr3_length": "int",
    "utr5_length": "int",
    "total_exon_count": "int",
    "ptc_pos_in_alt_transcript": "int",
    "ptc_exon_number": "int",
    "upstream_exon_count": "int",
    "downstream_exon_count": "int",
    "ptc_to_start_codon": "int",
    "ptc_less_than_150nt_to_start": "bool",
    "ptc_exon_length": "int",
    "annotated_stop_distance": "int",
    "ptc_to_exon_end": "int",
    "likely_misannotated": "bool",
}

# Kind of every column that evaluate_nmd_escape_rules returns, in output order
NMD_RULE_COLUMN_KINDS = {
    "nmd_last_exon_rule": "bool",
    "nmd_50nt_penultimate_rule": "bool",
    "nmd_long_exon_rule": "bool",
    "nmd_start_proximal_rule": "bool",
    "nmd_single_exon_rule": "bool",
    "nmd_escape": "bool",
}

# Kind of the column that add_features_and_rules adds after the rules. It says whether the NMD efficiency model
# can score the row (see MODEL_STATUSES).
MODEL_STATUS_COLUMN_KINDS = {"nmd_model_status": "nmd_model_status"}

# Kind of every output column, in output order. The kind gives the pandas dtype (KIND_DTYPES) and the
# Arrow and Parquet type (KIND_ARROW_TYPES, see to_arrow). Without a fixed schema, pandas and pyarrow infer each type from the
# data. A column with only missing values, or a table without rows, then gets a different type from run
# to run.
OUTPUT_COLUMN_KINDS = {
    **PTC_COLUMN_KINDS,
    **NMD_FEATURE_COLUMN_KINDS,
    **NMD_RULE_COLUMN_KINDS,
    **MODEL_STATUS_COLUMN_KINDS,
}

# The 4 columns that hold a sequence. They make up most of the table's size, in memory and on disk.
# annotate(..., sequences=False) and the CLI flag --no-sequences leave them out (see output_column_kinds).
SEQUENCE_COLUMNS = ("ref_cds_seq", "alt_cds_seq", "transcript_seq", "alt_transcript_seq")

# The columns of kind stop_codon_list. They hold lists of {"position": ..., "codon": ...} records, e.g.
# {"position": 5442, "codon": "TGA"}.
STOP_CODON_COLUMNS = ("ref_stop_codons", "alt_stop_codons", "alt_scan_stop_codons")

# The 19 inputs of the NMD efficiency model best_model.pkl (see scripts/train_new.ipynb), in the order that the
# model takes them. The model cannot score a row in which one of them is null.
MODEL_INPUTS = [
    "start_loss",
    "stop_loss",
    "total_exon_count",
    "ptc_less_than_150nt_to_start",
    "nmd_long_exon_rule",
    "nmd_start_proximal_rule",
    "nmd_single_exon_rule",
    "nmd_escape",
    "downstream_exon_count",
    "nmd_last_exon_rule",
    "ptc_to_start_codon",
    "annotated_stop_distance",
    "ptc_exon_length",
    "ptc_to_exon_end",
    "upstream_exon_count",
    "nmd_50nt_penultimate_rule",
    "utr5_length",
    "utr3_length",
    "transcript_length",
]


def apply_schema(table, column_kinds=OUTPUT_COLUMN_KINDS):
    """
    Return a copy of ``table`` with the columns of ``column_kinds``, in that order, and the dtype of each
    column's kind (see KIND_DTYPES). A missing value in an int, bool or string column becomes pd.NA, and one
    in a categorical column NaN. pandas raises if a value does not fit its column's dtype, e.g. 1.5 in an int
    column. A value outside the categories of a categorical column raises too.

    :param table: DataFrame with exactly the columns of ``column_kinds``, in any order
    :param column_kinds: dict of column name to kind, e.g. OUTPUT_COLUMN_KINDS, PTC_COLUMN_KINDS or
        output_column_kinds(sequences=False)
    :return: DataFrame with the schema of ``column_kinds``
    :raises ValueError: if ``table`` lacks a column of ``column_kinds`` or has a column it does not list, or if a
        categorical column holds a value outside its categories (see CATEGORIES)
    """

    missing = [column for column in column_kinds if column not in table.columns]
    unknown = [column for column in table.columns if column not in column_kinds]
    if missing or unknown:
        raise ValueError(
            f"The table does not match the schema. Missing columns: {missing}. Unknown columns: {unknown}."
        )
    # astype turns a value outside the categories into NaN, without an error
    for column, kind in column_kinds.items():
        if kind not in CATEGORIES:
            continue
        outside = set(table[column].dropna()) - set(CATEGORIES[kind])
        if outside:
            raise ValueError(
                f"The column {column} holds values outside its categories {CATEGORIES[kind]}: {sorted(outside, key=str)}."
            )

    return table[list(column_kinds)].astype({column: KIND_DTYPES[kind] for column, kind in column_kinds.items()})


def output_column_kinds(sequences=True):
    """
    Return the kind of every output column, in output order. With ``sequences=False``, the 4 columns of
    SEQUENCE_COLUMNS are left out, and the other 76 columns keep their order. Pass the result to
    ``apply_schema`` or ``empty_table``.

    :param sequences: whether to keep the columns of SEQUENCE_COLUMNS
    :return: a new dict of column name to kind: OUTPUT_COLUMN_KINDS, or OUTPUT_COLUMN_KINDS without SEQUENCE_COLUMNS
    """

    return {column: kind for column, kind in OUTPUT_COLUMN_KINDS.items() if sequences or column not in SEQUENCE_COLUMNS}


def empty_table(column_kinds=OUTPUT_COLUMN_KINDS):
    """
    Return a DataFrame without rows, with the columns and dtypes of ``column_kinds``.
    """

    return apply_schema(pd.DataFrame(columns=list(column_kinds)), column_kinds)


def to_arrow(results: pd.DataFrame) -> pa.Table:
    """
    Return ``results`` as a pyarrow Table, with the Arrow type of each column's kind (KIND_ARROW_TYPES).

    The types do not depend on the data. So they are the same for every input, also for zero rows or for a
    column with only missing values. A missing value becomes null. ``results`` is not changed. ``write_results``
    writes this table for a .parquet output, and ``pyarrow.parquet.write_table`` can write it too.

    :param results: DataFrame with output columns, e.g. from ``annotate``, with or without SEQUENCE_COLUMNS.
        Any subset of the output columns works, in any order.
    :return: pyarrow Table with the columns of ``results``, in that order, and without the index
    :raises KeyError: if ``results`` has a column that OUTPUT_COLUMN_KINDS does not list
    """

    return pa.Table.from_pandas(results, schema=_arrow_schema(results.columns), preserve_index=False)


def _arrow_schema(columns):
    """
    Return the pyarrow schema of ``columns``, with the Arrow type of each column's kind (KIND_ARROW_TYPES).
    A column that OUTPUT_COLUMN_KINDS does not list raises a KeyError.
    """

    return pa.schema([pa.field(column, KIND_ARROW_TYPES[OUTPUT_COLUMN_KINDS[column]]) for column in columns])
