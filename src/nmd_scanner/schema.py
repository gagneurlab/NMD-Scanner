"""
Schema of the result table: its columns, their order and their dtypes.

Each function that returns results gives them this schema, also for a table without rows. The section "Output
columns" of "Technical Notes.md" gives the meaning of each column and says when it is null. ``to_arrow`` converts
a result to a pyarrow Table with the same types for every input.
"""

import pandas as pd
import pyarrow as pa

# pandas dtype of each column kind. The int, bool and string dtypes are nullable, so a missing value is
# pd.NA, and a column with missing values keeps its dtype. The string dtype has "python" storage: it
# behaves the same on pandas 2 and 3, with or without pyarrow. The default string dtype of pandas 3,
# "str", marks a missing value as NaN instead. The list kinds hold Python lists in object columns:
# "pair_list" holds (exon_number, length) tuples, "int_list" exon numbers and "stop_codon_list"
# (position, codon) tuples.
KIND_DTYPES = {
    "string": pd.StringDtype("python"),
    "int": pd.Int64Dtype(),
    "bool": pd.BooleanDtype(),
    "pair_list": object,
    "int_list": object,
    "stop_codon_list": object,
}

# Arrow type of each column kind, which to_arrow and the Parquet output use. A pair_list becomes a list of
# [exon_number, length] lists, and a stop_codon_list a list of {"position": ..., "codon": ...} structs.
KIND_ARROW_TYPES = {
    "string": pa.string(),
    "int": pa.int64(),
    "bool": pa.bool_(),
    "pair_list": pa.list_(pa.list_(pa.int64())),
    "int_list": pa.list_(pa.int64()),
    "stop_codon_list": pa.list_(pa.struct([pa.field("position", pa.int64()), pa.field("codon", pa.string())])),
}

# Kind of every column that extract_ptc returns, in output order
PTC_COLUMN_KINDS = {
    "transcript_id": "string",
    "variant_id": "string",
    "ref_cds_start": "int",
    "ref_cds_stop": "int",
    "ref_cds_seq": "string",
    "ref_cds_len": "int",
    "alt_cds_start": "int",
    "alt_cds_stop": "int",
    "alt_cds_seq": "string",
    "alt_cds_len": "int",
    "chromosome": "string",
    "gene_id": "string",
    "strand": "string",
    "has_start_codon": "bool",
    "has_stop_codon": "bool",
    "cds_frame": "int",
    "ref": "string",
    "alt": "string",
    "start_variant": "int",
    "end_variant": "int",
    "ref_cds_info": "pair_list",
    "alt_cds_info": "pair_list",
    "cds_in_transcript": "bool",
    "ref_start_codon_pos": "int",
    "ref_start_codon_exon": "int",
    "ref_last_codon": "string",
    "ref_valid_stop": "bool",
    "ref_first_stop_codon": "string",
    "ref_first_stop_pos": "int",
    "ref_num_stop_codons": "int",
    "ref_all_stop_codons": "stop_codon_list",
    "ref_stop_codon_exons": "int_list",
    "ref_is_premature": "bool",
    "alt_start_codon_pos": "int",
    "alt_start_codon_exon": "int",
    "alt_last_codon": "string",
    "alt_valid_stop": "bool",
    "alt_first_stop_codon": "string",
    "alt_first_stop_pos": "int",
    "alt_num_stop_codons": "int",
    "alt_all_stop_codons": "stop_codon_list",
    "alt_stop_codon_exons": "int_list",
    "alt_is_premature": "bool",
    "start_loss": "bool",
    "stop_loss": "bool",
    "transcript_start": "int",
    "transcript_end": "int",
    "transcript_seq": "string",
    "transcript_length": "int",
    "cds_start_in_transcript": "int",
    "cds_end_in_transcript": "int",
    "alt_transcript_seq": "string",
    "alt_transcript_length": "int",
    "alt_cds_start_in_transcript": "int",
    "transcript_exon_info": "pair_list",
    "alt_transcript_exon_info": "pair_list",
    "transcript_start_codon_pos": "int",
    "transcript_start_codon_exon": "int",
    "transcript_last_codon": "string",
    "transcript_valid_stop": "bool",
    "transcript_first_stop_codon": "string",
    "transcript_first_stop_pos": "int",
    "transcript_num_stop_codons": "int",
    "transcript_all_stop_codons": "stop_codon_list",
    "transcript_stop_codon_exons": "int_list",
    "unknown_reason": "string",
}

# Kind of every column that add_nmd_features returns, in output order
NMD_FEATURE_COLUMN_KINDS = {
    "utr3_length": "int",
    "utr5_length": "int",
    "total_exon_count": "int",
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
MODEL_STATUS_COLUMN_KINDS = {"nmd_model_status": "string"}

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

# The columns of kind stop_codon_list. They hold lists of (position, codon) tuples, e.g. (5442, "TGA").
# pyarrow's pandas conversion treats each tuple as a flat list of one type: it infers int from the first
# field and then fails on the string. to_arrow turns the tuples into {"position": ..., "codon": ...}
# records instead, which fit the struct of KIND_ARROW_TYPES["stop_codon_list"].
STOP_CODON_COLUMNS = ("ref_all_stop_codons", "alt_all_stop_codons", "transcript_all_stop_codons")

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

# The values of nmd_model_status, in the order they are checked. A row gets the first value whose condition holds:
# - unknown_effect: unknown_reason is set. The alt transcript is unknown, so alt_is_premature and 15 of the 19 model
#   inputs are null. unknown_reason says why.
# - no_ptc: alt_is_premature is not True, so there is no PTC to score.
# - ref_ptc: ref_is_premature is True. The reference has a PTC already, so the variant does not create it.
# - no_annotated_stop: has_stop_codon is False, which makes annotated_stop_distance and utr3_length null.
# - no_annotated_start: has_start_codon is False, which makes ptc_to_start_codon null. The true start codon lies
#   upstream of the CDS, at an unknown distance (e.g. cds_start_NF).
# - start_lost: start_loss is True and ptc_to_start_codon is null. After a start loss, the scan of alt_transcript_seq
#   takes the next ATG, and ptc_to_start_codon runs from it to the PTC. A scanned row is a PTC row only if the scan
#   finds an ATG and a stop codon after it, so it always has a ptc_to_start_codon. Only a start-loss row without
#   alt_transcript_seq gets start_lost, e.g. one of a transcript without exon rows. It is not scanned and keeps the
#   PTC of the alt CDS.
# - missing_input: another model input is null.
# - ok: the variant creates the PTC, and no model input is null.
# no_annotated_stop and no_annotated_start hold for every variant of the transcript. So they go before start_lost and
# missing_input, which depend on the variant.
MODEL_STATUSES = (
    "unknown_effect",
    "no_ptc",
    "ref_ptc",
    "no_annotated_stop",
    "no_annotated_start",
    "start_lost",
    "missing_input",
    "ok",
)


def apply_schema(table, column_kinds=OUTPUT_COLUMN_KINDS):
    """
    Return a copy of ``table`` with the columns of ``column_kinds``, in that order, and the dtype of each
    column's kind (see KIND_DTYPES). A missing value in an int, bool or string column becomes pd.NA.
    pandas raises if a value does not fit its column's dtype, e.g. 1.5 in an int column.

    :param table: DataFrame with exactly the columns of ``column_kinds``, in any order
    :param column_kinds: dict of column name to kind, e.g. OUTPUT_COLUMN_KINDS, PTC_COLUMN_KINDS or
        output_column_kinds(sequences=False)
    :return: DataFrame with the schema of ``column_kinds``
    :raises ValueError: if ``table`` lacks a column of ``column_kinds`` or has a column it does not list
    """

    missing = [column for column in column_kinds if column not in table.columns]
    unknown = [column for column in table.columns if column not in column_kinds]
    if missing or unknown:
        raise ValueError(
            f"The table does not match the schema. Missing columns: {missing}. Unknown columns: {unknown}."
        )

    return table[list(column_kinds)].astype({column: KIND_DTYPES[kind] for column, kind in column_kinds.items()})


def output_column_kinds(sequences=True):
    """
    Return the kind of every output column, in output order. With ``sequences=False``, the 4 columns of
    SEQUENCE_COLUMNS are left out, and the other 80 columns keep their order. Pass the result to
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
    column with only missing values. A missing value becomes null. The stop codon columns (STOP_CODON_COLUMNS)
    become lists of {"position": ..., "codon": ...} records. ``results`` is not changed. ``write_results``
    writes this table for a .parquet output, and ``pyarrow.parquet.write_table`` can write it too.

    :param results: DataFrame with output columns, e.g. from ``annotate``, with or without SEQUENCE_COLUMNS.
        Any subset of the output columns works, in any order.
    :return: pyarrow Table with the columns of ``results``, in that order, and without the index
    :raises KeyError: if ``results`` has a column that OUTPUT_COLUMN_KINDS does not list
    """

    return pa.Table.from_pandas(
        _stop_codon_records(results), schema=_arrow_schema(results.columns), preserve_index=False
    )


def _arrow_schema(columns):
    """
    Return the pyarrow schema of ``columns``, with the Arrow type of each column's kind (KIND_ARROW_TYPES).
    A column that OUTPUT_COLUMN_KINDS does not list raises a KeyError.
    """

    return pa.schema([pa.field(column, KIND_ARROW_TYPES[OUTPUT_COLUMN_KINDS[column]]) for column in columns])


def _stop_codon_records(results):
    """
    Return ``results`` with the (position, codon) tuples of STOP_CODON_COLUMNS turned into
    {"position": ..., "codon": ...} records. Without a stop codon column, return ``results`` itself.
    Otherwise return a copy, and leave ``results`` unchanged.
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
