"""
Schema of the result table: its columns, their order and their dtypes.

Each function that returns results gives them this schema, also for a table without rows. The section "Output
columns" of "Technical Notes.md" gives the meaning of each column and says when it is null.
"""

import pandas as pd

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
    "has_stop_codon": "bool",
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
    "transcript_exon_info": "pair_list",
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
    "stop_codon_distance": "int",
    "ptc_to_intron": "int",
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

# Kind of every output column, in output order. The kind gives the pandas dtype (KIND_DTYPES) and the
# parquet type (cli.parquet_schema). Without a fixed schema, pandas and pyarrow infer each type from the
# data. A column with only missing values, or a table without rows, then gets a different type from run
# to run.
OUTPUT_COLUMN_KINDS = {**PTC_COLUMN_KINDS, **NMD_FEATURE_COLUMN_KINDS, **NMD_RULE_COLUMN_KINDS}


def apply_schema(table, column_kinds=OUTPUT_COLUMN_KINDS):
    """
    Return a copy of ``table`` with the columns of ``column_kinds``, in that order, and the dtype of each
    column's kind (see KIND_DTYPES). A missing value in an int, bool or string column becomes pd.NA.
    pandas raises if a value does not fit its column's dtype, e.g. 1.5 in an int column.

    :param table: DataFrame with exactly the columns of ``column_kinds``, in any order
    :param column_kinds: dict of column name to kind, e.g. OUTPUT_COLUMN_KINDS or PTC_COLUMN_KINDS
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


def empty_table(column_kinds=OUTPUT_COLUMN_KINDS):
    """
    Return a DataFrame without rows, with the columns and dtypes of ``column_kinds``.
    """

    return apply_schema(pd.DataFrame(columns=list(column_kinds)), column_kinds)
