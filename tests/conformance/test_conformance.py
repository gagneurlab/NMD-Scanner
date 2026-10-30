"""
Run every conformance case, also without the sequence columns and through to_arrow, and check that the cases cover
every output column, every value of the bool and categorical columns, every documented null case and every reason for
no result.
Check that the drawing of each case holds its rendered layout block.
"""

import importlib
import re
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from nmd_scanner import to_arrow
from nmd_scanner.schema import CATEGORIES, OUTPUT_COLUMN_KINDS

from .runner import NO_RESULT_REASONS, NoResult, Raises, case_params, check, expected_pair, render_case, run

CASE_MODULES = {
    path.stem: importlib.import_module(f".{path.stem}", __package__)
    for path in sorted(Path(__file__).parent.glob("cases_*.py"))
}
CASES = [case for module in CASE_MODULES.values() for case in module.CASES]
PARAMS = case_params(CASES)


@pytest.mark.parametrize(("case", "change", "strand"), PARAMS)
def test_case(case, change, strand, tmp_path):
    check(case, change, strand, tmp_path)


@pytest.mark.parametrize("case", CASES, ids=[case.name for case in CASES])
def test_drawing_holds_its_layout_block(case):
    block = render_case(case)
    lines = [line.rstrip() for line in textwrap.dedent(case.drawing).splitlines()]
    wanted = block.splitlines()
    held = any(lines[i : i + len(wanted)] == wanted for i in range(len(lines)))
    assert held, f"the drawing does not hold its layout block, at the indentation of its other lines:\n{block}"


# The 4 sequence columns, which annotate(..., sequences=False) leaves out ("Technical Notes.md", "Output columns")
SEQUENCE_COLUMNS = ("ref_cds_seq", "alt_cds_seq", "transcript_seq", "alt_transcript_seq")


@pytest.mark.parametrize(
    "case",
    [pytest.param(case, marks=[pytest.mark.xfail(reason=case.bug, strict=True)] if case.bug else []) for case in CASES],
    ids=[case.name for case in CASES],
)
def test_case_without_sequences(case, tmp_path):
    """
    With sequences=False, each case gives its expected result without the 4 sequence columns: its variant-transcript
    pairs, with the values and dtypes of the other columns, or zero rows with these columns, or its error. Once per
    case: on the plus strand, with the first description of the variant.
    """
    columns = {column: kind for column, kind in OUTPUT_COLUMN_KINDS.items() if column not in SEQUENCE_COLUMNS}
    assert len(columns) == len(OUTPUT_COLUMN_KINDS) - 4

    check(case, case.change, "+", tmp_path, sequences=False, columns=columns)


# The Arrow type of each column kind: see `schema.KIND_ARROW_TYPES`
ARROW_TYPES = {
    "int": pa.int64(),
    "bool": pa.bool_(),
    "string": pa.string(),
    "pair_list": pa.list_(pa.struct([("exon_number", pa.int64()), ("length", pa.int64())])),
    "int_list": pa.list_(pa.int64()),
    "stop_codon_list": pa.list_(pa.struct([("position", pa.int64()), ("codon", pa.string())])),
    **dict.fromkeys(CATEGORIES, pa.dictionary(pa.int8(), pa.string())),
}
# The list columns: lists of {"exon_number", "length"} records, of exon numbers, or of {"position", "codon"} records
LIST_COLUMNS = tuple(
    column for column, kind in OUTPUT_COLUMN_KINDS.items() if kind in ("pair_list", "int_list", "stop_codon_list")
)
CASES_WITHOUT_AN_ERROR = [case for case in CASES if not isinstance(case.expected, Raises)]
CASES_WITH_A_RESULT = [case for case in CASES if isinstance(case.expected, dict)]


@pytest.mark.parametrize("sequences", [True, False], ids=["sequences", "no_sequences"])
@pytest.mark.parametrize("case", CASES_WITHOUT_AN_ERROR, ids=[case.name for case in CASES_WITHOUT_AN_ERROR])
def test_to_arrow_of_the_case_result(case, sequences, tmp_path):
    """
    to_arrow gives each column of the case result the Arrow type of its kind: with rows, without rows, in a column
    with only nulls, and without the sequence columns. The list columns hold the expected values, with the records as
    structs. to_arrow leaves the result unchanged. Once per case: on the plus strand, with the first
    description of the variant.
    """
    results, _ = run(case, case.change, "+", tmp_path, sequences=sequences)
    original = results.copy()

    table = to_arrow(results)

    columns = [column for column in OUTPUT_COLUMN_KINDS if sequences or column not in SEQUENCE_COLUMNS]
    assert [(field.name, field.type) for field in table.schema] == [
        (column, ARROW_TYPES[OUTPUT_COLUMN_KINDS[column]]) for column in columns
    ]
    pd.testing.assert_frame_equal(results, original)
    if isinstance(case.expected, NoResult):
        assert table.num_rows == 0
        return
    expected_pairs = [expected_pair(case, "+", more) for more in ({}, *case.more_pairs)]
    by_key = {(pair["transcript_id"], pair["variant_id"]): pair for pair in expected_pairs}
    assert table.num_rows == len(by_key) == len(expected_pairs)
    keys = zip(table.column("transcript_id").to_pylist(), table.column("variant_id").to_pylist())
    expected = [by_key[key] for key in keys]
    for column in LIST_COLUMNS:
        assert table.column(column).to_pylist() == [pair[column] for pair in expected], column


@pytest.mark.parametrize("case", CASES_WITH_A_RESULT, ids=[case.name for case in CASES_WITH_A_RESULT])
def test_list_columns_keep_their_shape_through_parquet(case, tmp_path):
    """
    A list value of the case result and the same value after to_arrow, Parquet and pd.read_parquet with default
    arguments have elements of the same shape: a record stays a dict with the same keys and values, and an exon
    number stays an int. Only the container changes, from a list to a numpy array, and an int becomes a numpy int. So
    code that reads the fields by name works on both. polars reads the records as structs with the same field names.
    Once per case: on the plus strand, with the first description of the variant.
    """
    results, _ = run(case, case.change, "+", tmp_path)
    path = tmp_path / "results.parquet"
    pq.write_table(to_arrow(results), path)

    loaded = pd.read_parquet(path)

    for column in LIST_COLUMNS:
        for before, after in zip(results[column], loaded[column], strict=True):
            if pd.api.types.is_scalar(before) and pd.isna(before):
                assert after is None, column
                continue
            assert isinstance(before, list) and isinstance(after, np.ndarray), column
            assert len(after) == len(before), column
            for element_before, element_after in zip(before, after):
                if isinstance(element_before, dict):
                    assert type(element_after) is dict and element_after == element_before, column
                else:
                    assert type(element_before) is int and isinstance(element_after, np.integer), column
                    assert element_after == element_before, column
    polars = pytest.importorskip("polars")
    schema = polars.read_parquet(path).schema
    assert schema["transcript_exons"] == polars.List(
        polars.Struct({"exon_number": polars.Int64, "length": polars.Int64})
    )
    assert schema["ref_stop_codon_exons"] == polars.List(polars.Int64)
    assert schema["ref_stop_codons"] == polars.List(polars.Struct({"position": polars.Int64, "codon": polars.String}))


def test_to_arrow_types_a_column_with_only_nulls(tmp_path):
    """
    A column of each Arrow type can hold only nulls. The one variant-transcript pair of a destroyed splice site has a
    null in the alt CDS columns, and the one pair of a missense SNV has a null unknown_reason.
    """
    null_columns = {
        "snv_at_donor_plus_1_destroys_the_splice_site": (
            *("alt_cds_length", "start_loss", "alt_cds_seq"),
            *("alt_cds_exons", "alt_stop_codon_exons", "alt_stop_codons"),
        ),
        "missense_snv_inside_an_internal_coding_exon_changes_one_cds_base": ("unknown_reason",),
    }
    arrow_types = {ARROW_TYPES[OUTPUT_COLUMN_KINDS[column]] for columns in null_columns.values() for column in columns}
    assert arrow_types == set(ARROW_TYPES.values())

    for name, columns in null_columns.items():
        case = next(case for case in CASES if case.name == name)
        directory = tmp_path / name
        directory.mkdir()
        table = to_arrow(run(case, case.change, "+", directory)[0])

        for column in columns:
            assert table.column(column).null_count == table.num_rows == 1, column
            assert table.schema.field(column).type == ARROW_TYPES[OUTPUT_COLUMN_KINDS[column]], column


def test_case_names_are_unique():
    names = [case.name for case in CASES]
    assert sorted(name for name in set(names) if names.count(name) > 1) == []


def _expected_pairs():
    """The expected variant-transcript pair of each case that gives a result, on each strand."""
    return [expected_pair(case, strand) for case in CASES if isinstance(case.expected, dict) for strand in "+-"]


# Every value that annotate() can give in a bool or categorical column
VALUES = {
    "strand": {"+", "-"},
    "has_start_codon": {True, False},
    "has_stop_codon": {True, False},
    "cds_frame": {0, 1, 2},
    "cds_in_transcript": {True, False},
    # None: a ref CDS of fewer than 3 nt
    "ref_valid_stop": {True, False, None},
    "ref_first_stop_codon": {"TAA", "TAG", "TGA", None},
    "ref_has_ptc": {True, False, None},
    "alt_valid_stop": {True, False, None},
    "alt_first_stop_codon": {"TAA", "TAG", "TGA", None},
    "alt_has_ptc": {True, False, None},
    "start_loss": {True, False, None},
    "stop_loss": {True, False, None},
    # None: a variant-transcript pair with unknown_reason
    "stop_classification": {"alt_transcript", "start_loss_scan", "alt_cds", None},
    "alt_scan_first_stop_codon": {"TAA", "TAG", "TGA", None},
    "unknown_reason": {"splice_site_destroyed", "exon_boundary_ambiguous", None},
    "ptc_less_than_150nt_to_start": {True, False, None},
    "likely_misannotated": {True, False},
    "nmd_last_exon_rule": {True, False, None},
    "nmd_50nt_penultimate_rule": {True, False, None},
    "nmd_long_exon_rule": {True, False, None},
    "nmd_start_proximal_rule": {True, False, None},
    "nmd_single_exon_rule": {True, False, None},
    "nmd_escape": {True, False, None},
    "nmd_efficiency_prediction_status": {
        "unknown_effect",
        "no_ptc",
        "ref_ptc",
        "no_annotated_stop",
        "no_annotated_start",
        "start_loss",
        "missing_input",
        "ok",
    },
}


@pytest.mark.parametrize("column", sorted(VALUES))
def test_every_value_has_a_case(column):
    assert VALUES[column] - {pair[column] for pair in _expected_pairs()} == set()


def test_values_of_a_categorical_column_are_its_categories():
    categorical = {column: kind for column, kind in OUTPUT_COLUMN_KINDS.items() if kind in CATEGORIES}
    assert {column: VALUES[column] - {None} for column in categorical} == {
        column: set(CATEGORIES[kind]) for column, kind in categorical.items()
    }


def test_every_output_column_has_a_case():
    # The runner checks every column of every case, so a column has a case if a case gives a result
    assert {column for pair in _expected_pairs() for column in pair} == set(OUTPUT_COLUMN_KINDS)


@pytest.mark.parametrize("reason", NO_RESULT_REASONS)
def test_every_reason_for_no_result_has_a_case(reason):
    assert any(isinstance(case.expected, NoResult) and case.expected.reason == reason for case in CASES)


# The null cases: one entry per "Null when" clause of the column tables in "Technical Notes.md" and input_defects.md.
# A predicate tells whether an expected variant-transcript pair shows the constellation of its clause, without
# another clause of the column that would make the column null too.


def _known(pair):
    return pair["unknown_reason"] is None


def _scanned(pair):
    """Whether the scan of the alt transcript runs: after a start or stop loss, on an alt transcript of 3 nt or more."""
    seq = pair["alt_transcript_seq"]
    return bool(pair["start_loss"] or pair["stop_loss"]) and seq is not None and len(seq) >= 3


def _no_exon_records(pair):
    """A variant-transcript pair shows that its transcript has no exon records by a null transcript_exons."""
    return pair["transcript_exons"] is None


def _keeps_the_flags_from_the_cds(pair):
    """Whether the variant-transcript pair keeps alt_has_ptc and stop_loss from the codon scan of the alt CDS."""
    if not _known(pair):
        return False
    seq = pair["alt_transcript_seq"]
    if seq is None or len(seq) < 3 or pair["cds_start_in_transcript"] is None:
        return True
    if pair["start_loss"] or not pair["has_stop_codon"]:
        return False
    # The ref transcript, read in frame from the first complete codon, does not stop at the annotated stop codon
    seq = pair["transcript_seq"]
    start = pair["cds_start_in_transcript"] + pair["cds_frame"]
    stops = [i for i in range(start, len(seq) - 2, 3) if seq[i : i + 3] in {"TAA", "TAG", "TGA"}]
    return stops[:1] != [pair["cds_end_in_transcript"] - 3]


def _ptc_variant(pair):
    return pair["alt_has_ptc"] is True


def _no_orf_overlaps_the_cds_after_a_start_loss(pair):
    """A start loss whose scan found no ATG, or one downstream of the first base of the annotated stop codon."""
    if not (_known(pair) and pair["has_stop_codon"] and pair["start_loss"] and _scanned(pair)):
        return False
    # The variant changes the start codon only, so the annotated stop codon is the last codon of the alt CDS
    stop = pair["alt_cds_start_in_transcript"] + pair["alt_cds_length"] - 3
    return pair["alt_scan_start_codon_pos"] is None or pair["alt_scan_start_codon_pos"] > stop


UNKNOWN = "`unknown_reason` is set"
NO_START_CODON = "`has_start_codon` is False"
SHORT_CDS = "the CDS has fewer than 3 nt"
NO_IN_FRAME_STOP = "no in-frame stop codon"
NO_EXON_RECORDS = "the transcript has no exon records"
NOT_SCANNED = "not scanned"
NOT_A_PTC_VARIANT = "not a PTC variant"
NO_ALT_EXONS = "`alt_transcript_exons` is null"
CDS_START_IS_NULL = "`cds_start_in_transcript` is null"
# The 5 NMD rules, which nmd_escape combines
RULES = (
    "nmd_last_exon_rule",
    "nmd_50nt_penultimate_rule",
    "nmd_long_exon_rule",
    "nmd_start_proximal_rule",
    "nmd_single_exon_rule",
)
NULL_CASES = [
    # A variant-transcript pair does not show the ID of its record, and a null variant_id has no other cause
    ("variant_id", "the ID of the VCF record is `.`", lambda pair: True),
    # transcript_type and transcript_tags come from the transcript record. If one of them is set, the transcript
    # has a transcript record.
    (
        "transcript_type",
        "the GENCODE transcript has no transcript record",
        lambda pair: pair["transcript_tags"] is None,
    ),
    (
        "transcript_type",
        "the transcript record has no `transcript_type` (GENCODE) or `biotype` (Ensembl)",
        lambda pair: pair["transcript_tags"] is not None,
    ),
    ("transcript_tags", "the transcript has no tag", lambda pair: pair["transcript_type"] is not None),
    (
        "transcript_tags",
        "the GENCODE transcript has no transcript record",
        lambda pair: pair["transcript_type"] is None,
    ),
    *[
        (column, UNKNOWN, lambda pair: not _known(pair))
        for column in (
            *("alt_cds_seq", "alt_cds_length", "alt_cds_exons", "alt_last_codon", "alt_valid_stop"),
            *("alt_first_stop_codon", "alt_first_stop_pos", "alt_stop_codon_count", "alt_stop_codons"),
            *("alt_stop_codon_exons", "alt_has_ptc", "start_loss", "stop_loss", "stop_classification"),
            *("alt_transcript_seq", "alt_transcript_length", "alt_cds_start_in_transcript"),
            *("alt_transcript_exons", "annotated_stop_distance"),
        )
    ],
    ("start_codon_exon_number", NO_START_CODON, lambda pair: not pair["has_start_codon"]),
    ("start_codon_exon_number", SHORT_CDS, lambda pair: pair["has_start_codon"] and pair["ref_cds_length"] < 3),
    *[
        (column, SHORT_CDS, lambda pair: pair["ref_cds_length"] < 3)
        for column in (
            *("ref_last_codon", "ref_valid_stop", "ref_first_stop_codon", "ref_first_stop_pos"),
            *("ref_stop_codon_count", "ref_stop_codons", "ref_stop_codon_exons", "ref_has_ptc"),
        )
    ],
    *[
        (column, NO_IN_FRAME_STOP, lambda pair: pair["ref_cds_length"] >= 3 and pair["ref_stop_codon_count"] == 0)
        for column in ("ref_first_stop_codon", "ref_first_stop_pos")
    ],
    *[
        (column, SHORT_CDS, lambda pair: _known(pair) and pair["alt_cds_length"] < 3)
        for column in (
            *("alt_last_codon", "alt_valid_stop", "alt_first_stop_codon", "alt_first_stop_pos"),
            *("alt_stop_codon_count", "alt_stop_codons", "alt_stop_codon_exons"),
        )
    ],
    *[
        (
            column,
            NO_IN_FRAME_STOP,
            lambda pair: _known(pair) and pair["alt_cds_length"] >= 3 and pair["alt_stop_codon_count"] == 0,
        )
        for column in ("alt_first_stop_codon", "alt_first_stop_pos")
    ],
    (
        "alt_has_ptc",
        "the variant-transcript pair keeps the flags from the CDS, and the alt CDS has fewer than 3 nt",
        lambda pair: _keeps_the_flags_from_the_cds(pair) and pair["alt_cds_length"] < 3,
    ),
    # A transcript without exon records gives a variant-transcript pair with null transcript columns
    *[
        (column, NO_EXON_RECORDS, _no_exon_records)
        for column in (
            *("transcript_start", "transcript_end", "transcript_seq", "transcript_length"),
            *("cds_start_in_transcript", "cds_end_in_transcript", "transcript_exons", "total_exon_count"),
        )
    ],
    *[
        entry
        for column in (
            *("alt_transcript_seq", "alt_transcript_length", "alt_cds_start_in_transcript"),
            "alt_transcript_exons",
        )
        for entry in [
            (column, CDS_START_IS_NULL, lambda pair: _known(pair) and pair["cds_start_in_transcript"] is None),
            (
                column,
                "`transcript_seq` does not hold `ref_cds_seq` at `cds_start_in_transcript`, or the ref bases of a "
                "UTR change next to it",
                lambda pair: (
                    _known(pair) and pair["cds_start_in_transcript"] is not None and pair["alt_transcript_seq"] is None
                ),
            ),
        ]
    ],
    (
        "alt_transcript_exons",
        "the exon lengths do not add up to `alt_transcript_length`",
        lambda pair: pair["alt_transcript_seq"] is not None,
    ),
    *[
        entry
        for column in ("alt_scan_start_codon_pos", "alt_scan_start_codon_exon_number")
        for entry in [
            (column, NOT_SCANNED, lambda pair: not _scanned(pair)),
            (
                column,
                "after a start loss: no ATG found",
                lambda pair: (
                    _scanned(pair)
                    and pair["start_loss"]
                    and "ATG"
                    not in pair["alt_transcript_seq"][pair["alt_cds_start_in_transcript"] + pair["cds_frame"] :]
                ),
            ),
            (
                column,
                "after a stop loss: `has_start_codon` is False",
                lambda pair: _scanned(pair) and not pair["start_loss"] and not pair["has_start_codon"],
            ),
        ]
    ],
    (
        "alt_scan_start_codon_exon_number",
        NO_ALT_EXONS,
        lambda pair: (
            _scanned(pair) and pair["alt_scan_start_codon_pos"] is not None and pair["alt_transcript_exons"] is None
        ),
    ),
    *[
        (column, NOT_SCANNED, lambda pair: not _scanned(pair))
        for column in (
            *("alt_scan_first_stop_codon", "alt_scan_first_stop_pos", "alt_scan_stop_codon_count"),
            "alt_scan_stop_codons",
            "alt_scan_stop_codon_exons",
        )
    ],
    (
        "alt_scan_stop_codon_exons",
        NO_ALT_EXONS,
        lambda pair: _scanned(pair) and pair["alt_transcript_exons"] is None,
    ),
    *[
        (column, "no stop codon found", lambda pair: _scanned(pair) and pair["alt_scan_stop_codon_count"] == 0)
        for column in ("alt_scan_first_stop_codon", "alt_scan_first_stop_pos")
    ],
    ("unknown_reason", "the alt transcript is known", lambda pair: pair["alt_cds_seq"] is not None),
    (
        "utr3_length",
        "`has_stop_codon` is False",
        lambda pair: not pair["has_stop_codon"] and pair["cds_end_in_transcript"] is not None,
    ),
    ("utr3_length", "`cds_end_in_transcript` is null", lambda pair: pair["cds_end_in_transcript"] is None),
    ("utr5_length", CDS_START_IS_NULL, lambda pair: pair["cds_start_in_transcript"] is None),
    *[
        entry
        for column in ("upstream_exon_count", "downstream_exon_count", "ptc_exon_length", "ptc_to_exon_end")
        for entry in [
            (column, NOT_A_PTC_VARIANT, lambda pair: pair["alt_has_ptc"] is False),
            (column, NO_ALT_EXONS, lambda pair: _ptc_variant(pair) and pair["alt_transcript_exons"] is None),
        ]
    ],
    ("ptc_pos_in_alt_transcript", NOT_A_PTC_VARIANT, lambda pair: pair["alt_has_ptc"] is False),
    (
        "ptc_pos_in_alt_transcript",
        "`alt_transcript_seq` is null",
        lambda pair: _ptc_variant(pair) and pair["alt_transcript_seq"] is None,
    ),
    ("ptc_exon_number", NOT_A_PTC_VARIANT, lambda pair: pair["alt_has_ptc"] is False),
    ("ptc_exon_number", NO_ALT_EXONS, lambda pair: _ptc_variant(pair) and pair["alt_transcript_exons"] is None),
    ("ptc_to_start_codon", NOT_A_PTC_VARIANT, lambda pair: pair["alt_has_ptc"] is False),
    (
        "ptc_to_start_codon",
        "the transcript has no annotated start codon",
        lambda pair: _ptc_variant(pair) and not pair["has_start_codon"],
    ),
    (
        "ptc_to_start_codon",
        "the annotated start codon is a stop codon, such as TAG",
        lambda pair: (
            _ptc_variant(pair)
            and pair["has_start_codon"]
            and pair["ref_cds_length"] >= 3
            and pair["ref_cds_seq"][:3] in {"TAA", "TAG", "TGA"}
        ),
    ),
    (
        "ptc_to_start_codon",
        "after a start loss: the variant-transcript pair has no `alt_transcript_seq`",
        lambda pair: _ptc_variant(pair) and pair["start_loss"] and pair["alt_transcript_seq"] is None,
    ),
    *[
        (column, NOT_A_PTC_VARIANT, lambda pair: pair["alt_has_ptc"] is False)
        for column in (
            *("ptc_less_than_150nt_to_start", "nmd_last_exon_rule", "nmd_50nt_penultimate_rule"),
            *("nmd_long_exon_rule", "nmd_start_proximal_rule", "nmd_single_exon_rule", "nmd_escape"),
        )
    ],
    # For a PTC variant, a rule is null if one of its inputs is null
    *[
        (
            column,
            f"`{rule_input}` is null",
            lambda pair, rule_input=rule_input: _ptc_variant(pair) and pair[rule_input] is None,
        )
        for column, rule_input in (
            ("ptc_less_than_150nt_to_start", "ptc_to_start_codon"),
            ("nmd_last_exon_rule", "downstream_exon_count"),
            ("nmd_50nt_penultimate_rule", "alt_transcript_exons"),
            ("nmd_long_exon_rule", "ptc_exon_length"),
            ("nmd_start_proximal_rule", "ptc_to_start_codon"),
            ("nmd_single_exon_rule", "total_exon_count"),
        )
    ],
    (
        "nmd_escape",
        "no rule is True, and one is null",
        lambda pair: (
            _ptc_variant(pair) and True not in [pair[rule] for rule in RULES] and None in [pair[rule] for rule in RULES]
        ),
    ),
    ("annotated_stop_distance", "`has_stop_codon` is False", lambda pair: _known(pair) and not pair["has_stop_codon"]),
    (
        "annotated_stop_distance",
        "a nonstop: the alt transcript has no in-frame stop codon",
        lambda pair: (
            _known(pair)
            and pair["has_stop_codon"]
            and not pair["start_loss"]
            and pair["stop_loss"]
            and not _keeps_the_flags_from_the_cds(pair)
            and pair["alt_scan_stop_codon_count"] == 0
        ),
    ),
    (
        "annotated_stop_distance",
        "after a start loss: the scan found no ATG, or one downstream of the annotated stop codon",
        _no_orf_overlaps_the_cds_after_a_start_loss,
    ),
    (
        "annotated_stop_distance",
        "on a variant-transcript pair that keeps the flags from the CDS, the alt CDS has no in-frame stop codon",
        lambda pair: (
            _keeps_the_flags_from_the_cds(pair) and pair["has_stop_codon"] and pair["alt_stop_codon_count"] in (0, None)
        ),
    ),
]
# The clauses that annotate() cannot reach, as (column, clause, reason)
CDS_RECORD_OUTSIDE_ITS_EXONS = (
    "read_gff3() raises a ValueError for a CDS record outside the exon records of its transcript"
)
UNREACHABLE_NULL_CASES = [
    ("cds_start_in_transcript", "the 5' CDS base lies outside the exons", CDS_RECORD_OUTSIDE_ITS_EXONS),
    ("cds_end_in_transcript", "the 5' CDS base lies outside the exons", CDS_RECORD_OUTSIDE_ITS_EXONS),
]


@pytest.mark.parametrize(
    ("column", "clause", "predicate"), NULL_CASES, ids=[f"{column}: {clause}" for column, clause, _ in NULL_CASES]
)
def test_every_null_case_has_a_case(column, clause, predicate):
    assert any(pair[column] is None and predicate(pair) for pair in _expected_pairs())


def _null_when(path, width):
    """{column: "Null when" cell} of the column table rows in the file at path, rows of width cells."""
    null_when = {}
    for line in path.read_text().splitlines():
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) == width and cells[0].strip("`") in OUTPUT_COLUMN_KINDS:
            null_when[cells[0].strip("`")] = cells[-1]
    return null_when


def _documented_null_clauses():
    """(column, clause) of each "Null when" clause in "Technical Notes.md" and input_defects.md."""
    null_when = _null_when(Path(__file__).parents[2] / "Technical Notes.md", 4)
    assert sorted(null_when) == sorted(OUTPUT_COLUMN_KINDS)
    defects = _null_when(Path(__file__).parent / "input_defects.md", 2)

    def clauses(column):
        # "as `x`" stands for the clauses of column x in both files
        found = []
        for cell in (null_when[column], defects.get(column, "never")):
            for clause in cell.split("; ") if cell != "never" else []:
                match = re.fullmatch(r"as `(\w+)`", clause)
                found += clauses(match[1]) if match else [clause]
        return found

    return [(column, clause) for column in null_when for clause in clauses(column)]


def test_null_cases_name_every_documented_clause_once():
    named = [(column, clause) for column, clause, _ in NULL_CASES + UNREACHABLE_NULL_CASES]
    assert sorted(named) == sorted(_documented_null_clauses())
