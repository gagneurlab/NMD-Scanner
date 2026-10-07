"""
Run every conformance case, also without the sequence columns and through to_arrow, and check that the cases cover
every output column, every value of the bool and categorical columns, every documented null case and every reason for
no row.
Check that the drawing of each case holds its rendered layout block.
"""

import importlib
import re
import textwrap
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pytest

from nmd_scanner import to_arrow
from nmd_scanner.schema import OUTPUT_COLUMN_KINDS

from .runner import NO_ROW_REASONS, NoRow, Raises, case_params, check, expected_row, render_case, run

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
    With sequences=False, each case gives its expected result without the 4 sequence columns: its rows, with the
    values and dtypes of the other columns, or no row with these columns, or its error. Once per case: on the plus
    strand, with the first description of the variant.
    """
    columns = {column: kind for column, kind in OUTPUT_COLUMN_KINDS.items() if column not in SEQUENCE_COLUMNS}
    assert len(columns) == len(OUTPUT_COLUMN_KINDS) - 4

    check(case, case.change, "+", tmp_path, sequences=False, columns=columns)


# The Arrow type of each column kind: its Parquet type in "Technical Notes.md" ("Kinds and dtypes")
ARROW_TYPES = {
    "int": pa.int64(),
    "bool": pa.bool_(),
    "string": pa.string(),
    "pair_list": pa.list_(pa.list_(pa.int64())),
    "int_list": pa.list_(pa.int64()),
    "stop_codon_list": pa.list_(pa.struct([("position", pa.int64()), ("codon", pa.string())])),
}
# The columns of kind stop_codon_list: lists of (position, codon) tuples, which to_arrow turns into records
STOP_CODON_COLUMNS = ("ref_all_stop_codons", "alt_all_stop_codons", "transcript_all_stop_codons")
CASES_WITH_A_RESULT = [case for case in CASES if not isinstance(case.expected, Raises)]


@pytest.mark.parametrize("sequences", [True, False], ids=["sequences", "no_sequences"])
@pytest.mark.parametrize("case", CASES_WITH_A_RESULT, ids=[case.name for case in CASES_WITH_A_RESULT])
def test_to_arrow_of_the_case_result(case, sequences, tmp_path):
    """
    to_arrow gives each column of the case result the Arrow type of its kind: with rows, without rows, in a column
    with only nulls, and without the sequence columns. The stop codon columns hold {"position", "codon"} records with
    the expected values. to_arrow leaves the result unchanged. Once per case: on the plus strand, with the first
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
    if isinstance(case.expected, NoRow):
        assert table.num_rows == 0
        return
    rows = [expected_row(case, "+", more) for more in ({}, *case.more_rows)]
    by_key = {(row["transcript_id"], row["variant_id"]): row for row in rows}
    assert table.num_rows == len(by_key) == len(rows)
    keys = zip(table.column("transcript_id").to_pylist(), table.column("variant_id").to_pylist())
    expected = [by_key[key] for key in keys]
    for column in STOP_CODON_COLUMNS:
        records = [
            None if row[column] is None else [{"position": position, "codon": codon} for position, codon in row[column]]
            for row in expected
        ]
        assert table.column(column).to_pylist() == records, column


def test_to_arrow_types_a_column_with_only_nulls(tmp_path):
    """The one row of a destroyed splice site has a null in a column of each kind: the alt CDS columns."""
    case = next(case for case in CASES if case.name == "snv_at_donor_plus_1_destroys_the_splice_site")
    null_columns = (
        *("alt_cds_len", "start_loss", "alt_cds_seq"),
        *("alt_cds_info", "alt_stop_codon_exons", "alt_all_stop_codons"),
    )
    assert {OUTPUT_COLUMN_KINDS[column] for column in null_columns} == set(ARROW_TYPES)

    table = to_arrow(run(case, case.change, "+", tmp_path)[0])

    for column in null_columns:
        assert table.column(column).null_count == table.num_rows == 1, column
        assert table.schema.field(column).type == ARROW_TYPES[OUTPUT_COLUMN_KINDS[column]], column


def test_case_names_are_unique():
    names = [case.name for case in CASES]
    assert sorted(name for name in set(names) if names.count(name) > 1) == []


def _rows():
    """The expected row of each case that gives a row, on each strand."""
    return [expected_row(case, strand) for case in CASES if isinstance(case.expected, dict) for strand in "+-"]


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
    "ref_is_premature": {True, False, None},
    "alt_valid_stop": {True, False, None},
    "alt_first_stop_codon": {"TAA", "TAG", "TGA", None},
    "alt_is_premature": {True, False, None},
    "start_loss": {True, False, None},
    "stop_loss": {True, False, None},
    "transcript_valid_stop": {True, False, None},
    "transcript_first_stop_codon": {"TAA", "TAG", "TGA", None},
    "unknown_reason": {"splice_site_destroyed", "exon_boundary_ambiguous", None},
    "ptc_less_than_150nt_to_start": {True, False, None},
    "likely_misannotated": {True, False},
    "nmd_last_exon_rule": {True, False, None},
    "nmd_50nt_penultimate_rule": {True, False, None},
    "nmd_long_exon_rule": {True, False, None},
    "nmd_start_proximal_rule": {True, False, None},
    "nmd_single_exon_rule": {True, False, None},
    "nmd_escape": {True, False, None},
    "nmd_model_status": {
        "unknown_effect",
        "no_ptc",
        "ref_ptc",
        "no_annotated_stop",
        "no_annotated_start",
        "start_lost",
        "missing_input",
        "ok",
    },
}


@pytest.mark.parametrize("column", sorted(VALUES))
def test_every_value_has_a_case(column):
    assert VALUES[column] - {row[column] for row in _rows()} == set()


def test_every_output_column_has_a_case():
    # The runner checks every column of every case, so a column has a case if a case gives a row
    assert {column for row in _rows() for column in row} == set(OUTPUT_COLUMN_KINDS)


@pytest.mark.parametrize("reason", NO_ROW_REASONS)
def test_every_reason_for_no_row_has_a_case(reason):
    assert any(isinstance(case.expected, NoRow) and case.expected.reason == reason for case in CASES)


# The null cases: one entry per "Null when" clause of the column tables in "Technical Notes.md" and "Input Defects.md".
# A predicate tells whether an expected row shows the constellation of its clause, without another clause of the
# column that would make the column null too.


def _known(row):
    return row["unknown_reason"] is None


def _scanned(row):
    """Whether the scan of the alt transcript runs: after a start or stop loss, on an alt transcript of 3 nt or more."""
    seq = row["alt_transcript_seq"]
    return bool(row["start_loss"] or row["stop_loss"]) and seq is not None and len(seq) >= 3


def _no_exon_rows(row):
    """A row shows that its transcript has no exon rows by a null transcript_exon_info."""
    return row["transcript_exon_info"] is None


def _keeps_the_flags_from_the_cds(row):
    """Whether the row keeps alt_is_premature and stop_loss from the codon scan of the alt CDS."""
    if not _known(row):
        return False
    seq = row["alt_transcript_seq"]
    if seq is None or len(seq) < 3 or row["cds_start_in_transcript"] is None:
        return True
    if row["start_loss"] or not row["has_stop_codon"]:
        return False
    # The ref transcript, read in frame from the first complete codon, does not stop at the annotated stop codon
    seq = row["transcript_seq"]
    start = row["cds_start_in_transcript"] + row["cds_frame"]
    stops = [i for i in range(start, len(seq) - 2, 3) if seq[i : i + 3] in {"TAA", "TAG", "TGA"}]
    return stops[:1] != [row["cds_end_in_transcript"] - 3]


def _ptc_row(row):
    return row["alt_is_premature"] is True


def _no_orf_overlaps_the_cds_after_a_start_loss(row):
    """A start loss whose scan found no ATG, or one downstream of the first base of the annotated stop codon."""
    if not (_known(row) and row["has_stop_codon"] and row["start_loss"] and _scanned(row)):
        return False
    # The variant changes the start codon only, so the annotated stop codon is the last codon of the alt CDS
    stop = row["alt_cds_start_in_transcript"] + row["alt_cds_len"] - 3
    return row["transcript_start_codon_pos"] is None or row["transcript_start_codon_pos"] > stop


UNKNOWN = "`unknown_reason` is set"
NO_START_CODON = "`has_start_codon` is False"
SHORT_CDS = "the CDS has fewer than 3 nt"
NO_IN_FRAME_STOP = "no in-frame stop codon"
NO_EXON_ROWS = "the transcript has no exon rows"
NOT_SCANNED = "not scanned"
NOT_A_PTC_ROW = "not a PTC row"
NO_ALT_EXONS = "`alt_transcript_exon_info` is null"
CDS_START_IS_NULL = "`cds_start_in_transcript` is null"
NULL_CASES = [
    *[
        (column, UNKNOWN, lambda row: not _known(row))
        for column in (
            *("alt_cds_start", "alt_cds_stop", "alt_cds_seq", "alt_cds_len", "alt_cds_info"),
            *("alt_start_codon_pos", "alt_start_codon_exon", "alt_last_codon", "alt_valid_stop"),
            *("alt_first_stop_codon", "alt_first_stop_pos", "alt_num_stop_codons", "alt_all_stop_codons"),
            *("alt_stop_codon_exons", "alt_is_premature", "start_loss", "stop_loss"),
            *("alt_transcript_seq", "alt_transcript_length", "alt_cds_start_in_transcript"),
            "alt_transcript_exon_info",
            *("ptc_less_than_150nt_to_start", "annotated_stop_distance"),
            *("nmd_last_exon_rule", "nmd_50nt_penultimate_rule", "nmd_long_exon_rule", "nmd_start_proximal_rule"),
            *("nmd_single_exon_rule", "nmd_escape"),
        )
    ],
    *[
        entry
        for column in ("ref_start_codon_pos", "ref_start_codon_exon")
        for entry in [
            (column, NO_START_CODON, lambda row: not row["has_start_codon"]),
            (column, SHORT_CDS, lambda row: row["has_start_codon"] and row["ref_cds_len"] < 3),
        ]
    ],
    *[
        (column, SHORT_CDS, lambda row: row["ref_cds_len"] < 3)
        for column in (
            *("ref_last_codon", "ref_valid_stop", "ref_first_stop_codon", "ref_first_stop_pos"),
            *("ref_num_stop_codons", "ref_all_stop_codons", "ref_stop_codon_exons", "ref_is_premature"),
        )
    ],
    *[
        (column, NO_IN_FRAME_STOP, lambda row: row["ref_cds_len"] >= 3 and row["ref_num_stop_codons"] == 0)
        for column in ("ref_first_stop_codon", "ref_first_stop_pos")
    ],
    *[
        entry
        for column in ("alt_start_codon_pos", "alt_start_codon_exon")
        for entry in [
            (column, NO_START_CODON, lambda row: _known(row) and not row["has_start_codon"]),
            (
                column,
                SHORT_CDS,
                lambda row: _known(row) and row["has_start_codon"] and row["alt_cds_len"] < 3 and not row["start_loss"],
            ),
            (
                column,
                "the variant changes the start codon",
                lambda row: _known(row) and row["start_loss"] and row["alt_cds_len"] >= 3,
            ),
        ]
    ],
    *[
        (column, SHORT_CDS, lambda row: _known(row) and row["alt_cds_len"] < 3)
        for column in (
            *("alt_last_codon", "alt_valid_stop", "alt_first_stop_codon", "alt_first_stop_pos"),
            *("alt_num_stop_codons", "alt_all_stop_codons", "alt_stop_codon_exons"),
        )
    ],
    *[
        (
            column,
            NO_IN_FRAME_STOP,
            lambda row: _known(row) and row["alt_cds_len"] >= 3 and row["alt_num_stop_codons"] == 0,
        )
        for column in ("alt_first_stop_codon", "alt_first_stop_pos")
    ],
    (
        "alt_is_premature",
        "the row keeps the flags from the CDS, and the alt CDS has fewer than 3 nt",
        lambda row: _keeps_the_flags_from_the_cds(row) and row["alt_cds_len"] < 3,
    ),
    # A transcript without exon rows gives a row with null transcript columns (NU-08)
    *[
        (column, NO_EXON_ROWS, _no_exon_rows)
        for column in (
            *("transcript_start", "transcript_end", "transcript_seq", "transcript_length"),
            *("cds_start_in_transcript", "cds_end_in_transcript", "transcript_exon_info", "total_exon_count"),
        )
    ],
    *[
        entry
        for column in (
            *("alt_transcript_seq", "alt_transcript_length", "alt_cds_start_in_transcript"),
            "alt_transcript_exon_info",
        )
        for entry in [
            (column, CDS_START_IS_NULL, lambda row: _known(row) and row["cds_start_in_transcript"] is None),
            (
                column,
                "`transcript_seq` does not hold `ref_cds_seq` at `cds_start_in_transcript`, or the ref bases of a "
                "UTR change next to it",
                lambda row: (
                    _known(row) and row["cds_start_in_transcript"] is not None and row["alt_transcript_seq"] is None
                ),
            ),
        ]
    ],
    (
        "alt_transcript_exon_info",
        "the exon lengths do not add up to `alt_transcript_length`",
        lambda row: row["alt_transcript_seq"] is not None,
    ),
    *[
        entry
        for column in ("transcript_start_codon_pos", "transcript_start_codon_exon")
        for entry in [
            (column, NOT_SCANNED, lambda row: not _scanned(row)),
            (
                column,
                "after a start loss: no ATG found",
                lambda row: (
                    _scanned(row)
                    and row["start_loss"]
                    and "ATG" not in row["alt_transcript_seq"][row["alt_cds_start_in_transcript"] + row["cds_frame"] :]
                ),
            ),
            (
                column,
                "after a stop loss: `has_start_codon` is False",
                lambda row: _scanned(row) and not row["start_loss"] and not row["has_start_codon"],
            ),
        ]
    ],
    (
        "transcript_start_codon_exon",
        NO_ALT_EXONS,
        lambda row: (
            _scanned(row) and row["transcript_start_codon_pos"] is not None and row["alt_transcript_exon_info"] is None
        ),
    ),
    *[
        (column, NOT_SCANNED, lambda row: not _scanned(row))
        for column in (
            *("transcript_last_codon", "transcript_valid_stop", "transcript_first_stop_codon"),
            *("transcript_first_stop_pos", "transcript_num_stop_codons", "transcript_all_stop_codons"),
            "transcript_stop_codon_exons",
        )
    ],
    (
        "transcript_stop_codon_exons",
        NO_ALT_EXONS,
        lambda row: _scanned(row) and row["alt_transcript_exon_info"] is None,
    ),
    *[
        (column, "no stop codon found", lambda row: _scanned(row) and row["transcript_num_stop_codons"] == 0)
        for column in ("transcript_first_stop_codon", "transcript_first_stop_pos")
    ],
    ("unknown_reason", "the alt transcript is known", lambda row: row["alt_cds_seq"] is not None),
    (
        "utr3_length",
        "`has_stop_codon` is False",
        lambda row: not row["has_stop_codon"] and row["cds_end_in_transcript"] is not None,
    ),
    ("utr3_length", "`cds_end_in_transcript` is null", lambda row: row["cds_end_in_transcript"] is None),
    ("utr5_length", CDS_START_IS_NULL, lambda row: row["cds_start_in_transcript"] is None),
    *[
        entry
        for column in ("upstream_exon_count", "downstream_exon_count", "ptc_exon_length", "ptc_to_exon_end")
        for entry in [
            (column, NOT_A_PTC_ROW, lambda row: row["alt_is_premature"] is False),
            (column, NO_ALT_EXONS, lambda row: _ptc_row(row) and row["alt_transcript_exon_info"] is None),
        ]
    ],
    ("ptc_to_start_codon", NOT_A_PTC_ROW, lambda row: row["alt_is_premature"] is False),
    (
        "ptc_to_start_codon",
        "the transcript has no annotated start codon",
        lambda row: _ptc_row(row) and not row["has_start_codon"],
    ),
    (
        "ptc_to_start_codon",
        "the annotated start codon is a stop codon, such as TAG",
        lambda row: (
            _ptc_row(row) and row["ref_start_codon_pos"] == 0 and row["ref_cds_seq"][:3] in {"TAA", "TAG", "TGA"}
        ),
    ),
    (
        "ptc_to_start_codon",
        "after a start loss: the row has no `alt_transcript_seq`",
        lambda row: _ptc_row(row) and row["start_loss"] and row["alt_transcript_seq"] is None,
    ),
    ("annotated_stop_distance", "`has_stop_codon` is False", lambda row: _known(row) and not row["has_stop_codon"]),
    (
        "annotated_stop_distance",
        "a nonstop: the alt transcript has no in-frame stop codon",
        lambda row: (
            _known(row)
            and row["has_stop_codon"]
            and not row["start_loss"]
            and row["stop_loss"]
            and not _keeps_the_flags_from_the_cds(row)
            and row["transcript_num_stop_codons"] == 0
        ),
    ),
    (
        "annotated_stop_distance",
        "after a start loss: the scan found no ATG, or one downstream of the annotated stop codon",
        _no_orf_overlaps_the_cds_after_a_start_loss,
    ),
    (
        "annotated_stop_distance",
        "on a row that keeps the flags from the CDS, the alt CDS has no in-frame stop codon",
        lambda row: (
            _keeps_the_flags_from_the_cds(row) and row["has_stop_codon"] and row["alt_num_stop_codons"] in (0, None)
        ),
    ),
]
# The clauses that annotate() cannot reach, as (column, clause, reason)
CDS_ROW_OUTSIDE_ITS_EXONS = "read_gff3() raises a ValueError for a CDS row outside the exon rows of its transcript"
UNREACHABLE_NULL_CASES = [
    ("cds_start_in_transcript", "the 5' CDS base lies outside the exons", CDS_ROW_OUTSIDE_ITS_EXONS),
    ("cds_end_in_transcript", "the 5' CDS base lies outside the exons", CDS_ROW_OUTSIDE_ITS_EXONS),
]


@pytest.mark.parametrize(
    ("column", "clause", "predicate"), NULL_CASES, ids=[f"{column}: {clause}" for column, clause, _ in NULL_CASES]
)
def test_every_null_case_has_a_case(column, clause, predicate):
    assert any(row[column] is None and predicate(row) for row in _rows())


def _null_when(name, width):
    """{column: "Null when" cell} of the column table rows in the file name, rows of width cells."""
    null_when = {}
    for line in (Path(__file__).parents[2] / name).read_text().splitlines():
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) == width and cells[0].strip("`") in OUTPUT_COLUMN_KINDS:
            null_when[cells[0].strip("`")] = cells[-1]
    return null_when


def _documented_null_clauses():
    """(column, clause) of each "Null when" clause in "Technical Notes.md" and "Input Defects.md"."""
    null_when = _null_when("Technical Notes.md", 4)
    assert sorted(null_when) == sorted(OUTPUT_COLUMN_KINDS)
    defects = _null_when("Input Defects.md", 2)

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
