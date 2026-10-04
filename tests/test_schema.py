import numpy as np
import pandas as pd
import pytest

from nmd_scanner import cli
from nmd_scanner.extra_features import add_nmd_features, evaluate_nmd_escape_rules
from nmd_scanner.rules import extract_ptc
from nmd_scanner.schema import (
    KIND_DTYPES,
    NMD_FEATURE_COLUMN_KINDS,
    NMD_RULE_COLUMN_KINDS,
    OUTPUT_COLUMN_KINDS,
    PTC_COLUMN_KINDS,
    apply_schema,
    empty_table,
)


def assert_schema(table, column_kinds):
    """Assert that ``table`` has the columns of ``column_kinds``, in that order, and their dtypes."""

    assert list(table.columns) == list(column_kinds)
    expected = {column: pd.api.types.pandas_dtype(KIND_DTYPES[kind]) for column, kind in column_kinds.items()}
    assert dict(table.dtypes) == expected


def test_every_kind_has_a_dtype():
    assert set(OUTPUT_COLUMN_KINDS.values()) <= set(KIND_DTYPES)


def test_output_columns_are_the_ptc_feature_and_rule_columns_in_order():
    parts = [*PTC_COLUMN_KINDS, *NMD_FEATURE_COLUMN_KINDS, *NMD_RULE_COLUMN_KINDS]
    assert list(OUTPUT_COLUMN_KINDS) == parts
    assert len(set(parts)) == len(parts)


def test_cli_output_column_kinds_is_the_schema():
    assert cli.OUTPUT_COLUMN_KINDS is OUTPUT_COLUMN_KINDS


def test_apply_schema_orders_columns_and_sets_dtypes():
    column_kinds = {"name": "string", "count": "int", "flag": "bool", "exons": "int_list"}
    table = pd.DataFrame(
        {
            "exons": [[1, 2], None],
            "flag": [True, None],
            "count": [3.0, np.nan],
            "name": ["a", None],
        }
    )

    result = apply_schema(table, column_kinds)

    assert_schema(result, column_kinds)
    assert result["count"].tolist() == [3, pd.NA]
    assert result["flag"].tolist() == [True, pd.NA]
    assert result["name"].tolist() == ["a", pd.NA]
    assert result["exons"].tolist() == [[1, 2], None]


@pytest.mark.parametrize(
    "columns, message",
    [(["name"], "Missing columns: \\['count'\\]"), (["name", "count", "extra"], "Unknown columns: \\['extra'\\]")],
)
def test_apply_schema_rejects_a_missing_or_unknown_column(columns, message):
    table = pd.DataFrame({column: [1] for column in columns})

    with pytest.raises(ValueError, match=message):
        apply_schema(table, {"name": "string", "count": "int"})


def test_apply_schema_rejects_a_fraction_in_an_int_column():
    with pytest.raises((TypeError, ValueError)):
        apply_schema(pd.DataFrame({"count": [1.5]}), {"count": "int"})


def test_empty_table_has_the_schema():
    table = empty_table()

    assert len(table) == 0
    assert_schema(table, OUTPUT_COLUMN_KINDS)


@pytest.fixture
def run_main(tmp_path, monkeypatch):
    """
    Return a function that runs main() on a VCF with the chr18 annotation. It returns the table that
    extract_ptc returned to main() and the table that main() wrote.
    """

    def run(vcf_path):
        tables = {}

        def recording_extract_ptc(*args, **kwargs):
            tables["extract_ptc"] = extract_ptc(*args, **kwargs)
            return tables["extract_ptc"]

        def recording_write_results(results, output):
            tables["written"] = results

        monkeypatch.setattr(cli, "extract_ptc", recording_extract_ptc)
        monkeypatch.setattr(cli, "write_results", recording_write_results)
        cli.main(vcf_path, "resources/chr18.gtf.gz", "resources/chr18.fa.gz", str(tmp_path / "out.parquet"))
        return tables["extract_ptc"], tables["written"]

    return run


@pytest.mark.parametrize(
    "vcf, has_rows",
    [
        ("resources/test_files/variants.vcf", True),
        ("intergenic_vcf", False),
        ("reference_mismatch_vcf", False),
    ],
)
def test_empty_and_nonempty_results_have_the_same_schema(run_main, request, vcf, has_rows):
    vcf_path = request.getfixturevalue(vcf) if vcf.endswith("_vcf") else vcf

    ptc_table, written = run_main(vcf_path)

    assert (len(ptc_table) > 0) == has_rows
    assert (len(written) > 0) == has_rows
    assert_schema(ptc_table, PTC_COLUMN_KINDS)
    assert_schema(written, OUTPUT_COLUMN_KINDS)


def test_feature_and_rule_columns_match_the_schema(run_main):
    _, written = run_main("resources/test_files/variants.vcf")
    row = written.iloc[0]

    assert list(add_nmd_features(row)) == list(NMD_FEATURE_COLUMN_KINDS)
    assert list(evaluate_nmd_escape_rules(row)) == list(NMD_RULE_COLUMN_KINDS)


def schema_row(**values):
    """Return one row of the result table with ``values`` set and every other column missing."""

    table = pd.DataFrame([values]).reindex(columns=list(OUTPUT_COLUMN_KINDS))
    return apply_schema(table, OUTPUT_COLUMN_KINDS).iloc[0]


def none_row(**values):
    """Return ``values`` as a dict row with None for every other output column."""

    return {**dict.fromkeys(OUTPUT_COLUMN_KINDS), **values}


@pytest.mark.parametrize(
    "values",
    [
        {},
        {"alt_is_premature": True, "alt_first_stop_pos": 30, "cds_in_transcript": True, "ref_valid_stop": True},
        {
            "alt_is_premature": True,
            "alt_start_codon_pos": 0,
            "alt_first_stop_pos": 30,
            "alt_stop_codon_exons": [1],
            "alt_cds_info": [(1, 60), (2, 90)],
            "transcript_exon_info": [(1, 100), (2, 120)],
            "cds_in_transcript": True,
            "ref_start_codon_pos": 0,
            "ref_valid_stop": False,
        },
    ],
)
def test_features_and_rules_treat_pd_na_like_none(values):
    features = add_nmd_features(schema_row(**values))

    assert features == add_nmd_features(none_row(**values))
    assert evaluate_nmd_escape_rules(schema_row(**values, **features)) == evaluate_nmd_escape_rules(
        none_row(**values, **features)
    )


def test_likely_misannotated_reads_numpy_bools():
    # .iloc gives numpy.bool_ values, and `numpy.False_ is False` is False
    row = schema_row(cds_in_transcript=False, ref_start_codon_pos=0, ref_valid_stop=True)
    assert isinstance(row["cds_in_transcript"], np.bool_)

    assert add_nmd_features(row)["likely_misannotated"] is True
