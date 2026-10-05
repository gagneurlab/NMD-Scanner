import numpy as np
import pandas as pd
import pytest

from nmd_scanner import cli
from nmd_scanner.extra_features import (
    add_features_and_rules,
    add_nmd_features,
    evaluate_nmd_escape_rules,
    nmd_model_status,
)
from nmd_scanner.rules import extract_ptc
from nmd_scanner.schema import (
    KIND_DTYPES,
    MODEL_INPUTS,
    MODEL_STATUS_COLUMN_KINDS,
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


def test_output_columns_are_the_ptc_feature_rule_and_status_columns_in_order():
    parts = [*PTC_COLUMN_KINDS, *NMD_FEATURE_COLUMN_KINDS, *NMD_RULE_COLUMN_KINDS, *MODEL_STATUS_COLUMN_KINDS]
    assert list(OUTPUT_COLUMN_KINDS) == parts
    assert len(set(parts)) == len(parts)


def test_model_inputs_are_19_distinct_output_columns():
    assert len(MODEL_INPUTS) == 19
    assert len(set(MODEL_INPUTS)) == 19
    assert set(MODEL_INPUTS) <= set(NMD_FEATURE_COLUMN_KINDS) | set(NMD_RULE_COLUMN_KINDS) | set(PTC_COLUMN_KINDS)


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
    extract_ptc returned to main() and the table that main() wrote. After each call, the attribute
    ptc_table_before_main of the function holds a copy of the first table, taken before main() used it.
    """

    def run(vcf_path):
        tables = {}

        def recording_extract_ptc(*args, **kwargs):
            tables["extract_ptc"] = extract_ptc(*args, **kwargs)
            run.ptc_table_before_main = tables["extract_ptc"].copy()
            return tables["extract_ptc"]

        def recording_write_results(results, output):
            tables["written"] = results

        monkeypatch.setattr(cli, "extract_ptc", recording_extract_ptc)
        monkeypatch.setattr(cli, "write_results", recording_write_results)
        cli.main(vcf_path, "resources/chr18.gff3.gz", "resources/chr18.fa.gz", str(tmp_path / "out.parquet"))
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


@pytest.mark.parametrize("vcf", ["resources/test_files/variants.vcf", "intergenic_vcf", "reference_mismatch_vcf"])
def test_add_features_and_rules_returns_the_schema_for_zero_and_many_rows(run_main, request, vcf):
    vcf_path = request.getfixturevalue(vcf) if vcf.endswith("_vcf") else vcf
    ptc_table, _ = run_main(vcf_path)

    result = add_features_and_rules(ptc_table)

    assert len(result) == len(ptc_table)
    assert_schema(result, OUTPUT_COLUMN_KINDS)


def test_add_features_and_rules_gives_the_same_dtypes_for_zero_and_many_rows(run_main, intergenic_vcf):
    many = add_features_and_rules(run_main("resources/test_files/variants.vcf")[0])
    zero = add_features_and_rules(run_main(intergenic_vcf)[0])

    assert len(many) > 0
    assert len(zero) == 0
    assert list(zero.columns) == list(many.columns)
    assert dict(zero.dtypes) == dict(many.dtypes)


def test_add_features_and_rules_equals_what_main_writes(run_main):
    ptc_table, written = run_main("resources/test_files/test_variants.vcf")
    # main() already passed ptc_table to add_features_and_rules, so a copy taken now could miss a change
    unchanged = run_main.ptc_table_before_main

    result = add_features_and_rules(ptc_table)

    pd.testing.assert_frame_equal(result, written)
    pd.testing.assert_frame_equal(ptc_table, unchanged)


def test_add_features_and_rules_equals_the_row_functions_applied_one_by_one(run_main):
    ptc_table, _ = run_main("resources/test_files/test_variants.vcf")

    features = ptc_table.apply(add_nmd_features, axis=1, result_type="expand")
    table = pd.concat([ptc_table, features], axis=1)
    rules = table.apply(evaluate_nmd_escape_rules, axis=1, result_type="expand")
    table = pd.concat([table, rules], axis=1)
    expected = apply_schema(table.assign(nmd_model_status=nmd_model_status(table)), OUTPUT_COLUMN_KINDS)

    pd.testing.assert_frame_equal(add_features_and_rules(ptc_table), expected)


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
            "alt_cds_start_in_transcript": 40,
            "transcript_exon_info": [(1, 100), (2, 120)],
            "alt_transcript_exon_info": [(1, 100), (2, 120)],
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
