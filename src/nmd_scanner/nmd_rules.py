from dataclasses import dataclass

import numpy as np
import pandas as pd

from nmd_scanner.nmd_features import (
    FeatureRow,
    add_nmd_features,
    calculate_ptc_to_start_distance,
    ptc_in_alt_transcript,
)
from nmd_scanner.schema import (
    KIND_DTYPES,
    MODEL_INPUTS,
    NMD_RULE_COLUMN_KINDS,
    OUTPUT_COLUMN_KINDS,
    PREDICTION_STATUSES,
    apply_schema,
    record,
    records,
)


@dataclass(frozen=True, slots=True)
class RuleRow:
    """
    A variant-transcript pair with the columns of extract_ptc that evaluate_nmd_escape_rules reads, and 3 columns that
    add_nmd_features adds, see schema.records.
    """

    has_start_codon: bool
    ref_cds_seq: str
    alt_cds_seq: str | None
    alt_first_stop_pos: int | None
    alt_has_ptc: bool | None
    start_loss: bool | None
    alt_scan_start_codon_pos: int | None
    alt_scan_first_stop_pos: int | None
    alt_cds_start_in_transcript: int | None
    alt_transcript_exons: list[dict] | None
    total_exon_count: int | None
    downstream_exon_count: int | None
    ptc_exon_length: int | None


def evaluate_nmd_escape_rules(row):
    """
    Evaluate whether a premature stop codon in a transcript is likely to escape nonsense-mediated decay (NMD) based on
    established biological rules. This function applies five NMD escape rules to determine if a premature termination
    codon (PTC) is likely to escape degradation:
    1. Last exon rule: The PTC is in the last exon
    2. 50nt penultimate rule: The PTC is within 50 nucleotides upstream of the last exon junction
    3. Long exon rule: The PTC is in an exon with >407 nucleotides
    4. Start proximal rule: The PTC is within 150 nucleotides of the start codon
    5. Single exon rule: The transcript where the PTC lays consists only of a single exon
    A PTC is considered to escape NMD if it satisfies any of the above rules. "Technical Notes.md" has figures of the
    rules.

    A rule is None for a variant-transcript pair that is not a PTC variant, also for one with unknown_reason, because
    it does not apply there. For a PTC variant, a rule is None if one of its inputs is None, because its value is
    unknown then. nmd_escape is True if one rule is True, None if no rule is True and one is None, and False otherwise.

    :param row: A variant-transcript pair of extract_ptc with the features of add_nmd_features: a RuleRow, or a
                pandas Series or a dict with a value for each field of RuleRow. These are alt_has_ptc (bool), the
                columns that ptc_in_alt_transcript reads (the PTC position and alt_transcript_exons),
                total_exon_count, downstream_exon_count and ptc_exon_length (see add_nmd_features), and the columns
                that calculate_ptc_to_start_distance reads.
    :return: A dictionary with a flag or None for each rule and for the overall NMD escape
    """

    if not isinstance(row, RuleRow):
        row = record(row, RuleRow)

    # Only relevant for premature stop codons. A variant-transcript pair with unknown_reason has no known PTC either.
    if not row.alt_has_ptc:
        return dict.fromkeys(NMD_RULE_COLUMN_KINDS)

    # Extract relevant data
    ptc = ptc_in_alt_transcript(row)

    total_exons = row.total_exon_count
    downstream_exons = row.downstream_exon_count
    ptc_exon_length = row.ptc_exon_length

    # Single exon rule
    rule_single_exon = None if total_exons is None else total_exons == 1

    # Last exon rule
    rule_last_exon = None if downstream_exons is None else downstream_exons == 0

    # 1 to 50 nt upstream of the last exon junction of the alt transcript, i.e. the 3' end of its penultimate exon.
    # The junction lies past the CDS end if the last exon holds no CDS. A transcript of one exon has no junction.
    rule_50nt_penultimate = None if ptc is None else False
    if ptc is not None and len(ptc[1]) >= 2:
        stop_pos, exons, _ = ptc
        last_junction = exons[-2][2]
        rule_50nt_penultimate = last_junction - 50 <= stop_pos < last_junction

    # Long exon rule (with exon longer than >407nt)
    rule_long_exon = None if ptc_exon_length is None else ptc_exon_length > 407

    # Start-proximal rule (closer than 150nt from the start codon). Without a known start codon, the distance and the
    # rule are unknown. After a start loss, the start codon is the ATG of the scan.
    ptc_to_start_codon = calculate_ptc_to_start_distance(row)
    rule_start_proximal = None if ptc_to_start_codon is None else ptc_to_start_codon < 150

    # NMD escape if any rule is true, unknown if no rule is true and one is unknown (three-valued OR)
    rules = [rule_last_exon, rule_50nt_penultimate, rule_long_exon, rule_start_proximal, rule_single_exon]
    escape = True if True in rules else None if None in rules else False

    return {
        "nmd_last_exon_rule": rule_last_exon,
        "nmd_50nt_penultimate_rule": rule_50nt_penultimate,
        "nmd_long_exon_rule": rule_long_exon,
        "nmd_start_proximal_rule": rule_start_proximal,
        "nmd_single_exon_rule": rule_single_exon,
        "nmd_escape": escape,
    }


def add_features_and_rules(results):
    """
    Add the NMD features and the NMD escape rules to a result of ``extract_ptc``.

    This runs ``add_nmd_features`` on the FeatureRow of each variant-transcript pair, then ``evaluate_nmd_escape_rules``
    on the RuleRow of each pair (it reads columns that the features add). Then it adds the column
    nmd_efficiency_prediction_status (see ``nmd_efficiency_prediction_status``) and applies the output schema.
    The result has the columns, column order and dtypes of OUTPUT_COLUMN_KINDS (see nmd_scanner.schema), also for zero
    rows. ``results`` is not changed.

    :param results: DataFrame returned by ``extract_ptc``, also an empty one
    :return: DataFrame with the columns and dtypes of OUTPUT_COLUMN_KINDS
    """

    if results.empty:
        # Without rows, the lists of features and rules below give DataFrames without columns
        return apply_schema(results.reindex(columns=list(OUTPUT_COLUMN_KINDS)), OUTPUT_COLUMN_KINDS)

    features = pd.DataFrame([add_nmd_features(row) for row in records(results, FeatureRow)], index=results.index)
    results = pd.concat([results, features], axis=1)
    rules = pd.DataFrame([evaluate_nmd_escape_rules(row) for row in records(results, RuleRow)], index=results.index)
    results = pd.concat([results, rules], axis=1)
    results["nmd_efficiency_prediction_status"] = nmd_efficiency_prediction_status(results)
    return apply_schema(results, OUTPUT_COLUMN_KINDS)


def nmd_efficiency_prediction_status(results):
    """
    Return the nmd_efficiency_prediction_status of each variant-transcript pair: "ok" if the NMD efficiency model can
    predict the pair, else the first reason why it cannot. PREDICTION_STATUSES in nmd_scanner.schema lists the values
    and their conditions, in the order they are checked.

    :param results: DataFrame with unknown_reason, alt_has_ptc, ref_has_ptc, has_stop_codon, has_start_codon
                    and the columns of MODEL_INPUTS
    :return: Series of the values of PREDICTION_STATUSES, with the index of ``results`` and the categorical dtype of
             nmd_efficiency_prediction_status
    """

    def is_true(column):
        return results[column].astype("boolean").fillna(False).to_numpy(dtype=bool)

    def is_false(column):
        return ~results[column].astype("boolean").fillna(True).to_numpy(dtype=bool)

    conditions = [
        results["unknown_reason"].notna().to_numpy(),
        ~is_true("alt_has_ptc"),
        is_true("ref_has_ptc"),
        is_false("has_stop_codon"),
        is_false("has_start_codon"),
        is_true("start_loss"),
        results[MODEL_INPUTS].isna().any(axis=1).to_numpy(),
    ]
    status = np.select(conditions, list(PREDICTION_STATUSES[:-1]), default=PREDICTION_STATUSES[-1])
    return pd.Series(status, index=results.index, dtype=KIND_DTYPES["nmd_efficiency_prediction_status"])
