import math

import numpy as np
import pandas as pd

from nmd_scanner.rules import (
    annotated_stop_distance,
    classify_rescued_orf,
    ends_at_annotated_stop,
    first_stop_codon,
    starts_with_annotated_start_codon,
)
from nmd_scanner.schema import (
    KIND_DTYPES,
    MODEL_INPUTS,
    MODEL_STATUSES,
    NMD_RULE_COLUMN_KINDS,
    OUTPUT_COLUMN_KINDS,
    apply_schema,
)


def _plain_values(row):
    """
    Return the values of ``row`` as a dict of plain Python values, with None for a missing value.

    A row of the result table holds pd.NA for a missing value in an int, bool or string column. A row of a
    table without that schema can hold NaN instead, and a row taken with ``.iloc`` holds numpy scalars. The
    feature and rule functions test for a missing value with ``is None`` and for False with ``is False``, so
    they need plain values. Lists, e.g. ``transcript_exons``, stay as they are.

    :param row: A DataFrame row (pandas Series) or a dict
    :return: A dict of column name to value
    """

    values = {}
    for column, value in row.items():
        if isinstance(value, np.generic):
            value = value.item()
        if value is pd.NA or (isinstance(value, float) and math.isnan(value)):
            value = None
        values[column] = value
    return values


def ptc_in_alt_transcript(row):
    """
    Locate the PTC in the exons of the alt transcript, the mRNA.

    The PTC is the first in-frame stop codon of the alt CDS, at alt_first_stop_pos in alt CDS coordinates. Its position
    in the alt transcript adds alt_cds_start_in_transcript. After a start loss, translation starts at the ATG of the
    scan, and the PTC is the first in-frame stop codon after it, at alt_scan_first_stop_pos in the alt transcript.
    alt_transcript_exons gives the exons of the alt transcript, with the length changes of the variant in the CDS
    and in the UTRs. An exon that the variant deletes has length 0. It is not in the mRNA, so it is left out.

    :param row: A dict of plain values (see _plain_values)
    :return: Tuple (position, exons, index): the position of the PTC in the alt transcript, (exon_number, start, end)
             of each exon of the alt transcript with bases, 5' to 3', in alt transcript positions (end is the position
             after the last base), and the index of the PTC exon in this list. None if the row is not a PTC row, or if
             the PTC position or alt_transcript_exons is null.
    """

    if not row.get("alt_has_ptc"):
        return None
    if row.get("start_loss"):
        position = row.get("alt_scan_first_stop_pos")
    else:
        stop = row.get("alt_first_stop_pos")
        cds_start = row.get("alt_cds_start_in_transcript")
        position = None if stop is None or cds_start is None else cds_start + stop
    exon_info = row.get("alt_transcript_exons")
    if position is None or not exon_info:
        return None

    exons = []
    end = 0
    for exon in exon_info:
        length = int(exon["length"])
        if length > 0:
            exons.append((int(exon["exon_number"]), end, end + length))
        end += length
    index = next((i for i, (_, exon_start, exon_end) in enumerate(exons) if exon_start <= position < exon_end), None)
    return None if index is None else (position, exons, index)


def add_nmd_features(row):
    """
    Compute additional features which might be relevant for analyzing nonsense-mediated decay (NMD) behavior,
    inspired by benchmark datasets from nmd_eff. These features include UTR lengths, exon structure and positional information
    of the premature termination codon (PTC).

    :param row: A DataFrame row with annotated transcript information
    :return: A dictionary with additional NMD related features.
    """

    row = _plain_values(row)

    # Without an alt transcript (unknown_reason), only the features of the reference are known
    if row.get("unknown_reason") is not None:
        return {
            **calculate_utr_lengths(row),
            "total_exon_count": calculate_exon_features(row)["total_exon_count"],
            "upstream_exon_count": None,
            "downstream_exon_count": None,
            "ptc_to_start_codon": None,
            "ptc_less_than_150nt_to_start": None,
            "ptc_exon_length": None,
            "annotated_stop_distance": None,
            "ptc_to_exon_end": None,
            "likely_misannotated": add_likely_misannotated_flag(row),
        }

    # 5' and 3' UTR lengths
    utr_lengths = calculate_utr_lengths(row)
    utr3_length = utr_lengths["utr3_length"]
    utr5_length = utr_lengths["utr5_length"]

    # Total, Upstream and Downstream exon count
    exon_features = calculate_exon_features(row)
    total_exon_count = exon_features["total_exon_count"]
    upstream_exon_count = exon_features["upstream_exon_count"]
    downstream_exon_count = exon_features["downstream_exon_count"]

    # Distance between PTC to start codon
    ptc_to_start_codon = calculate_ptc_to_start_distance(row)
    # PTC location < 150nt to start codon. Null without a distance, e.g. on a row that is not a PTC row.
    ptc_less_than_150nt_to_start = None if ptc_to_start_codon is None else ptc_to_start_codon < 150

    # PTC exon length
    ptc_exon_length = calculate_ptc_exon_length(row)

    # Distance PTC to normal stop codon
    stop_distance = calculate_stop_codon_dist(row)

    # Distance PTC to the 3' end of the PTC exon: the downstream exon junction, or the transcript end
    ptc_to_exon_end = calculate_ptc_to_downstream_ej(row)

    # Add likely_misannotated flag
    likely_misannotated = add_likely_misannotated_flag(row)

    return {
        "utr3_length": utr3_length,
        "utr5_length": utr5_length,
        "total_exon_count": total_exon_count,
        "upstream_exon_count": upstream_exon_count,
        "downstream_exon_count": downstream_exon_count,
        # "ptc_pos_codon": ptc_pos_codon,
        "ptc_to_start_codon": ptc_to_start_codon,
        "ptc_less_than_150nt_to_start": ptc_less_than_150nt_to_start,
        "ptc_exon_length": ptc_exon_length,
        "annotated_stop_distance": stop_distance,
        "ptc_to_exon_end": ptc_to_exon_end,
        "likely_misannotated": likely_misannotated,
    }


def calculate_utr_lengths(row):
    """
    Calculate the 5' and 3' UTR lengths of the reference transcript from the position of the CDS in it.
    The 3'UTR starts after the stop codon, so its length is None without an annotated stop codon (has_stop_codon False).

    :param row: A row of the DataFrame including cds_start_in_transcript and cds_end_in_transcript
                (coding region, from cds_range_in_transcript), transcript_exons and has_stop_codon
    :return: A dictionary with utr5_length and utr3_length, both None if the CDS position or the exons are unknown
    """

    cds_start = row.get("cds_start_in_transcript")
    cds_end = row.get("cds_end_in_transcript")
    transcript_exons = row.get("transcript_exons") or []
    has_stop_codon = bool(row["has_stop_codon"])

    if cds_start is None or cds_end is None or not transcript_exons:
        return {"utr5_length": None, "utr3_length": None}

    utr3 = sum(int(exon["length"]) for exon in transcript_exons) - cds_end
    return {"utr5_length": cds_start, "utr3_length": utr3 if utr3 >= 0 and has_stop_codon else None}


def calculate_exon_features(row):
    """
    Calculate exon-related features:
    - total_exon_count: the number of exons of the ref transcript (transcript_exons)
    - upstream_exon_count / downstream_exon_count: the number of exons of the alt transcript upstream and downstream of
      the PTC exon, only on a PTC row (see ptc_in_alt_transcript). An exon that the variant deletes is not in the mRNA,
      so it does not count there.
    """

    total_exons = len(row.get("transcript_exons") or [])
    ptc = ptc_in_alt_transcript(row)
    if ptc is None:
        return {
            "total_exon_count": total_exons if total_exons > 0 else None,
            "upstream_exon_count": None,
            "downstream_exon_count": None,
        }

    _, exons, index = ptc
    return {
        "total_exon_count": total_exons if total_exons > 0 else None,
        "upstream_exon_count": index,
        "downstream_exon_count": len(exons) - 1 - index,
    }


def calculate_ptc_to_start_distance(row):
    """
    Calculate the distance in nt from the start codon to the PTC.
    The start codon is the annotated one at CDS position 0, and the PTC lies at alt_first_stop_pos, a CDS position. The
    alt CDS must start with the annotated start codon (see rules.starts_with_annotated_start_codon). The distance is None if the transcript has no annotated start codon (e.g.
    cds_start_NF: the true start lies upstream of the CDS, at an unknown distance).
    After a start loss, translation starts at the ATG that the scan of the alt transcript found. So the distance runs
    from that ATG (alt_scan_start_codon_pos) to the first in-frame stop codon after it (alt_scan_first_stop_pos),
    both alt transcript positions. Such a row is a PTC row only if the scan found both (see rules.classify_rescued_orf).
    The distance is None, too, if the PTC is the annotated start codon itself, i.e. the annotated start codon is a stop
    codon such as TAG. Translation cannot start on a stop codon, so such a start codon is a misannotation.
    """

    if not row.get("alt_has_ptc"):
        return None

    if row.get("start_loss"):
        start = row.get("alt_scan_start_codon_pos")
        stop = row.get("alt_scan_first_stop_pos")
    else:
        start = 0 if starts_with_annotated_start_codon(row, "alt") else None
        stop = row.get("alt_first_stop_pos")

    if start is None or stop is None:
        return None

    # PTC codon position: is PTC_to_start_codon / 3 --> leave it out
    # offset = stop - start
    # return offset // 3 if offset >= 0 else None

    # Without a start loss, the start codon lies at CDS position 0, so only a stop codon as start codon gives a PTC at
    # the start. After a start loss, the scan reads the stop codons after the ATG.
    if stop <= start:
        return None

    return stop - start  # distance between the PTC to start codon in nt


def calculate_ptc_exon_length(row):
    """
    Return the length of the PTC exon in the alt transcript, as in the mRNA, UTR included (see ptc_in_alt_transcript).
    An indel in the PTC exon changes it. None if there is no PTC exon.
    """

    ptc = ptc_in_alt_transcript(row)
    if ptc is None:
        return None
    _, exons, index = ptc
    _, start, end = exons[index]
    return end - start


def calculate_stop_codon_dist(row):
    """
    Calculate the distance in nt between the reference stop codon and the alternative stop codon.
    Positive means the PTC is upstream of the reference stop codon, 0 means the alternative stop codon is the
    reference stop codon, and negative means it lies downstream (stop loss). Without an alternative stop codon, e.g.
    for a nonstop, the distance is None.

    The reference stop codon is the annotated one. The first in-frame stop of the reference is not used, because it
    can be an internal one, e.g. a selenocysteine TGA. With an alt transcript, the distance is the one by which
    analyze_transcript classifies the first in-frame stop codon of the alt transcript (see annotated_stop_distance).
    A row that keeps the flags from the CDS there takes alt_first_stop_pos, the first stop codon of the alt CDS.
    Without an alt transcript or alt_cds_start_in_transcript, both positions are in alt CDS coordinates, and the reference
    stop codon is the last codon of the alt coding region, at alt_cds_length - 3. That holds only for a variant upstream of the stop codon: an
    insertion inside it (TAA>TGAA) lengthens the alt coding region but leaves the stop codon in place. An indel
    upstream of the PTC shifts both positions by the same amount, so the distance is the same as in ref CDS
    coordinates.
    Without an annotated stop codon (has_stop_codon False), there is no reference stop codon and the distance is None.
    After a start loss, the alternative stop codon is the first in-frame stop codon after the ATG of the scan (see
    rules.classify_rescued_orf). The distance is None if the scan found no ATG, or an ATG downstream of the reference
    stop codon.
    """

    if not row["has_stop_codon"]:
        return None

    alt_seq = row.get("alt_transcript_seq")
    alt_stop = row.get("alt_first_stop_pos")
    alt_cds_start = row.get("alt_cds_start_in_transcript")
    if isinstance(alt_seq, str) and alt_cds_start is not None and row.get("start_loss"):
        return classify_rescued_orf(row, row.get("alt_scan_start_codon_pos"), row.get("alt_scan_first_stop_pos"))[2]
    if not isinstance(alt_seq, str) or alt_cds_start is None:
        alt_cds_length = row.get("alt_cds_length")
        if alt_cds_length is None or alt_stop is None:
            return None
        return alt_cds_length - 3 - alt_stop

    if ends_at_annotated_stop(row):
        first_stop = first_stop_codon(alt_seq, alt_cds_start + int(row["cds_frame"]))
    else:
        first_stop = None if alt_stop is None else alt_cds_start + alt_stop
    return annotated_stop_distance(row, first_stop)


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

    A rule is None on a row that is not a PTC row, also on one with unknown_reason, because it does not apply there.
    On a PTC row, a rule is None if one of its inputs is None, because its value is unknown then. nmd_escape is
    True if one rule is True, None if no rule is True and one is None, and False otherwise.

    :param row: A row of the DataFrame including alt_has_ptc (bool), the columns that ptc_in_alt_transcript reads
                (the PTC position and alt_transcript_exons), total_exon_count, downstream_exon_count and
                ptc_exon_length (see add_nmd_features), and the columns that calculate_ptc_to_start_distance reads
    :return: A dictionary with a flag or None for each rule and for the overall NMD escape
    """

    row = _plain_values(row)

    # Only relevant for premature stop codons. A row with unknown_reason has no known PTC either.
    if not row.get("alt_has_ptc"):
        return dict.fromkeys(NMD_RULE_COLUMN_KINDS)

    # Extract relevant data
    ptc = ptc_in_alt_transcript(row)

    total_exons = row.get("total_exon_count")
    downstream_exons = row.get("downstream_exon_count")
    ptc_exon_length = row.get("ptc_exon_length")

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


def calculate_ptc_to_downstream_ej(row):
    """
    Calculate the distance from the PTC to the 3' end of the PTC exon in the alt transcript (see ptc_in_alt_transcript).
    For an internal exon, that end is the downstream exon junction. For the last exon, it is the transcript end, so the
    distance is the length of the 3' UTR that the PTC creates. Returns None if not applicable. "Technical Notes.md" has a
    figure of each case.
    """

    ptc = ptc_in_alt_transcript(row)
    if ptc is None:
        return None
    position, exons, index = ptc
    return exons[index][2] - position


def add_likely_misannotated_flag(row):
    """
    Flag rows that look inconsistent between CDS and transcript annotations and might be likely misannotated.
    A row is flagged as likely misannotated if any of these conditions apply:
        cds_in_transcript = False (the assembled CDS is not found in the transcript sequence)
        the reference CDS does not start with an annotated start codon, e.g. cds_start_NF (see
        rules.starts_with_annotated_start_codon)
        ref_valid_stop is False (the reference does not end in a valid annotated stop codon, e.g. for a transcript
        without stop_codon rows such as one tagged cds_end_NF)

    :return: A boolean flag. True if any condition above is met and thus the row is likely misannotated, False otherwise.
    """

    # Add likely_misannotated flag: when
    # "cds_in_transcript" is FALSE
    # the ref CDS does not start with an annotated start codon
    # "ref_valid_stop" is FALSE

    cds_in_transcript = row.get("cds_in_transcript")
    ref_valid_stop = row.get("ref_valid_stop")

    # if any of these are missing entirely, flag as likely misannotated
    if cds_in_transcript is None or ref_valid_stop is None:
        return True

    flag = (
        (cds_in_transcript is False) or not starts_with_annotated_start_codon(row, "ref") or (ref_valid_stop is False)
    )

    return flag


def add_features_and_rules(results):
    """
    Add the NMD features and the NMD escape rules to a result of ``extract_ptc``.

    This runs ``add_nmd_features`` on each row, then ``evaluate_nmd_escape_rules`` (it reads columns that the
    features add). Then it adds the column nmd_model_status (see ``nmd_model_status``) and applies the output schema.
    The result has the columns, column order and dtypes of OUTPUT_COLUMN_KINDS (see nmd_scanner.schema), also for zero
    rows. ``results`` is not changed.

    :param results: DataFrame returned by ``extract_ptc``, also an empty one
    :return: DataFrame with the columns and dtypes of OUTPUT_COLUMN_KINDS
    """

    if results.empty:
        # DataFrame.apply with result_type="expand" returns the input columns again for zero rows
        return apply_schema(results.reindex(columns=list(OUTPUT_COLUMN_KINDS)), OUTPUT_COLUMN_KINDS)

    extra_features = results.apply(add_nmd_features, axis=1, result_type="expand")
    results = pd.concat([results, extra_features], axis=1)
    nmd_results = results.apply(evaluate_nmd_escape_rules, axis=1, result_type="expand")
    results = pd.concat([results, nmd_results], axis=1)
    results["nmd_model_status"] = nmd_model_status(results)
    return apply_schema(results, OUTPUT_COLUMN_KINDS)


def nmd_model_status(results):
    """
    Return the nmd_model_status of each row: "ok" if the NMD efficiency model can score the row, else the first
    reason why it cannot. MODEL_STATUSES in nmd_scanner.schema lists the values and their conditions, in the order
    they are checked.

    :param results: DataFrame with unknown_reason, alt_has_ptc, ref_has_ptc, has_stop_codon, has_start_codon
                    and the columns of MODEL_INPUTS
    :return: Series of the values of MODEL_STATUSES, with the index of ``results`` and the string dtype
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
        is_true("start_loss") & results["ptc_to_start_codon"].isna().to_numpy(),
        results[MODEL_INPUTS].isna().any(axis=1).to_numpy(),
    ]
    status = np.select(conditions, list(MODEL_STATUSES[:-1]), default=MODEL_STATUSES[-1])
    return pd.Series(status, index=results.index, dtype=KIND_DTYPES["string"])
