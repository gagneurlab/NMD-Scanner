import math

import numpy as np
import pandas as pd

from nmd_scanner.schema import OUTPUT_COLUMN_KINDS, apply_schema


def _plain_values(row):
    """
    Return the values of ``row`` as a dict of plain Python values, with None for a missing value.

    A row of the result table holds pd.NA for a missing value in an int, bool or string column. A row of a
    table without that schema can hold NaN instead, and a row taken with ``.iloc`` holds numpy scalars. The
    feature and rule functions test for a missing value with ``is None`` and for False with ``is False``, so
    they need plain values. Lists, e.g. ``transcript_exon_info``, stay as they are.

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


def add_nmd_features(row):
    """
    Compute additional features which might be relevant for analyzing nonsense-mediated decay (NMD) behavior,
    inspired by benchmark datasets from nmd_eff. These features include UTR lengths, exon structure and positional information
    of the premature termination codon (PTC).

    :param row: A DataFrame row with annotated transcript information
    :return: A dictionary with additional NMD related features.
    """

    row = _plain_values(row)

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
    # PTC location < 150nt to start codon
    ptc_less_than_150nt_to_start = ptc_to_start_codon is not None and ptc_to_start_codon < 150

    # PTC exon length
    ptc_exon_length = calculate_ptc_exon_length(row)

    # Distance PTC to normal stop codon
    stop_codon_distance = calculate_stop_codon_dist(row)

    # Distance PTC to downstream exon junction
    ptc_to_intron = calculate_ptc_to_downstream_ej(row)

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
        "stop_codon_distance": stop_codon_distance,
        "ptc_to_intron": ptc_to_intron,
        "likely_misannotated": likely_misannotated,
    }


def calculate_utr_lengths(row):
    """
    Calculate the 5' and 3' UTR lengths of the reference transcript from the position of the CDS in it.
    The 3'UTR starts after the stop codon, so its length is None without an annotated stop codon (has_stop_codon False).

    :param row: A row of the DataFrame including cds_start_in_transcript and cds_end_in_transcript
                (coding region, from cds_range_in_transcript), transcript_exon_info and has_stop_codon
    :return: A dictionary with utr5_length and utr3_length, both None if the CDS position or the exons are unknown
    """

    cds_start = row.get("cds_start_in_transcript")
    cds_end = row.get("cds_end_in_transcript")
    transcript_exon_info = row.get("transcript_exon_info") or []
    has_stop_codon = bool(row["has_stop_codon"])

    if cds_start is None or cds_end is None or not transcript_exon_info:
        return {"utr5_length": None, "utr3_length": None}

    utr3 = sum(int(length) for _, length in transcript_exon_info) - cds_end
    return {"utr5_length": cds_start, "utr3_length": utr3 if utr3 >= 0 and has_stop_codon else None}


def calculate_exon_features(row):
    """
    Calculate exon-related features:
    - total_exon_count: always computed if transcript_exon_info is available
    - upstream_exon_count / downstream_exon_count: only computed if a PTC exists
    """

    exon_info = row.get("transcript_exon_info") or []
    stop_exons = row.get("alt_stop_codon_exons") or []

    total_exons = len(exon_info)

    if not row.get("alt_is_premature") or not stop_exons or not exon_info:
        return {
            "total_exon_count": total_exons if total_exons > 0 else None,
            "upstream_exon_count": None,
            "downstream_exon_count": None,
        }

    # Take the PTC exon closest to CDS start
    ptc_exon = min(int(e) for e in stop_exons)

    # get exon numbers from transcript_exon_info
    exon_numbers = [int(e[0]) for e in exon_info]

    # If the PTC exon is not in transcript → cannot compute
    if ptc_exon not in exon_numbers:
        return {"total_exon_count": int(total_exons), "upstream_exon_count": None, "downstream_exon_count": None}

    upstream = sum(1 for e in exon_numbers if e < ptc_exon)
    downstream = sum(1 for e in exon_numbers if e > ptc_exon)

    return {
        "total_exon_count": int(total_exons),
        "upstream_exon_count": int(upstream),
        "downstream_exon_count": int(downstream),
    }


def calculate_ptc_to_start_distance(row):

    if not row.get("alt_is_premature"):
        return None

    start = row.get("alt_start_codon_pos")
    stop = row.get("alt_first_stop_pos")

    if start is None or stop is None:
        return None

    # PTC codon position: is PTC_to_start_codon / 3 --> leave it out
    # offset = stop - start
    # return offset // 3 if offset >= 0 else None

    if stop <= start:
        return None

    return stop - start  # distance between the PTC to start codon in nt


def calculate_ptc_exon_length(row):
    """
    Return the length of the exon containing the first premature stop codon (PTC).
    """

    if not row.get("alt_is_premature"):
        return None

    stop_exons = row.get("alt_stop_codon_exons") or []
    exon_info = row.get("transcript_exon_info") or []

    if not stop_exons or not exon_info:
        return None

    # First PTC exon = smallest exon number (transcript-order, strand-corrected)
    ptc_exon = min(int(e) for e in stop_exons)

    exon_dict = {int(e): int(length) for e, length in exon_info}
    return exon_dict.get(ptc_exon)


def calculate_stop_codon_dist(row):
    """
    Calculate the distance in nt between the reference stop codon and the alternative stop codon (alt_first_stop_pos).
    Positive means the PTC is upstream of the reference stop codon, 0 means the alternative stop codon is the
    reference stop codon.

    Both positions are in alt CDS coordinates. The reference stop codon is the annotated one: the last codon of the
    alt coding region, at alt_cds_len - 3. The first in-frame stop of the reference is not used, because it can be an
    internal one, e.g. a selenocysteine TGA. An indel upstream of the PTC shifts both positions by the same amount,
    so the distance is the same as in ref CDS coordinates.
    Without an annotated stop codon (has_stop_codon False), there is no reference stop codon and the distance is None.
    """

    if not row["has_stop_codon"]:
        return None

    alt_cds_len = row.get("alt_cds_len")
    alt_stop = row.get("alt_first_stop_pos")

    if alt_cds_len is None or alt_stop is None:
        return None

    return alt_cds_len - 3 - alt_stop


def exon_end_in_alt_cds(row, exon):
    """
    Return where a transcript exon ends in alt CDS coordinates (as alt_first_stop_pos), or None if the CDS position
    in the transcript is unknown. The value is the CDS position of the first base after the exon, i.e. of its
    downstream exon junction. It is negative for an exon upstream of the CDS.

    Exon numbers follow transcript order. The exon end in transcript coordinates, minus cds_start_in_transcript,
    gives the position in ref CDS coordinates. The length change of the CDS up to this exon converts it to alt CDS
    coordinates.
    """

    cds_start = row.get("cds_start_in_transcript")
    tx_exons = row.get("transcript_exon_info") or []
    alt_cds = {int(e): int(length) for e, length in row.get("alt_cds_info") or []}
    ref_cds = {int(e): int(length) for e, length in row.get("ref_cds_info") or []}

    if cds_start is None or not tx_exons or alt_cds.keys() != ref_cds.keys():
        return None

    exon_end = sum(int(length) for e, length in tx_exons if int(e) <= exon)
    cds_length_change = sum(alt_cds[e] - ref_cds[e] for e in alt_cds if e <= exon)
    return exon_end - cds_start + cds_length_change


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
    A PTC is considered to escape NMD if it satisfies any of the above rules.

    :param row: A row of the DataFrame including alt_is_premature (bool), alt_first_stop_pos (int),
                alt_stop_codon_exons (list[int]), transcript_exon_info (list[tuple[exon_number (int), exon_length (int)]]),
                alt_cds_info and ref_cds_info (same format, CDS part per exon), cds_start_in_transcript (int),
                alt_start_codon_pos (int)
    :return: A dictionary with boolean flags for each rule and overall NMD escape
    """

    row = _plain_values(row)

    # Only relevant for premature stop codons
    if not row.get("alt_is_premature"):
        return {
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": False,
        }

    # Extract relevant data
    stop_pos = row.get("alt_first_stop_pos")
    start_pos = row.get("alt_start_codon_pos")
    tx_exon_nums = sorted(int(e) for e, _ in row.get("transcript_exon_info") or [])

    total_exons = row.get("total_exon_count")
    downstream_exons = row.get("downstream_exon_count")
    ptc_exon_length = row.get("ptc_exon_length")

    # Single exon rule
    rule_single_exon = total_exons == 1

    # Last exon rule
    rule_last_exon = downstream_exons == 0 if downstream_exons is not None else False

    # 50nt upstream of the last exon junction, i.e. the 3' end of the penultimate exon (CDS-relative, matching stop_pos).
    # The junction lies past the CDS end if the last exon holds no CDS.
    last_junction = exon_end_in_alt_cds(row, tx_exon_nums[-2]) if len(tx_exon_nums) >= 2 else None
    rule_50nt_penultimate = (
        last_junction is not None
        and stop_pos is not None
        and (stop_pos >= last_junction - 50)
        and (stop_pos < last_junction)
    )

    # Long exon rule (with exon longer than >407nt)
    # rule_long_exon = any(exon_length_map.get(exon, 0) > 407 for exon in stop_exons) # old code
    rule_long_exon = ptc_exon_length is not None and ptc_exon_length > 407

    # Start-proximal rule (closer than 150nt from the start codon)
    rule_start_proximal = (
        start_pos is not None and stop_pos is not None and (stop_pos - start_pos) < 150 and (stop_pos - start_pos) >= 0
    )

    # NMD escape if any rule is true
    escape = rule_last_exon or rule_50nt_penultimate or rule_long_exon or rule_start_proximal or rule_single_exon

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
    Calculate distance from PTC to the downstream exon junction, i.e. the 3' end of the PTC exon.
    Returns None if not applicable. This includes a PTC in the last exon, which has no downstream junction.
    """

    # only calculate if we have PTC
    if not row.get("alt_is_premature"):
        return None

    stop_exons = row.get("alt_stop_codon_exons") or []
    tx_exon_nums = [int(e) for e, _ in row.get("transcript_exon_info") or []]
    ptc_pos = row.get("alt_first_stop_pos")

    if not stop_exons or not tx_exon_nums or ptc_pos is None:
        return None

    # Choose the PTC exon (smallest number, closer to start)
    ptc_exon = min(stop_exons)

    if ptc_exon >= max(tx_exon_nums):
        return None

    # The PTC exon can go on with 3' UTR, so its junction can lie past the CDS end
    exon_end = exon_end_in_alt_cds(row, ptc_exon)
    return exon_end - ptc_pos if exon_end is not None else None


def add_likely_misannotated_flag(row):
    """
    Flag rows that look inconsistent between CDS and transcript annotations and might be likely misannotated.
    A row is flagged as likely misannotated if any of these conditions apply:
        cds_in_transcript = False (the assembled CDS is not found in the transcript sequence)
        ref_start_codon_pos is defined and not 0 (reference CDS has a start codon not at the very start)
        ref_valid_stop is False (the reference does not end in a valid annotated stop codon, e.g. for a transcript
        without stop_codon rows such as one tagged cds_end_NF)

    :return: A boolean flag. True if any condition above is met and thus the row is likely misannotated, False otherwise.
    """

    # Add likely_misannotated flag: when
    # "cds_in_transcript" is FALSE
    # "ref_start_codon_pos" is not 0
    # "ref_valid_stop" is FALSE

    cds_in_transcript = row.get("cds_in_transcript")
    ref_start_codon_pos = row.get("ref_start_codon_pos")
    ref_valid_stop = row.get("ref_valid_stop")

    # if any of these are missing entirely, flag as likely misannotated
    if cds_in_transcript is None or ref_start_codon_pos is None or ref_valid_stop is None:
        return True

    flag = (
        (cds_in_transcript is False)
        or ((ref_start_codon_pos is not None) and (ref_start_codon_pos != 0))
        or (ref_valid_stop is False)
    )

    return flag


def add_features_and_rules(results):
    """
    Add the NMD features and the NMD escape rules to a result of ``extract_ptc``.

    This runs ``add_nmd_features`` on each row, then ``evaluate_nmd_escape_rules`` (it reads columns that the
    features add), and applies the output schema. The result has the columns, column order and dtypes of
    OUTPUT_COLUMN_KINDS (see nmd_scanner.schema), also for zero rows. ``results`` is not changed.

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
    return apply_schema(results, OUTPUT_COLUMN_KINDS)
