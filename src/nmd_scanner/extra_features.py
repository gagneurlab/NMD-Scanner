import math

import numpy as np
import pandas as pd

from nmd_scanner.rules import annotated_stop_distance, classify_rescued_orf, ends_at_annotated_stop, first_stop_codon
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


def ptc_in_alt_transcript(row):
    """
    Locate the PTC in the exons of the alt transcript, the mRNA.

    The PTC is the first in-frame stop codon of the alt CDS, at alt_first_stop_pos in alt CDS coordinates. Its position
    in the alt transcript adds alt_cds_start_in_transcript. After a start loss, translation starts at the ATG of the
    scan, and the PTC is the first in-frame stop codon after it, at transcript_first_stop_pos in the alt transcript.
    alt_transcript_exon_info gives the exons of the alt transcript, with the length changes of the variant in the CDS
    and in the UTRs. An exon that the variant deletes has length 0. It is not in the mRNA, so it is left out.

    :param row: A dict of plain values (see _plain_values)
    :return: Tuple (position, exons, index): the position of the PTC in the alt transcript, (exon_number, start, end)
             of each exon of the alt transcript with bases, 5' to 3', in alt transcript positions (end is the position
             after the last base), and the index of the PTC exon in this list. None if the row is not a PTC row, or if
             the PTC position or alt_transcript_exon_info is null.
    """

    if not row.get("alt_is_premature"):
        return None
    if row.get("start_loss"):
        position = row.get("transcript_first_stop_pos")
    else:
        stop = row.get("alt_first_stop_pos")
        cds_start = row.get("alt_cds_start_in_transcript")
        position = None if stop is None or cds_start is None else cds_start + stop
    exon_info = row.get("alt_transcript_exon_info")
    if position is None or not exon_info:
        return None

    exons = []
    end = 0
    for exon_number, length in exon_info:
        if int(length) > 0:
            exons.append((int(exon_number), end, end + int(length)))
        end += int(length)
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
            "stop_codon_distance": None,
            "ptc_to_intron": None,
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
    # PTC location < 150nt to start codon
    ptc_less_than_150nt_to_start = ptc_to_start_codon is not None and ptc_to_start_codon < 150

    # PTC exon length
    ptc_exon_length = calculate_ptc_exon_length(row)

    # Distance PTC to normal stop codon
    stop_codon_distance = calculate_stop_codon_dist(row)

    # Distance PTC to the 3' end of the PTC exon: the downstream exon junction, or the transcript end
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
    - total_exon_count: the number of exons of the ref transcript (transcript_exon_info)
    - upstream_exon_count / downstream_exon_count: the number of exons of the alt transcript upstream and downstream of
      the PTC exon, only on a PTC row (see ptc_in_alt_transcript). An exon that the variant deletes is not in the mRNA,
      so it does not count there.
    """

    total_exons = len(row.get("transcript_exon_info") or [])
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
    The start codon is the annotated one at CDS position 0 (alt_start_codon_pos), and the PTC lies at
    alt_first_stop_pos, both CDS positions. The distance is None if the transcript has no annotated start codon (e.g.
    cds_start_NF: the true start lies upstream of the CDS, at an unknown distance).
    After a start loss, translation starts at the ATG that the scan of the alt transcript found. So the distance runs
    from that ATG (transcript_start_codon_pos) to the first in-frame stop codon after it (transcript_first_stop_pos),
    both alt transcript positions. Such a row is a PTC row only if the scan found both (see rules.classify_rescued_orf).
    """

    if not row.get("alt_is_premature"):
        return None

    if row.get("start_loss"):
        start = row.get("transcript_start_codon_pos")
        stop = row.get("transcript_first_stop_pos")
    else:
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
    Return the length of the PTC exon in the ref transcript (transcript_exon_info), UTR included. The PTC exon is the
    exon of the alt transcript that holds the PTC (see ptc_in_alt_transcript). None if there is none.
    """

    ptc = ptc_in_alt_transcript(row)
    if ptc is None:
        return None
    _, exons, index = ptc
    ref_lengths = {int(exon_number): int(length) for exon_number, length in row.get("transcript_exon_info") or []}
    return ref_lengths.get(exons[index][0])


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
    stop codon is the last codon of the alt coding region, at alt_cds_len - 3. That holds only for a variant upstream of the stop codon: an
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
        return classify_rescued_orf(row, row.get("transcript_start_codon_pos"), row.get("transcript_first_stop_pos"))[2]
    if not isinstance(alt_seq, str) or alt_cds_start is None:
        alt_cds_len = row.get("alt_cds_len")
        if alt_cds_len is None or alt_stop is None:
            return None
        return alt_cds_len - 3 - alt_stop

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

    :param row: A row of the DataFrame including alt_is_premature (bool), the columns that ptc_in_alt_transcript reads
                (the PTC position and alt_transcript_exon_info), total_exon_count, downstream_exon_count and
                ptc_exon_length (see add_nmd_features), and the columns that calculate_ptc_to_start_distance reads
    :return: A dictionary with boolean flags for each rule and overall NMD escape
    """

    row = _plain_values(row)

    # Unknown without an alt transcript
    if row.get("unknown_reason") is not None:
        return {
            "nmd_last_exon_rule": None,
            "nmd_50nt_penultimate_rule": None,
            "nmd_long_exon_rule": None,
            "nmd_start_proximal_rule": None,
            "nmd_single_exon_rule": None,
            "nmd_escape": None,
        }

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
    ptc = ptc_in_alt_transcript(row)

    total_exons = row.get("total_exon_count")
    downstream_exons = row.get("downstream_exon_count")
    ptc_exon_length = row.get("ptc_exon_length")

    # Single exon rule
    rule_single_exon = total_exons == 1

    # Last exon rule
    rule_last_exon = downstream_exons == 0 if downstream_exons is not None else False

    # 1 to 50 nt upstream of the last exon junction of the alt transcript, i.e. the 3' end of its penultimate exon.
    # The junction lies past the CDS end if the last exon holds no CDS.
    rule_50nt_penultimate = False
    if ptc is not None and len(ptc[1]) >= 2:
        stop_pos, exons, _ = ptc
        last_junction = exons[-2][2]
        rule_50nt_penultimate = last_junction - 50 <= stop_pos < last_junction

    # Long exon rule (with exon longer than >407nt)
    rule_long_exon = ptc_exon_length is not None and ptc_exon_length > 407

    # Start-proximal rule (closer than 150nt from the start codon). Without a known start codon, the rule does not
    # apply. After a start loss, the start codon is the ATG of the scan.
    ptc_to_start_codon = calculate_ptc_to_start_distance(row)
    rule_start_proximal = ptc_to_start_codon is not None and ptc_to_start_codon < 150

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
        ref_start_codon_pos is None or not 0 (the reference CDS does not start with an annotated start codon, e.g.
        cds_start_NF)
        ref_valid_stop is False (the reference does not end in a valid annotated stop codon, e.g. for a transcript
        without stop_codon rows such as one tagged cds_end_NF)

    :return: A boolean flag. True if any condition above is met and thus the row is likely misannotated, False otherwise.
    """

    # Add likely_misannotated flag: when
    # "cds_in_transcript" is FALSE
    # "ref_start_codon_pos" is not 0 (None: no annotated start codon)
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
