from dataclasses import dataclass

from nmd_scanner.ptc import (
    annotated_stop_distance,
    classify_rescued_orf,
    ends_at_annotated_stop,
    first_stop_codon,
    starts_with_annotated_start_codon,
)
from nmd_scanner.schema import record


@dataclass(frozen=True, slots=True)
class FeatureRow:
    """A variant-transcript pair with the columns of extract_ptc that add_nmd_features reads, see schema.records."""

    unknown_reason: str | None
    has_start_codon: bool
    has_stop_codon: bool
    cds_frame: int
    ref_cds_seq: str
    alt_cds_seq: str | None
    alt_cds_length: int | None
    cds_in_transcript: bool
    ref_valid_stop: bool | None
    alt_first_stop_pos: int | None
    alt_has_ptc: bool | None
    start_loss: bool | None
    alt_scan_start_codon_pos: int | None
    alt_scan_first_stop_pos: int | None
    transcript_seq: str | None
    cds_start_in_transcript: int | None
    cds_end_in_transcript: int | None
    transcript_exons: list[dict] | None
    alt_transcript_seq: str | None
    alt_cds_start_in_transcript: int | None
    alt_transcript_exons: list[dict] | None


def ptc_pos_in_alt_transcript(row):
    """
    Return the position of the PTC in alt_transcript_seq: the position of its first base.

    The PTC is the first in-frame stop codon of the alt CDS, at alt_first_stop_pos in alt CDS coordinates. Its position
    in the alt transcript adds alt_cds_start_in_transcript. After a start loss, translation starts at the ATG of the
    scan, and the PTC is the first in-frame stop codon after it, at alt_scan_first_stop_pos in the alt transcript.

    :param row: A FeatureRow or an nmd_rules.RuleRow
    :return: The position, or None if the variant-transcript pair is not a PTC variant or has no alt_transcript_seq
    """

    if not row.alt_has_ptc:
        return None
    if row.start_loss:
        return row.alt_scan_first_stop_pos
    stop = row.alt_first_stop_pos
    cds_start = row.alt_cds_start_in_transcript
    return None if stop is None or cds_start is None else cds_start + stop


def ptc_in_alt_transcript(row):
    """
    Locate the PTC in the exons of the alt transcript, the mRNA.

    The PTC position comes from ptc_pos_in_alt_transcript. alt_transcript_exons gives the exons of the alt transcript,
    with the length changes of the variant in the CDS and in the UTRs. An exon that the variant deletes has length 0.
    It is not in the mRNA, so it is left out.

    :param row: A FeatureRow or an nmd_rules.RuleRow
    :return: Tuple (position, exons, index): the position of the PTC in the alt transcript, (exon_number, start, end)
             of each exon of the alt transcript with bases, 5' to 3', in alt transcript positions (end is the position
             after the last base), and the index of the PTC exon in this list. None if the variant-transcript pair is
             not a PTC variant, or if the PTC position or alt_transcript_exons is null.
    """

    position = ptc_pos_in_alt_transcript(row)
    exon_info = row.alt_transcript_exons
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

    :param row: A variant-transcript pair of extract_ptc: a FeatureRow, or a pandas Series or a dict with a value for
                each field of FeatureRow
    :return: A dictionary with additional NMD related features.
    """

    if not isinstance(row, FeatureRow):
        row = record(row, FeatureRow)

    # Without an alt transcript (unknown_reason), only the features of the reference are known
    if row.unknown_reason is not None:
        return {
            **calculate_utr_lengths(row),
            "total_exon_count": calculate_exon_features(row)["total_exon_count"],
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "upstream_exon_count": None,
            "downstream_exon_count": None,
            "ptc_to_start_codon": None,
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

    # The PTC in alt_transcript_seq, and the exon number of its exon
    ptc_position = ptc_pos_in_alt_transcript(row)
    ptc_exon_number = calculate_ptc_exon_number(row)

    # Distance between PTC to start codon
    ptc_to_start_codon = calculate_ptc_to_start_distance(row)

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
        "ptc_pos_in_alt_transcript": ptc_position,
        "ptc_exon_number": ptc_exon_number,
        "upstream_exon_count": upstream_exon_count,
        "downstream_exon_count": downstream_exon_count,
        "ptc_to_start_codon": ptc_to_start_codon,
        "ptc_exon_length": ptc_exon_length,
        "annotated_stop_distance": stop_distance,
        "ptc_to_exon_end": ptc_to_exon_end,
        "likely_misannotated": likely_misannotated,
    }


def calculate_utr_lengths(row):
    """
    Calculate the 5' and 3' UTR lengths of the reference transcript from the position of the CDS in it.
    The 3'UTR starts after the stop codon, so its length is None without an annotated stop codon (has_stop_codon False).

    :param row: A FeatureRow: it has cds_start_in_transcript and cds_end_in_transcript (coding region, from
                cds_range_in_transcript), transcript_exons and has_stop_codon
    :return: A dictionary with utr5_length and utr3_length, both None if the CDS position or the exons are unknown
    """

    cds_start = row.cds_start_in_transcript
    cds_end = row.cds_end_in_transcript
    transcript_exons = row.transcript_exons or []
    has_stop_codon = bool(row.has_stop_codon)

    if cds_start is None or cds_end is None or not transcript_exons:
        return {"utr5_length": None, "utr3_length": None}

    utr3 = sum(int(exon["length"]) for exon in transcript_exons) - cds_end
    return {"utr5_length": cds_start, "utr3_length": utr3 if utr3 >= 0 and has_stop_codon else None}


def calculate_exon_features(row):
    """
    Calculate exon-related features:
    - total_exon_count: the number of exons of the ref transcript (transcript_exons)
    - upstream_exon_count / downstream_exon_count: the number of exons of the alt transcript upstream and downstream of
      the PTC exon, only for a PTC variant (see ptc_in_alt_transcript). An exon that the variant deletes is not in the
      mRNA, so it does not count there.
    """

    total_exons = len(row.transcript_exons or [])
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
    alt CDS must start with the annotated start codon (see ptc.starts_with_annotated_start_codon). The distance is None if the transcript has no annotated start codon (e.g.
    cds_start_NF: the true start lies upstream of the CDS, at an unknown distance).
    After a start loss, translation starts at the ATG that the scan of the alt transcript found. So the distance runs
    from that ATG (alt_scan_start_codon_pos) to the first in-frame stop codon after it (alt_scan_first_stop_pos),
    both alt transcript positions. Such a variant-transcript pair is a PTC variant only if the scan found both (see
    ptc.classify_rescued_orf).
    The distance is None, too, if the PTC is the annotated start codon itself, i.e. the annotated start codon is a stop
    codon such as TAG. Translation cannot start on a stop codon, so such a start codon is a misannotation.
    """

    if not row.alt_has_ptc:
        return None

    if row.start_loss:
        start = row.alt_scan_start_codon_pos
        stop = row.alt_scan_first_stop_pos
    else:
        start = 0 if starts_with_annotated_start_codon(row, "alt") else None
        stop = row.alt_first_stop_pos

    if start is None or stop is None:
        return None

    # Without a start loss, the start codon lies at CDS position 0, so only a stop codon as start codon gives a PTC at
    # the start. After a start loss, the scan reads the stop codons after the ATG.
    if stop <= start:
        return None

    return stop - start  # distance between the PTC to start codon in nt


def calculate_ptc_exon_number(row):
    """
    Return the exon number of the PTC exon, from alt_transcript_exons (see ptc_in_alt_transcript). None if there is no
    PTC exon.
    """

    ptc = ptc_in_alt_transcript(row)
    if ptc is None:
        return None
    _, exons, index = ptc
    return exons[index][0]


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
    A variant-transcript pair that keeps the flags from the CDS there takes alt_first_stop_pos, the first stop codon
    of the alt CDS.
    Without an alt transcript or alt_cds_start_in_transcript, both positions are in alt CDS coordinates, and the reference
    stop codon is the last codon of the alt coding region, at alt_cds_length - 3. That holds only for a variant upstream of the stop codon: an
    insertion inside it (TAA>TGAA) lengthens the alt coding region but leaves the stop codon in place. An indel
    upstream of the PTC shifts both positions by the same amount, so the distance is the same as in ref CDS
    coordinates.
    Without an annotated stop codon (has_stop_codon False), there is no reference stop codon and the distance is None.
    After a start loss, the alternative stop codon is the first in-frame stop codon after the ATG of the scan (see
    ptc.classify_rescued_orf). The distance is None if the scan found no ATG, or an ATG downstream of the reference
    stop codon.
    """

    if not row.has_stop_codon:
        return None

    alt_seq = row.alt_transcript_seq
    alt_stop = row.alt_first_stop_pos
    alt_cds_start = row.alt_cds_start_in_transcript
    if isinstance(alt_seq, str) and alt_cds_start is not None and row.start_loss:
        return classify_rescued_orf(row, row.alt_scan_start_codon_pos, row.alt_scan_first_stop_pos)[2]
    if not isinstance(alt_seq, str) or alt_cds_start is None:
        alt_cds_length = row.alt_cds_length
        if alt_cds_length is None or alt_stop is None:
            return None
        return alt_cds_length - 3 - alt_stop

    if ends_at_annotated_stop(row):
        first_stop = first_stop_codon(alt_seq, alt_cds_start + int(row.cds_frame))
    else:
        first_stop = None if alt_stop is None else alt_cds_start + alt_stop
    return annotated_stop_distance(row, first_stop)


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
    Flag variant-transcript pairs that look inconsistent between CDS and transcript annotations and might be likely
    misannotated. A variant-transcript pair is flagged as likely misannotated if any of these conditions apply:
        cds_in_transcript = False (the assembled CDS is not found in the transcript sequence)
        the reference CDS does not start with an annotated start codon, e.g. cds_start_NF (see
        ptc.starts_with_annotated_start_codon)
        ref_valid_stop is False (the reference does not end in a valid annotated stop codon, e.g. for a transcript
        without stop_codon records such as one tagged cds_end_NF)

    :return: A boolean flag. True if any condition above is met and thus the variant-transcript pair is likely
             misannotated, False otherwise.
    """

    # Add likely_misannotated flag: when
    # "cds_in_transcript" is FALSE
    # the ref CDS does not start with an annotated start codon
    # "ref_valid_stop" is FALSE

    cds_in_transcript = row.cds_in_transcript
    ref_valid_stop = row.ref_valid_stop

    # if any of these are missing entirely, flag as likely misannotated
    if cds_in_transcript is None or ref_valid_stop is None:
        return True

    flag = (
        (cds_in_transcript is False) or not starts_with_annotated_start_codon(row, "ref") or (ref_valid_stop is False)
    )

    return flag
