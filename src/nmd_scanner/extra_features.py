"""
Compute additional features which might be relevant for analyzing nonsense-mediated decay (NMD) behavior,
inspired by benchmark datasets from nmd_eff. These features include UTR lengths, exon structure and positional information
of the premature termination codon (PTC).

:param row: A DataFrame row with annotated transcript information
:return: A dictionary with additional NMD related features.
"""

def add_nmd_features(row):

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
    ptc_less_than_150nt_to_start = (
            ptc_to_start_codon is not None and ptc_to_start_codon < 150
    )

    # PTC exon length
    ptc_exon_length = calculate_ptc_exon_length(row)

    # Distance PTC to normal stop codon
    stop_codon_distance = calculate_stop_codon_dist(row)

    # Distance PTC to downstream exon junction
    ptc_to_intron = calculate_ptc_to_downstream_ej(row)

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
    }

def calculate_utr_lengths(row):
    strand = row.get("strand")
    ref_cds_info = row.get("ref_cds_info") or []
    transcript_exon_info = row.get("transcript_exon_info") or []

    if not ref_cds_info or not transcript_exon_info:
        return {"utr5_length": None, "utr3_length": None}

    # Convert to dicts for easier lookup
    transcript_exon_dict = {int(k): int(v) for k, v in transcript_exon_info}
    cds_exons_dict = {int(exon): int(length) for exon, length in ref_cds_info}

    # Handle single exon
    if len(transcript_exon_dict) == 1:
        if strand == "+":
            utr5 = row["ref_cds_start"] - row["transcript_start"]
            utr3 = row["transcript_end"] - row["ref_cds_stop"]
        else:
            utr5 = row["transcript_end"] - row["ref_cds_stop"]
            utr3 = row["ref_cds_start"] - row["transcript_start"]

        utr5 = utr5 if utr5 >= 0 else None
        utr3 = utr3 if utr3 >= 0 else None

        return {"utr5_length": utr5, "utr3_length": utr3}

    cds_exon_nums = sorted(cds_exons_dict.keys())
    tx_exon_nums = sorted(transcript_exon_dict.keys())

    utr5 = 0
    utr3 = 0

    for exon in tx_exon_nums:
        exon_len = transcript_exon_dict[exon]
        cds_len = cds_exons_dict.get(exon, 0)
        utr_len = exon_len - cds_len

        if utr_len <= 0:
            continue

        if exon in cds_exons_dict:
            # Exon overlaps CDS, partial UTR
            if strand == "+":
                if exon == cds_exon_nums[0]:
                    utr5 += utr_len
                elif exon == cds_exon_nums[-1]:
                    utr3 += utr_len
            else:
                if exon == cds_exon_nums[0]:
                    utr3 += utr_len
                elif exon == cds_exon_nums[-1]:
                    utr5 += utr_len
        else:
            # Exon is outside CDS
            if strand == "+":
                if exon < cds_exon_nums[0]:
                    utr5 += exon_len
                elif exon > cds_exon_nums[-1]:
                    utr3 += exon_len
            else:
                if exon > cds_exon_nums[-1]:
                    utr5 += exon_len
                elif exon < cds_exon_nums[0]:
                    utr3 += exon_len

    utr5 = utr5 if utr5 >= 0 else None
    utr3 = utr3 if utr3 >= 0 else None
    return {"utr5_length": utr5, "utr3_length": utr3}

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
            "downstream_exon_count": None
        }

    # Take the PTC exon closest to CDS start
    ptc_exon = min(int(e) for e in stop_exons)

    # get exon numbers from transcript_exon_info
    exon_numbers = [int(e[0]) for e in exon_info]

    # If the PTC exon is not in transcript → cannot compute
    if ptc_exon not in exon_numbers:
        return {
            "total_exon_count": int(total_exons),
            "upstream_exon_count": None,
            "downstream_exon_count": None
        }

    upstream = sum(1 for e in exon_numbers if e < ptc_exon)
    downstream = sum(1 for e in exon_numbers if e > ptc_exon)

    return {
        "total_exon_count": int(total_exons),
        "upstream_exon_count": int(upstream),
        "downstream_exon_count": int(downstream)
    }

def calculate_ptc_to_start_distance(row):

    if not row.get("alt_is_premature"):
        return None

    start = row.get("alt_start_codon_pos")
    stop = row.get("alt_first_stop_pos")

    if start is None or stop is None:
        return None

    # PTC codon position: is PTC_to_start_codon / 3 --> leave it out
    #offset = stop - start
    #return offset // 3 if offset >= 0 else None

    if stop <= start:
        return None

    return stop-start # distance between the PTC to start codon in nt

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
    Calculate the distance between the reference stop codon and the alternative stop codon.
    Positive means the PTC is upstream of the reference stop codon.
    """

    ref_stop = row.get("ref_first_stop_pos")
    alt_stop = row.get("alt_first_stop_pos")

    if ref_stop is None or alt_stop is None:
        return None

    return ref_stop - alt_stop

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
                alt_start_codon_pos (int)
    :return: A dictionary with boolean flags for each rule and overall NMD escape
    """

    # Only relevant for premature stop codons
    if not row.get("alt_is_premature"):
        return {
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": False
        }

    # Extract relevant data
    stop_pos = row.get("alt_first_stop_pos")
    stop_exons = row.get("alt_stop_codon_exons") or []
    exon_info = row.get("transcript_exon_info") or []
    start_pos = row.get("alt_start_codon_pos")

    total_exons = row.get("total_exon_count")
    downstream_exons = row.get("downstream_exon_count")
    ptc_exon_length = row.get("ptc_exon_length")

    # Preprocess exon info
    sorted_exons = sorted(exon_info, key=lambda x: x[0])
    exon_length_map = {exon_num: length for exon_num, length in sorted_exons}
    exon_offsets = {}
    offset = 0
    for exon_num, length in sorted_exons:
        exon_offsets[exon_num] = (offset, offset + length)
        offset += length

    # Single exon rule
    # rule_single_exon = len(sorted_exons) == 1 # old code
    rule_single_exon = total_exons == 1

    # Last exon rule
    #if not rule_single_exon:
    #    last_exon = sorted_exons[-1][0] if sorted_exons else None
    #    rule_last_exon = bool(stop_exons and max(stop_exons) == last_exon)
    #else: rule_last_exon = False
    rule_last_exon = downstream_exons == 0 if downstream_exons is not None else False

    # 50nt from penultimate exon end
    if len(sorted_exons) >= 2:
        penultimate_exon_num, penultimate_len = sorted_exons[-2]
        pen_start, pen_end = exon_offsets.get(penultimate_exon_num, (None, None))
        rule_50nt_penultimate = pen_end is not None and (stop_pos >= pen_end - 50) and (stop_pos < pen_end)
    else:
        rule_50nt_penultimate = False

    # Long exon rule (with exon longer than >407nt)
    # rule_long_exon = any(exon_length_map.get(exon, 0) > 407 for exon in stop_exons) # old code
    rule_long_exon = ptc_exon_length is not None and ptc_exon_length > 407

    # Start-proximal rule (closer than 150nt from the start codon)
    rule_start_proximal = start_pos is not None and stop_pos is not None and (stop_pos - start_pos) < 150 and (stop_pos - start_pos) >= 0

    # NMD escape if any rule is true
    escape = rule_last_exon or rule_50nt_penultimate or rule_long_exon or rule_start_proximal or rule_single_exon

    return {
        "nmd_last_exon_rule": rule_last_exon,
        "nmd_50nt_penultimate_rule": rule_50nt_penultimate,
        "nmd_long_exon_rule": rule_long_exon,
        "nmd_start_proximal_rule": rule_start_proximal,
        "nmd_single_exon_rule": rule_single_exon,
        "nmd_escape": escape
    }

def calculate_ptc_to_downstream_ej(row):

    """
    Calculate distance from PTC to the downstream exon junction (next exon start/end depending on strand).
    Returns None if not applicable.
    """

    # only calculate if we have PTC
    if not row.get("alt_is_premature"):
        return None

    stop_exons = row.get("alt_stop_codon_exons") or []
    exon_info = row.get("transcript_exon_info") or []
    ptc_pos = row.get("alt_first_stop_pos")

    if not stop_exons or not exon_info or ptc_pos is None:
        return None

    # Choose the PTC exon (smallest number, closer to start)
    ptc_exon = min(stop_exons)

    # Sum lengths of exons up to and including ptc_exon
    cumulative_length = 0
    for exon_num, length in exon_info:
        cumulative_length += length
        if exon_num == ptc_exon:
            break

    # Distance from PTC to downstream exon junction
    distance = cumulative_length - ptc_pos
    return distance


def calculate_exon_features_old(row):

    """
    Calculate exon-related features:
    - total_exon_count: always computed if transcript_exon_info is available
    - upstream_exon_count / downstream_exon_count: only computed if a PTC exists
    """

    exon_info = row.get("transcript_exon_info") or []
    stop_exons = row.get("alt_stop_codon_exons") or []
    strand = row.get("strand")

    # Total exon count
    total_exons = len(exon_info)


    if not row.get("alt_is_premature") or not stop_exons:
        return {
            "total_exon_count": total_exons if total_exons > 0 else None,
            "upstream_exon_count": None,
            "downstream_exon_count": None
        }

    if not exon_info:
        return {
            "total_exon_count": None,
            "upstream_exon_count": None,
            "downstream_exon_count": None
        }

    # Order exons based on strand
    sorted_exons = sorted(exon_info, key=lambda x: int(x[0]), reverse=(strand == "-"))
    exon_order = [int(e[0]) for e in sorted_exons]

    ptc_exon = int(stop_exons[0])  # Take first exon containing PTC
    if ptc_exon not in exon_order:
        return {
            "total_exon_count": total_exons,
            "upstream_exon_count": None,
            "downstream_exon_count": None
        }

    ptc_index = exon_order.index(ptc_exon)
    upstream = ptc_index
    downstream = len(exon_order) - ptc_index - 1

    return {
        "total_exon_count": total_exons,
        "upstream_exon_count": upstream,
        "downstream_exon_count": downstream
    }