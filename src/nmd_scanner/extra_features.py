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

    return {
        "3UTR_length": utr3_length,
        "5UTR_length": utr5_length,
        "total_exon_count": total_exon_count,
        "upstream_exon_count": upstream_exon_count,
        "downstream_exon_count": downstream_exon_count,
        # "ptc_pos_codon": ptc_pos_codon,
        "ptc_to_start_codon": ptc_to_start_codon,
        "ptc_less_than_150nt_to_start": ptc_less_than_150nt_to_start,
        "ptc_exon_length": ptc_exon_length,
        "stop_codon_distance": stop_codon_distance
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
            utr5 = max(0, row["ref_cds_start"] - row["transcript_start"])
            utr3 = max(0, row["transcript_end"] - row["ref_cds_stop"])
        else:
            utr5 = max(0, row["transcript_end"] - row["ref_cds_stop"])
            utr3 = max(0, row["ref_cds_start"] - row["transcript_start"])
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

    return {"utr5_length": utr5, "utr3_length": utr3}

def calculate_exon_features(row):

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
        return None  # TODO: How to handle if PTC before start codon?

    return stop-start # distance between the PTC to start codon in nt

def calculate_ptc_exon_length(row):

    if not row.get("alt_is_premature"):
        return None

    stop_exons = row.get("alt_stop_codon_exons") or []
    exon_info = row.get("transcript_exon_info") or []

    if not stop_exons or not exon_info:
        return None

    ptc_exon = int(stop_exons[0])
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