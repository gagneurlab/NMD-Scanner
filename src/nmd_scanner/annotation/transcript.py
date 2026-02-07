"""
Transcript-level feature calculation for NMD analysis
"""

from .models import CDSAnnotation, TranscriptFeatures


def calculate_transcript_features(cds: CDSAnnotation) -> TranscriptFeatures:
    """
    Calculate transcript-level features for NMD analysis.
    
    :param cds: CDSAnnotation object with CDS and transcript information
    :return: TranscriptFeatures object with computed features
    """
    
    # 5' and 3' UTR lengths
    utr_lengths = _calculate_utr_lengths(
        strand=cds.strand,
        ref_cds_info=cds.ref_cds_info,
        transcript_exon_info=cds.transcript_exon_info,
        ref_cds_start=cds.ref_cds_start,
        transcript_start=cds.transcript_start,
        ref_cds_stop=cds.ref_cds_stop,
        transcript_end=cds.transcript_end
    )
    
    # Total, Upstream and Downstream exon count
    exon_features = _calculate_exon_features(
        transcript_exon_info=cds.transcript_exon_info,
        alt_stop_codon_exons=cds.alt_stop_codon_exons,
        alt_is_premature=cds.alt_is_premature
    )
    
    # Distance between PTC to start codon
    ptc_to_start_codon = _calculate_ptc_to_start_distance(
        alt_is_premature=cds.alt_is_premature,
        alt_start_codon_pos=cds.alt_start_codon_pos,
        alt_first_stop_pos=cds.alt_first_stop_pos
    )
    
    # PTC location < 150nt to start codon
    ptc_less_than_150nt_to_start = (
        ptc_to_start_codon is not None and ptc_to_start_codon < 150
    )
    
    # PTC exon length
    ptc_exon_length = _calculate_ptc_exon_length(
        alt_is_premature=cds.alt_is_premature,
        alt_stop_codon_exons=cds.alt_stop_codon_exons,
        transcript_exon_info=cds.transcript_exon_info
    )
    
    # Distance PTC to normal stop codon
    stop_codon_distance = _calculate_stop_codon_dist(
        ref_first_stop_pos=cds.ref_first_stop_pos,
        alt_first_stop_pos=cds.alt_first_stop_pos
    )
    
    # Distance PTC to downstream exon junction
    ptc_to_intron = _calculate_ptc_to_downstream_ej(
        alt_is_premature=cds.alt_is_premature,
        alt_stop_codon_exons=cds.alt_stop_codon_exons,
        alt_cds_info=cds.alt_cds_info,
        alt_first_stop_pos=cds.alt_first_stop_pos
    )
    
    # Add likely_misannotated flag
    likely_misannotated = _add_likely_misannotated_flag(
        cds_in_transcript=cds.cds_in_transcript,
        ref_start_codon_pos=cds.ref_start_codon_pos,
        ref_valid_stop=cds.ref_valid_stop
    )
    
    return TranscriptFeatures(
        utr3_length=utr_lengths["utr3_length"],
        utr5_length=utr_lengths["utr5_length"],
        total_exon_count=exon_features["total_exon_count"],
        upstream_exon_count=exon_features["upstream_exon_count"],
        downstream_exon_count=exon_features["downstream_exon_count"],
        ptc_to_start_codon=ptc_to_start_codon,
        ptc_less_than_150nt_to_start=ptc_less_than_150nt_to_start,
        ptc_exon_length=ptc_exon_length,
        stop_codon_distance=stop_codon_distance,
        ptc_to_intron=ptc_to_intron,
        likely_misannotated=likely_misannotated
    )


def _calculate_utr_lengths(strand, ref_cds_info, transcript_exon_info, ref_cds_start, 
                           transcript_start, ref_cds_stop, transcript_end):
    """Calculate 5' and 3' UTR lengths"""
    ref_cds_info = ref_cds_info or []
    transcript_exon_info = transcript_exon_info or []
    
    if not ref_cds_info or not transcript_exon_info:
        return {"utr5_length": None, "utr3_length": None}
    
    # Convert to dicts for easier lookup
    transcript_exon_dict = {int(k): int(v) for k, v in transcript_exon_info}
    cds_exons_dict = {int(exon): int(length) for exon, length in ref_cds_info}
    
    # Handle single exon
    if len(transcript_exon_dict) == 1:
        if strand == "+":
            utr5 = ref_cds_start - transcript_start
            utr3 = transcript_end - ref_cds_stop
        else:
            utr5 = transcript_end - ref_cds_stop
            utr3 = ref_cds_start - transcript_start
        
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


def _calculate_exon_features(transcript_exon_info, alt_stop_codon_exons, alt_is_premature):
    """Calculate exon-related features"""
    exon_info = transcript_exon_info or []
    stop_exons = alt_stop_codon_exons or []
    
    total_exons = len(exon_info)
    
    if not alt_is_premature or not stop_exons or not exon_info:
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


def _calculate_ptc_to_start_distance(alt_is_premature, alt_start_codon_pos, alt_first_stop_pos):
    """Calculate distance from PTC to start codon"""
    if not alt_is_premature:
        return None
    
    start = alt_start_codon_pos
    stop = alt_first_stop_pos
    
    if start is None or stop is None:
        return None
    
    if stop <= start:
        return None
    
    return stop - start


def _calculate_ptc_exon_length(alt_is_premature, alt_stop_codon_exons, transcript_exon_info):
    """Return the length of the exon containing the first premature stop codon (PTC)"""
    if not alt_is_premature:
        return None
    
    stop_exons = alt_stop_codon_exons or []
    exon_info = transcript_exon_info or []
    
    if not stop_exons or not exon_info:
        return None
    
    # First PTC exon = smallest exon number (transcript-order, strand-corrected)
    ptc_exon = min(int(e) for e in stop_exons)
    
    exon_dict = {int(e): int(length) for e, length in exon_info}
    return exon_dict.get(ptc_exon)


def _calculate_stop_codon_dist(ref_first_stop_pos, alt_first_stop_pos):
    """
    Calculate the distance between the reference stop codon and the alternative stop codon.
    Positive means the PTC is upstream of the reference stop codon.
    """
    ref_stop = ref_first_stop_pos
    alt_stop = alt_first_stop_pos
    
    if ref_stop is None or alt_stop is None:
        return None
    
    return ref_stop - alt_stop


def _calculate_ptc_to_downstream_ej(alt_is_premature, alt_stop_codon_exons, alt_cds_info, alt_first_stop_pos):
    """Calculate distance from PTC to the downstream exon junction"""
    if not alt_is_premature:
        return None
    
    stop_exons = alt_stop_codon_exons or []
    exon_info = alt_cds_info or []
    ptc_pos = alt_first_stop_pos
    
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


def _add_likely_misannotated_flag(cds_in_transcript, ref_start_codon_pos, ref_valid_stop):
    """Flag rows that look inconsistent between CDS and transcript annotations"""
    # if any of these are missing entirely, flag as likely misannotated
    if cds_in_transcript is None or ref_start_codon_pos is None or ref_valid_stop is None:
        return True
    
    flag = (
        (cds_in_transcript is False) or
        ((ref_start_codon_pos is not None) and (ref_start_codon_pos != 0)) or
        (ref_valid_stop is False)
    )
    
    return flag
