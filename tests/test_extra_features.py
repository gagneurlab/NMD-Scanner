"""
The transcript drawings show a transcript 5' to 3' and are not to scale. `u` is UTR, `=` is CDS, `[...]` is an exon,
`|` between two exons is an exon junction, and `*` is the PTC. The numbers under a transcript are positions in CDS
coordinates, as alt_first_stop_pos: 0 is the first base of the start codon. A row labelled tx gives transcript
coordinates instead. The numbers under an alt transcript are in alt CDS coordinates. `*--->|` is the distance from the
PTC to an exon junction or to the transcript end, and `<--->` is a length. "Technical Notes.md" defines the features
and the NMD escape rules with figures in the same style.
"""

import pandas as pd
import pytest

from nmd_scanner.extra_features import (
    add_likely_misannotated_flag,
    calculate_exon_features,
    calculate_ptc_exon_length,
    calculate_ptc_to_downstream_ej,
    calculate_ptc_to_start_distance,
    calculate_stop_codon_dist,
    calculate_utr_lengths,
    nmd_model_status,
)
from nmd_scanner.rules import analyze_sequence
from nmd_scanner.schema import MODEL_INPUTS, MODEL_STATUSES


def test_calculate_utr_lengths():
    # Example 1: - strand, CDS spans exon 1 to 8 (TXNL1). Exon numbers follow transcript order.
    # Exon 1 has 250 - 98 = 152 nt of 5'UTR, exon 8 has 5848 - 30 = 5818 nt of 3'UTR.
    row1 = {
        "strand": "-",
        "has_stop_codon": True,
        "ref_cds_exons": [
            {"exon_number": 8, "length": 30},
            {"exon_number": 7, "length": 105},
            {"exon_number": 6, "length": 173},
            {"exon_number": 5, "length": 70},
            {"exon_number": 4, "length": 123},
            {"exon_number": 3, "length": 174},
            {"exon_number": 2, "length": 97},
            {"exon_number": 1, "length": 98},
        ],
        "transcript_exons": [
            {"exon_number": "1", "length": 250},
            {"exon_number": "2", "length": 97},
            {"exon_number": "3", "length": 174},
            {"exon_number": "4", "length": 123},
            {"exon_number": "5", "length": 70},
            {"exon_number": "6", "length": 173},
            {"exon_number": "7", "length": 105},
            {"exon_number": "8", "length": 5848},
        ],
        "cds_start_in_transcript": 152,
        "cds_end_in_transcript": 152 + 870,
    }
    result1 = calculate_utr_lengths(row1)
    assert result1["utr5_length"] == 152
    assert result1["utr3_length"] == 5818

    # Example 2: + strand, CDS from exon 3 to 5
    row2 = {
        "strand": "+",
        "has_stop_codon": True,
        "ref_cds_exons": [
            {"exon_number": 3, "length": 50},
            {"exon_number": 4, "length": 120},
            {"exon_number": 5, "length": 80},
        ],
        "transcript_exons": [
            {"exon_number": "1", "length": 200},
            {"exon_number": "2", "length": 150},
            {"exon_number": "3", "length": 100},
            {"exon_number": "4", "length": 120},
            {"exon_number": "5", "length": 80},
            {"exon_number": "6", "length": 300},
        ],
        "cds_start_in_transcript": 400,
        "cds_end_in_transcript": 650,
    }
    result2 = calculate_utr_lengths(row2)
    # 5'UTR: exons 1 & 2 (200 + 150) + exon 3 (100 - 50)
    # 3'UTR: exon 6 (300)
    assert result2["utr5_length"] == 200 + 150 + 50
    assert result2["utr3_length"] == 300

    # Example 3: single exon, CDS fully inside it
    row3 = {
        "strand": "-",
        "has_stop_codon": True,
        "ref_cds_exons": [{"exon_number": 1, "length": 60}],
        "transcript_exons": [{"exon_number": "1", "length": 150}],
        "cds_start_in_transcript": 40,
        "cds_end_in_transcript": 100,
    }
    result3 = calculate_utr_lengths(row3)
    assert result3["utr5_length"] == 40
    assert result3["utr3_length"] == 50

    # Example 4: no annotated stop codon (cds_end_NF): the 3'UTR starts after the stop codon, so its length is unknown
    row = {**row2, "has_stop_codon": False}
    result = calculate_utr_lengths(row)
    assert result["utr5_length"] == 400
    assert result["utr3_length"] is None
    # single exon
    row = {**row3, "has_stop_codon": False}
    result = calculate_utr_lengths(row)
    assert result["utr5_length"] == 40
    assert result["utr3_length"] is None

    # Example 5: missing information
    # Example 5.1: CDS position in the transcript unknown
    row = {
        "strand": "+",
        "has_stop_codon": True,
        "ref_cds_exons": [{"exon_number": 1, "length": 100}, {"exon_number": 2, "length": 150}],
        "transcript_exons": [{"exon_number": "1", "length": 200}, {"exon_number": "2", "length": 300}],
        "cds_start_in_transcript": None,
        "cds_end_in_transcript": None,
    }
    result = calculate_utr_lengths(row)
    assert result["utr5_length"] is None
    assert result["utr3_length"] is None
    # Example 5.2: missing transcript_exons
    row = {
        "strand": "-",
        "has_stop_codon": True,
        "ref_cds_exons": [{"exon_number": 1, "length": 100}, {"exon_number": 2, "length": 150}],
        "cds_start_in_transcript": 0,
        "cds_end_in_transcript": 250,
    }
    result = calculate_utr_lengths(row)
    assert result["utr5_length"] is None
    assert result["utr3_length"] is None


def test_calculate_utr_lengths_cds_inside_one_exon():
    # Exons of 100/300/100 nt; the CDS with stop codon lies inside exon 2, at transcript positions 150 to 330.
    # The non-CDS part of exon 2 splits into 50 nt of 5'UTR and 70 nt of 3'UTR.
    row = {
        "strand": "+",
        "has_stop_codon": True,
        "ref_cds_exons": [{"exon_number": 2, "length": 180}],
        "transcript_exons": [
            {"exon_number": "1", "length": 100},
            {"exon_number": "2", "length": 300},
            {"exon_number": "3", "length": 100},
        ],
        "cds_start_in_transcript": 150,
        "cds_end_in_transcript": 330,
    }
    result = calculate_utr_lengths(row)
    assert result["utr5_length"] == 150
    assert result["utr3_length"] == 170


# Exons of a TCGA example (TXNL1), in transcript order. The CDS starts at transcript position 152, in exon 1.
TXNL1_EXONS = [
    {"exon_number": 1, "length": 250},
    {"exon_number": 2, "length": 97},
    {"exon_number": 3, "length": 174},
    {"exon_number": 4, "length": 123},
    {"exon_number": 5, "length": 70},
    {"exon_number": 6, "length": 173},
    {"exon_number": 7, "length": 105},
    {"exon_number": 8, "length": 5848},
]


def ptc_row(ptc_pos, exons=TXNL1_EXONS, alt_exons=None, cds_start=152):
    """A PTC row without a start loss: the PTC at alt CDS position ``ptc_pos``, and the alt exons ``alt_exons``."""

    return {
        "alt_has_ptc": True,
        "alt_first_stop_pos": ptc_pos,
        "alt_cds_start_in_transcript": cds_start,
        "transcript_exons": exons,
        "alt_transcript_exons": exons if alt_exons is None else alt_exons,
    }


def test_calculate_exon_features():
    # PTC in exon 1, at transcript position 152 + 50 = 202
    assert calculate_exon_features(ptc_row(50)) == {
        "total_exon_count": 8,
        "upstream_exon_count": 0,
        "downstream_exon_count": 7,
    }

    # PTC in exon 2, at transcript position 152 + 100 = 252: exon 1 ends at 250
    assert calculate_exon_features(ptc_row(100)) == {
        "total_exon_count": 8,
        "upstream_exon_count": 1,
        "downstream_exon_count": 6,
    }

    # PTC in the last exon, at transcript position 152 + 900 = 1052: exon 7 ends at 992
    assert calculate_exon_features(ptc_row(900)) == {
        "total_exon_count": 8,
        "upstream_exon_count": 7,
        "downstream_exon_count": 0,
    }

    # Single exon transcript
    row = ptc_row(90, exons=[{"exon_number": 1, "length": 500}], cds_start=100)
    assert calculate_exon_features(row) == {"total_exon_count": 1, "upstream_exon_count": 0, "downstream_exon_count": 0}

    # After a start loss, the PTC lies at alt_scan_first_stop_pos, in alt transcript positions: 260 is in exon 2
    row = {**ptc_row(None), "start_loss": True, "alt_scan_first_stop_pos": 260}
    assert calculate_exon_features(row) == {"total_exon_count": 8, "upstream_exon_count": 1, "downstream_exon_count": 6}

    # The counts take the alt exons. After a 3 nt deletion in exon 1, transcript position 152 + 95 = 247 is the first
    # base of exon 2. With the ref exon lengths, it would lie in exon 1.
    alt_exons = [{"exon_number": 1, "length": 247}, *TXNL1_EXONS[1:]]
    assert calculate_exon_features(ptc_row(95, alt_exons=alt_exons)) == {
        "total_exon_count": 8,
        "upstream_exon_count": 1,
        "downstream_exon_count": 6,
    }

    # An exon that the variant deletes (length 0) is not in the mRNA, so it is no upstream or downstream exon. It
    # still counts in total_exon_count, the exons of the ref transcript.
    alt_exons = [
        {"exon_number": 1, "length": 250},
        {"exon_number": 2, "length": 97},
        {"exon_number": 3, "length": 0},
        *TXNL1_EXONS[3:],
    ]
    assert calculate_exon_features(ptc_row(100, alt_exons=alt_exons)) == {
        "total_exon_count": 8,
        "upstream_exon_count": 1,
        "downstream_exon_count": 5,
    }

    # Without the PTC position, or without the alt exons (no alt transcript): no counts
    for row in [ptc_row(None), {**ptc_row(50), "alt_transcript_exons": None}]:
        assert calculate_exon_features(row) == {
            "total_exon_count": 8,
            "upstream_exon_count": None,
            "downstream_exon_count": None,
        }

    # Not a PTC row
    row = {**ptc_row(50), "alt_has_ptc": False}
    assert calculate_exon_features(row) == {
        "total_exon_count": 8,
        "upstream_exon_count": None,
        "downstream_exon_count": None,
    }

    # Transcript without exons
    row = {"alt_has_ptc": False, "transcript_exons": [], "alt_transcript_exons": None}
    assert calculate_exon_features(row) == {
        "total_exon_count": None,
        "upstream_exon_count": None,
        "downstream_exon_count": None,
    }


def test_calculate_ptc_to_start_distance():
    # The alt CDS starts with the annotated start codon, at CDS position 0, so the distance is alt_first_stop_pos
    start = {"has_start_codon": True, "ref_cds_seq": "ATGAAATAA", "alt_cds_seq": "ATGTAATAA"}
    row = {**start, "alt_has_ptc": True, "alt_first_stop_pos": 215}
    assert calculate_ptc_to_start_distance(row) == 215

    # The PTC is the annotated start codon itself, a stop codon such as TAG → None
    row2 = {**start, "alt_has_ptc": True, "alt_first_stop_pos": 0}
    assert calculate_ptc_to_start_distance(row2) is None

    # PTC not premature → None
    row4 = {**start, "alt_has_ptc": False, "alt_first_stop_pos": 300}
    assert calculate_ptc_to_start_distance(row4) is None

    # Missing alt_first_stop_pos → None
    row5 = {**start, "alt_has_ptc": True}
    assert calculate_ptc_to_start_distance(row5) is None

    # Without has_start_codon and alt_cds_seq, the alt CDS has no annotated start codon → None
    row6 = {
        "alt_has_ptc": True,
        "alt_first_stop_pos": 200,
    }
    assert calculate_ptc_to_start_distance(row6) is None

    # All fields missing → None
    row7 = {}
    assert calculate_ptc_to_start_distance(row7) is None


def test_calculate_ptc_exon_length():
    # PTC in exon 1, at transcript position 152 + 50 = 202
    assert calculate_ptc_exon_length(ptc_row(50)) == 250

    # PTC in exon 2, at transcript position 152 + 100 = 252
    assert calculate_ptc_exon_length(ptc_row(100)) == 97

    # The alt exons locate the PTC: after a 3 nt deletion in exon 1, transcript position 152 + 95 = 247 is the first
    # base of exon 2
    alt_exons = [{"exon_number": 1, "length": 247}, *TXNL1_EXONS[1:]]
    assert calculate_ptc_exon_length(ptc_row(95, alt_exons=alt_exons)) == 97

    # The length is the one in the alt transcript: a 3 nt deletion in exon 2 shortens it to 94 nt
    alt_exons = [{"exon_number": 1, "length": 250}, {"exon_number": 2, "length": 94}, *TXNL1_EXONS[2:]]
    assert calculate_ptc_exon_length(ptc_row(100, alt_exons=alt_exons)) == 94

    # Not a PTC row, no PTC position
    assert calculate_ptc_exon_length({**ptc_row(100), "alt_has_ptc": False}) is None
    assert calculate_ptc_exon_length(ptc_row(None)) is None


def _analyzed(ref_cds_seq, alt_cds_seq):
    """analyze_sequence row of a single exon CDS with an annotated stop codon."""
    df = pd.DataFrame(
        [
            {
                "ref_cds_seq": ref_cds_seq,
                "alt_cds_seq": alt_cds_seq,
                "ref_cds_length": len(ref_cds_seq),
                "alt_cds_length": len(alt_cds_seq),
                "ref_cds_exons": [{"exon_number": 1, "length": len(ref_cds_seq)}],
                "alt_cds_exons": [{"exon_number": 1, "length": len(alt_cds_seq)}],
                "has_start_codon": True,
                "has_stop_codon": True,
                "cds_frame": 0,
            }
        ]
    )
    return analyze_sequence(df).iloc[0]


def test_has_stop_codon_is_required():
    row = {
        "strand": "+",
        "ref_cds_exons": [{"exon_number": 1, "length": 60}],
        "transcript_exons": [{"exon_number": "1", "length": 100}],
    }
    with pytest.raises(KeyError, match="has_stop_codon"):
        calculate_utr_lengths(row)
    row = {"alt_cds_length": 903, "alt_first_stop_pos": 900, "alt_has_ptc": False}
    with pytest.raises(KeyError, match="has_stop_codon"):
        calculate_stop_codon_dist(row)


def test_calculate_stop_codon_dist():
    # Positions are in alt CDS coordinates: the reference stop codon is the last codon of the alt CDS.
    # Case 1: PTC upstream of reference stop
    row1 = {"alt_cds_length": 1003, "alt_first_stop_pos": 800, "alt_has_ptc": True, "has_stop_codon": True}
    assert calculate_stop_codon_dist(row1) == 200

    # Case 2: no PTC, the first stop codon of the alt is the reference stop codon
    row2 = {"alt_cds_length": 903, "alt_first_stop_pos": 900, "alt_has_ptc": False, "has_stop_codon": True}
    assert calculate_stop_codon_dist(row2) == 0

    # Case 3: Missing alt stop codon
    row3 = {"alt_cds_length": 903, "alt_first_stop_pos": None, "alt_has_ptc": False, "has_stop_codon": True}
    assert calculate_stop_codon_dist(row3) is None

    # Case 4: no annotated stop codon (cds_end_NF): there is no reference stop codon
    row4 = {"alt_cds_length": 903, "alt_first_stop_pos": 600, "alt_has_ptc": True, "has_stop_codon": False}
    assert calculate_stop_codon_dist(row4) is None

    # Case 5: internal in-frame TGA in the reference (selenocysteine): the reference stop codon is the annotated TAA
    #         ATG AAA TGA AAA CCC AAA TAA, PTC from AAA>TAA in codon 2
    row5 = _analyzed("ATGAAATGAAAACCCAAATAA", "ATGTAATGAAAACCCAAATAA")
    assert row5["ref_first_stop_pos"] == 6
    assert calculate_stop_codon_dist(row5) == 15

    # Case 6: frameshift deletion upstream of the PTC: deleting the C of CTG shifts the PTC by -1 in the alt CDS
    #         ref ATG AAA CTG ACC CCC TAA, alt ATG AAA TGA CCC CCT AA: PTC at ref position 7, alt position 6
    row6 = _analyzed("ATGAAACTGACCCCCTAA", "ATGAAATGACCCCCTAA")
    assert row6["alt_first_stop_pos"] == 6
    assert calculate_stop_codon_dist(row6) == 8

    # Case 7: TAA>TGAA, an insertion inside the stop codon: TGA sits at the position of the reference stop codon,
    #         although the alt CDS is 1 nt longer. Without the alt transcript or the CDS position in the transcript, the
    #         last codon of the alt CDS stands in.
    row7 = {
        "alt_cds_length": 13,
        "alt_first_stop_pos": 9,
        "has_stop_codon": True,
        "transcript_seq": "CCATGAAACCCTAAGG",
        "alt_transcript_seq": "CCATGAAACCCTGAAGG",
        "cds_start_in_transcript": 2,
        "alt_cds_start_in_transcript": 2,
        "cds_frame": 0,
        "cds_end_in_transcript": 14,
    }
    assert calculate_stop_codon_dist(row7) == 0
    assert calculate_stop_codon_dist({**row7, "alt_transcript_seq": None}) == 1
    assert (
        calculate_stop_codon_dist({**row7, "cds_start_in_transcript": None, "alt_cds_start_in_transcript": None}) == 1
    )

    # Case 8: TA>T deletes one A of the stop codon TAA before a 3'UTR A. The alt CDS ends in TA and has no stop codon,
    #         but the alt transcript still reads TAA at the annotated position.
    row8 = {
        "alt_cds_length": 11,
        "alt_first_stop_pos": None,
        "has_stop_codon": True,
        "transcript_seq": "CCATGAAATGGTAAACTGG",
        "alt_transcript_seq": "CCATGAAATGGTAACTGG",
        "cds_start_in_transcript": 2,
        "alt_cds_start_in_transcript": 2,
        "cds_frame": 0,
        "cds_end_in_transcript": 14,
    }
    assert calculate_stop_codon_dist(row8) == 0


def test_calculate_ptc_to_downstream_ej():
    # Case 1: PTC in exon 2, simple transcript
    #           100 nt            200 nt              150 nt
    #     5' [==========]|[===============*====]|[===============] 3'
    #        0           100              250   300              450
    #                                     *---->|  ptc_to_exon_end = 50
    row1 = {
        "alt_has_ptc": True,
        "alt_first_stop_pos": 250,  # PTC position
        "alt_transcript_exons": [
            {"exon_number": "1", "length": 100},
            {"exon_number": "2", "length": 200},
            {"exon_number": "3", "length": 150},
        ],
        "alt_cds_start_in_transcript": 0,
    }
    # End of exon 2: 100 + 200 = 300, distance = 300 - 250 = 50
    assert calculate_ptc_to_downstream_ej(row1) == 50

    # Case 2: PTC in first exon
    #           100 nt            200 nt              150 nt
    #     5' [======*===]|[====================]|[===============] 3'
    #        0      60   100                    300              450
    #               *--->|  ptc_to_exon_end = 40
    row2 = {
        "alt_has_ptc": True,
        "alt_first_stop_pos": 60,
        "alt_transcript_exons": [
            {"exon_number": "1", "length": 100},
            {"exon_number": "2", "length": 200},
            {"exon_number": "3", "length": 150},
        ],
        "alt_cds_start_in_transcript": 0,
    }
    # End of exon 1: 100, distance = 100 - 60 = 40
    assert calculate_ptc_to_downstream_ej(row2) == 40

    # Case 3: PTC in the last exon, which ends at the transcript end
    #           100 nt            200 nt                 200 nt
    #     5' [==========]|[====================]|[=============*======] 3'
    #        0           100                    300            430    500
    #                                                          *----->|  ptc_to_exon_end = 70
    row3 = {
        "alt_has_ptc": True,
        "alt_first_stop_pos": 430,
        "alt_transcript_exons": [
            {"exon_number": "1", "length": 100},
            {"exon_number": "2", "length": 200},
            {"exon_number": "3", "length": 200},
        ],
        "alt_cds_start_in_transcript": 0,
    }
    # The transcript ends at 100 + 200 + 200 = 500, distance = 500 - 430 = 70
    assert calculate_ptc_to_downstream_ej(row3) == 70

    # Case 4: Not premature → should return None
    row4 = {**row1, "alt_has_ptc": False}
    assert calculate_ptc_to_downstream_ej(row4) is None

    # Case 5: PTC position missing → should return None
    row5 = {**row1, "alt_first_stop_pos": None}
    assert calculate_ptc_to_downstream_ej(row5) is None

    # Case 6: PTC in the last CDS exon, followed by a UTR-only exon
    #                200 nt            100 nt         60 nt          300 nt
    #     5' [uuuu================]|[==========]|[====*===uuuu]|[uuuuuuuuuuuuuuu] 3'
    #        -40  0                160          260   280 300  320              620
    #                                                 *------->|  ptc_to_exon_end = 40
    row6 = {
        "alt_has_ptc": True,
        "alt_first_stop_pos": 280,
        "alt_transcript_exons": [
            {"exon_number": "1", "length": 200},
            {"exon_number": "2", "length": 100},
            {"exon_number": "3", "length": 60},
            {"exon_number": "4", "length": 300},
        ],
        "alt_cds_start_in_transcript": 40,
    }
    # The CDS ends at 300, and exon 3 goes on with 20 nt of 3'UTR: the junction is at 320, distance = 320 - 280 = 40
    assert calculate_ptc_to_downstream_ej(row6) == 40

    # Case 7: same transcript, a 1 nt deletion in exon 1 moves the junctions 1 nt upstream in alt CDS coordinates
    #                200 nt                 100 nt              60 nt          300 nt
    # ref 5' [uuuu================]|[====================]|[========uuuu]|[uuuuuuuuuuuuuuu] 3'
    #        -40  0                160                    260       300  320
    # alt 5' [uuuu================]|[=================*==]|[========uuuu]|[uuuuuuuuuuuuuuu] 3'
    #        -40  0                159                250 259       299  319
    #                                                 *-->|  ptc_to_exon_end = 9
    row7 = {
        **row6,
        "alt_first_stop_pos": 250,
        "alt_transcript_exons": [
            {"exon_number": "1", "length": 199},
            {"exon_number": "2", "length": 100},
            {"exon_number": "3", "length": 60},
            {"exon_number": "4", "length": 300},
        ],
    }
    # Exon 2 ends at alt transcript position 299, the PTC lies at 40 + 250 = 290: distance = 299 - 290 = 9
    assert calculate_ptc_to_downstream_ej(row7) == 9

    # Case 8: same transcript, a 3 nt deletion in the 3'UTR of exon 3 moves the end of exon 3 3 nt upstream
    #                200 nt            100 nt         60 nt          300 nt
    # ref 5' [uuuu================]|[==========]|[====*===uuuu]|[uuuuuuuuuuuuuuu] 3'
    #        -40  0                160          260   280 300  320
    # alt 5' [uuuu================]|[==========]|[====*===uuu]|[uuuuuuuuuuuuuuu] 3'
    #        -40  0                160          260   280 300 317
    #                                                 *------>|  ptc_to_exon_end = 37
    row8 = {
        **row6,
        "alt_transcript_exons": [
            {"exon_number": "1", "length": 200},
            {"exon_number": "2", "length": 100},
            {"exon_number": "3", "length": 57},
            {"exon_number": "4", "length": 300},
        ],
    }
    # Exon 3 ends at alt transcript position 357, the PTC lies at 40 + 280 = 320: distance = 357 - 320 = 37
    assert calculate_ptc_to_downstream_ej(row8) == 37

    # Case 9: alt transcript exons unknown
    row9 = {**row1, "alt_transcript_exons": None}
    assert calculate_ptc_to_downstream_ej(row9) is None

    # Case 10: PTC in the last exon, which goes on with 60 nt of 3'UTR after the CDS
    #           100 nt            200 nt                    260 nt
    #     5' [==========]|[====================]|[=============*======uuuuuu] 3'
    #        0           100                    300            430    500   560
    #                                                          *----------->|  ptc_to_exon_end = 130
    #                                                          *----->|  70 to the CDS end, not measured
    row10 = {
        **row3,
        "alt_transcript_exons": [
            {"exon_number": "1", "length": 100},
            {"exon_number": "2", "length": 200},
            {"exon_number": "3", "length": 260},
        ],
    }
    # The transcript ends at 560, distance = 560 - 430 = 130: the 3'UTR that the PTC creates, not the 70 nt to the CDS end
    assert calculate_ptc_to_downstream_ej(row10) == 130

    # Case 11: same transcript, a 1 nt deletion in exon 1 moves the transcript end 1 nt upstream in alt CDS coordinates
    #           100 nt            200 nt                    260 nt
    # ref 5' [==========]|[====================]|[====================uuuuuu] 3'
    #        0           100                    300                   500   560
    # alt 5' [==========]|[====================]|[=============*======uuuuuu] 3'
    #        0           99                     299            429    499   559
    #                                                          *----------->|  ptc_to_exon_end = 130
    row11 = {
        **row10,
        "alt_first_stop_pos": 429,
        "alt_transcript_exons": [
            {"exon_number": "1", "length": 99},
            {"exon_number": "2", "length": 200},
            {"exon_number": "3", "length": 260},
        ],
    }
    # The transcript ends at 560 - 1 = 559, distance = 559 - 429 = 130
    assert calculate_ptc_to_downstream_ej(row11) == 130

    # Case 12: single exon transcript of 300 nt, with 50 nt of 5'UTR and a CDS of 200 nt
    #                     300 nt
    #     5' [uuuuu==========*=========uuuuu] 3'
    #        -50   0         100       200  250
    #                        *------------->|  ptc_to_exon_end = 150
    row12 = {
        "alt_has_ptc": True,
        "alt_first_stop_pos": 100,
        "alt_transcript_exons": [{"exon_number": "1", "length": 300}],
        "alt_cds_start_in_transcript": 50,
    }
    # The transcript ends at 300 - 50 = 250 in CDS coordinates, distance = 250 - 100 = 150
    assert calculate_ptc_to_downstream_ej(row12) == 150


def test_add_likely_misannotated_flag():

    # Baseline: "good" annotation, not misannotated
    start = {"has_start_codon": True, "ref_cds_seq": "ATGAAATAA"}
    row0 = {"cds_in_transcript": True, **start, "ref_valid_stop": True}
    assert add_likely_misannotated_flag(row0) is False

    # CDS not in transcript
    row1 = {"cds_in_transcript": False, **start, "ref_valid_stop": True}
    assert add_likely_misannotated_flag(row1) is True

    # The CDS does not start with an annotated start codon
    row2 = {"cds_in_transcript": True, **start, "has_start_codon": False, "ref_valid_stop": True}
    assert add_likely_misannotated_flag(row2) is True

    # No valid stop codon at the end of the CDS sequence
    row3 = {"cds_in_transcript": True, **start, "ref_valid_stop": False}
    assert add_likely_misannotated_flag(row3) is True

    # Multiple conditions that point to a likely misannotation
    row4 = {"cds_in_transcript": True, **start, "has_start_codon": False, "ref_valid_stop": False}
    assert add_likely_misannotated_flag(row4) is True

    # Missing information
    row5 = {}
    assert add_likely_misannotated_flag(row5) is True

    # A CDS of fewer than 3 nt holds no start codon
    row5 = {"cds_in_transcript": True, "has_start_codon": True, "ref_cds_seq": "AT", "ref_valid_stop": True}
    assert add_likely_misannotated_flag(row5) is True


def _status_row(**values):
    """A PTC row that the model can score: a new PTC and all model inputs set. ``values`` overrides columns."""

    row = {name: 1 for name in MODEL_INPUTS}
    row.update(
        unknown_reason=None,
        alt_has_ptc=True,
        ref_has_ptc=False,
        has_stop_codon=True,
        has_start_codon=True,
        start_loss=False,
    )
    row.update(values)
    return row


# The null model inputs of a row with an unknown alt transcript: 15 of the 19 model inputs
_UNKNOWN_INPUTS = {
    name: None
    for name in MODEL_INPUTS
    if name not in ("total_exon_count", "utr5_length", "utr3_length", "transcript_length")
}


@pytest.mark.parametrize(
    ("values", "status"),
    [
        ({}, "ok"),
        ({"unknown_reason": "splice_site_destroyed", "alt_has_ptc": None, **_UNKNOWN_INPUTS}, "unknown_effect"),
        ({"alt_has_ptc": False}, "no_ptc"),
        ({"alt_has_ptc": None}, "no_ptc"),
        ({"ref_has_ptc": True}, "ref_ptc"),
        ({"has_stop_codon": False, "annotated_stop_distance": None, "utr3_length": None}, "no_annotated_stop"),
        ({"has_start_codon": False, "ptc_to_start_codon": None}, "no_annotated_start"),
        ({"start_loss": True, "ptc_to_start_codon": None}, "start_lost"),
        ({"start_loss": True}, "ok"),
        ({"ptc_to_start_codon": None}, "missing_input"),
        ({"ptc_to_exon_end": None}, "missing_input"),
        ({"ref_has_ptc": None}, "ok"),
        # overlapping reasons: the first one in MODEL_STATUSES wins
        (
            {"unknown_reason": "exon_boundary_ambiguous", "alt_has_ptc": None, "has_stop_codon": False},
            "unknown_effect",
        ),
        ({"alt_has_ptc": False, "ref_has_ptc": True, "ptc_to_start_codon": None}, "no_ptc"),
        ({"ref_has_ptc": True, "has_stop_codon": False, "utr3_length": None}, "ref_ptc"),
        ({"has_stop_codon": False, "utr3_length": None, "ptc_to_start_codon": None}, "no_annotated_stop"),
        (
            {"has_stop_codon": False, "has_start_codon": False, "utr3_length": None, "ptc_to_start_codon": None},
            "no_annotated_stop",
        ),
        ({"has_start_codon": False, "ptc_to_start_codon": None, "ptc_to_exon_end": None}, "no_annotated_start"),
        ({"start_loss": True, "ptc_to_start_codon": None, "ptc_to_exon_end": None}, "start_lost"),
    ],
)
def test_nmd_model_status(values, status):
    table = pd.DataFrame([_status_row(), _status_row(**values)], index=[7, 3])

    result = nmd_model_status(table)

    assert result.tolist() == ["ok", status]
    assert result.index.tolist() == [7, 3]
    assert result.dtype == pd.StringDtype("python")


def test_nmd_model_status_values_are_the_documented_ones():
    assert MODEL_STATUSES == (
        "unknown_effect",
        "no_ptc",
        "ref_ptc",
        "no_annotated_stop",
        "no_annotated_start",
        "start_lost",
        "missing_input",
        "ok",
    )
