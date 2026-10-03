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
    add_nmd_features,
    calculate_exon_features,
    calculate_ptc_exon_length,
    calculate_ptc_to_downstream_ej,
    calculate_ptc_to_start_distance,
    calculate_stop_codon_dist,
    calculate_utr_lengths,
    evaluate_nmd_escape_rules,
)
from nmd_scanner.rules import analyze_sequence


def test_calculate_utr_lengths():
    # Example 1: - strand, CDS spans exon 1 to 8 (TXNL1). Exon numbers follow transcript order.
    # Exon 1 has 250 - 98 = 152 nt of 5'UTR, exon 8 has 5848 - 30 = 5818 nt of 3'UTR.
    row1 = {
        "strand": "-",
        "has_stop_codon": True,
        "ref_cds_info": [(8, 30), (7, 105), (6, 173), (5, 70), (4, 123), (3, 174), (2, 97), (1, 98)],
        "transcript_exon_info": [
            ("1", 250),
            ("2", 97),
            ("3", 174),
            ("4", 123),
            ("5", 70),
            ("6", 173),
            ("7", 105),
            ("8", 5848),
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
        "ref_cds_info": [(3, 50), (4, 120), (5, 80)],
        "transcript_exon_info": [("1", 200), ("2", 150), ("3", 100), ("4", 120), ("5", 80), ("6", 300)],
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
        "ref_cds_info": [(1, 60)],
        "transcript_exon_info": [("1", 150)],
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
        "ref_cds_info": [(1, 100), (2, 150)],
        "transcript_exon_info": [("1", 200), ("2", 300)],
        "cds_start_in_transcript": None,
        "cds_end_in_transcript": None,
    }
    result = calculate_utr_lengths(row)
    assert result["utr5_length"] is None
    assert result["utr3_length"] is None
    # Example 5.2: missing transcript_exon_info
    row = {
        "strand": "-",
        "has_stop_codon": True,
        "ref_cds_info": [(1, 100), (2, 150)],
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
        "ref_cds_info": [(2, 180)],
        "transcript_exon_info": [("1", 100), ("2", 300), ("3", 100)],
        "cds_start_in_transcript": 150,
        "cds_end_in_transcript": 330,
    }
    result = calculate_utr_lengths(row)
    assert result["utr5_length"] == 150
    assert result["utr3_length"] == 170


def test_calculate_exon_features():

    # Example from tcga
    row = {
        "alt_is_premature": True,
        "strand": "+",
        "alt_stop_codon_exons": [1, 1],
        "transcript_exon_info": [
            ("1", 250),
            ("2", 97),
            ("3", 174),
            ("4", 123),
            ("5", 70),
            ("6", 173),
            ("7", 105),
            ("8", 5848),
        ],
    }

    result = calculate_exon_features(row)
    assert result["total_exon_count"] == 8
    assert result["upstream_exon_count"] == 0
    assert result["downstream_exon_count"] == 7

    # Example 1: + strand, PTC in middle exon
    row1 = {
        "alt_is_premature": True,
        "strand": "+",
        "transcript_exon_info": [("1", 100), ("2", 150), ("3", 120), ("4", 110)],
        "alt_stop_codon_exons": [2],
    }
    result1 = calculate_exon_features(row1)
    assert result1 == {"total_exon_count": 4, "upstream_exon_count": 1, "downstream_exon_count": 2}

    # Example 2: - strand, PTC in exon 2 (which is second in reverse)
    row2 = {
        "alt_is_premature": True,
        "strand": "-",
        "transcript_exon_info": [("1", 100), ("2", 150), ("3", 120), ("4", 110)],
        "alt_stop_codon_exons": [2],
    }
    result2 = calculate_exon_features(row2)
    assert result2 == {
        "total_exon_count": 4,
        "upstream_exon_count": 1,  # reversed: [4,3,2,1] → index 1
        "downstream_exon_count": 2,
    }

    # Example 3: PTC in first exon on + strand
    row3 = {
        "alt_is_premature": True,
        "strand": "+",
        "transcript_exon_info": [("1", 100), ("2", 150), ("3", 120)],
        "alt_stop_codon_exons": [1],
    }
    result3 = calculate_exon_features(row3)
    assert result3 == {"total_exon_count": 3, "upstream_exon_count": 0, "downstream_exon_count": 2}

    # Example 4: PTC in last exon on - strand
    row4 = {
        "alt_is_premature": True,
        "strand": "-",
        "transcript_exon_info": [("1", 100), ("2", 150), ("3", 120)],
        "alt_stop_codon_exons": [1],
    }
    result4 = calculate_exon_features(row4)
    assert result4 == {"total_exon_count": 3, "upstream_exon_count": 0, "downstream_exon_count": 2}

    # Example 5: Single exon transcript
    row5 = {"alt_is_premature": True, "strand": "+", "transcript_exon_info": [("1", 500)], "alt_stop_codon_exons": [1]}
    result5 = calculate_exon_features(row5)
    assert result5 == {"total_exon_count": 1, "upstream_exon_count": 0, "downstream_exon_count": 0}

    # Example 7: PTC exon not present in transcript
    row6 = {
        "alt_is_premature": True,
        "strand": "-",
        "transcript_exon_info": [("1", 100), ("2", 200)],
        "alt_stop_codon_exons": [99],
    }
    result6 = calculate_exon_features(row6)
    assert result6 == {"total_exon_count": 2, "upstream_exon_count": None, "downstream_exon_count": None}

    # Example 7: Missing stop codon exons
    row7 = {
        "alt_is_premature": True,
        "strand": "+",
        "transcript_exon_info": [("1", 100), ("2", 200)],
        "alt_stop_codon_exons": [],
    }
    result7 = calculate_exon_features(row7)
    assert result7 == {"total_exon_count": 2, "upstream_exon_count": None, "downstream_exon_count": None}

    # Example 8: is not PTC
    row8 = {
        "alt_is_premature": False,
        "strand": "-",
        "transcript_exon_info": [("1", 100), ("2", 200)],
        "alt_stop_codon_exons": [99],
    }
    result8 = calculate_exon_features(row8)
    assert result8 == {"total_exon_count": 2, "upstream_exon_count": None, "downstream_exon_count": None}

    # Example 9: Transcript with no exons
    row9 = {"alt_is_premature": False, "strand": "+", "transcript_exon_info": [], "alt_stop_codon_exons": []}
    result9 = calculate_exon_features(row9)
    assert result9 == {"total_exon_count": None, "upstream_exon_count": None, "downstream_exon_count": None}


def test_calculate_ptc_to_start_distance():
    row = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 215,
        "alt_start_codon_pos": 50,
    }
    assert calculate_ptc_to_start_distance(row) == 165

    # PTC before start codon
    row2 = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 30,
        "alt_start_codon_pos": 100,
    }
    assert calculate_ptc_to_start_distance(row2) is None

    # Same position → distance 0 # can not happen
    row3 = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 120,
        "alt_start_codon_pos": 120,
    }
    assert calculate_ptc_to_start_distance(row3) is None

    # PTC not premature → None
    row4 = {
        "alt_is_premature": False,
        "alt_first_stop_pos": 300,
        "alt_start_codon_pos": 200,
    }
    assert calculate_ptc_to_start_distance(row4) is None

    # Missing alt_first_stop_pos → None
    row5 = {
        "alt_is_premature": True,
        "alt_start_codon_pos": 200,
    }
    assert calculate_ptc_to_start_distance(row5) is None

    # Missing alt_start_codon_pos → None
    row6 = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 200,
    }
    assert calculate_ptc_to_start_distance(row6) is None

    # All fields missing → None
    row7 = {}
    assert calculate_ptc_to_start_distance(row7) is None


def test_calculate_ptc_exon_length():
    row = {
        "alt_is_premature": True,
        "alt_stop_codon_exons": [1, 1],
        "transcript_exon_info": [
            ("1", 250),
            ("2", 97),
            ("3", 174),
            ("4", 123),
            ("5", 70),
            ("6", 173),
            ("7", 105),
            ("8", 5848),
        ],
    }
    assert calculate_ptc_exon_length(row) == 250

    row2 = {
        "alt_is_premature": True,
        "alt_stop_codon_exons": [2, 3, 4],
        "transcript_exon_info": [
            ("1", 250),
            ("2", 97),
            ("3", 174),
            ("4", 123),
            ("5", 70),
            ("6", 173),
            ("7", 105),
            ("8", 5848),
        ],
    }
    assert calculate_ptc_exon_length(row2) == 97

    row3 = {
        "alt_is_premature": False,
        "alt_stop_codon_exons": [2, 3, 4],
        "transcript_exon_info": [
            ("1", 250),
            ("2", 97),
            ("3", 174),
            ("4", 123),
            ("5", 70),
            ("6", 173),
            ("7", 105),
            ("8", 5848),
        ],
    }
    assert calculate_ptc_exon_length(row3) is None

    row4 = {
        "alt_is_premature": True,
        "alt_stop_codon_exons": [],
        "transcript_exon_info": [
            ("1", 250),
            ("2", 97),
            ("3", 174),
            ("4", 123),
            ("5", 70),
            ("6", 173),
            ("7", 105),
            ("8", 5848),
        ],
    }
    assert calculate_ptc_exon_length(row4) is None


def _analyzed(ref_cds_seq, alt_cds_seq):
    """analyze_sequence row of a single exon CDS with an annotated stop codon."""
    df = pd.DataFrame(
        [
            {
                "ref_cds_seq": ref_cds_seq,
                "alt_cds_seq": alt_cds_seq,
                "ref_cds_len": len(ref_cds_seq),
                "alt_cds_len": len(alt_cds_seq),
                "ref_cds_info": [(1, len(ref_cds_seq))],
                "alt_cds_info": [(1, len(alt_cds_seq))],
                "has_stop_codon": True,
                "cds_frame": 0,
            }
        ]
    )
    return analyze_sequence(df).iloc[0]


def test_has_stop_codon_is_required():
    row = {"strand": "+", "ref_cds_info": [(1, 60)], "transcript_exon_info": [("1", 100)]}
    with pytest.raises(KeyError, match="has_stop_codon"):
        calculate_utr_lengths(row)
    row = {"alt_cds_len": 903, "alt_first_stop_pos": 900, "alt_is_premature": False}
    with pytest.raises(KeyError, match="has_stop_codon"):
        calculate_stop_codon_dist(row)


def test_calculate_stop_codon_dist():
    # Positions are in alt CDS coordinates: the reference stop codon is the last codon of the alt CDS.
    # Case 1: PTC upstream of reference stop
    row1 = {"alt_cds_len": 1003, "alt_first_stop_pos": 800, "alt_is_premature": True, "has_stop_codon": True}
    assert calculate_stop_codon_dist(row1) == 200

    # Case 2: no PTC, the first stop codon of the alt is the reference stop codon
    row2 = {"alt_cds_len": 903, "alt_first_stop_pos": 900, "alt_is_premature": False, "has_stop_codon": True}
    assert calculate_stop_codon_dist(row2) == 0

    # Case 3: Missing alt stop codon
    row3 = {"alt_cds_len": 903, "alt_first_stop_pos": None, "alt_is_premature": False, "has_stop_codon": True}
    assert calculate_stop_codon_dist(row3) is None

    # Case 4: no annotated stop codon (cds_end_NF): there is no reference stop codon
    row4 = {"alt_cds_len": 903, "alt_first_stop_pos": 600, "alt_is_premature": True, "has_stop_codon": False}
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
        "alt_cds_len": 13,
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
        "alt_cds_len": 11,
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


def test_evaluate_nmd_escape_rules():

    # Example 1: Single exon rule
    #                100 nt
    #     5' [=================*==] 3'
    #        0                 90 100
    #        <----------------->  90 nt from the start codon, < 150
    row2 = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 90,
        "alt_stop_codon_exons": [1],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [(1, 100)],
        "alt_cds_info": [(1, 100)],
        "total_exon_count": 1,
        "downstream_exon_count": 0,
        "ptc_exon_length": 100,
    }

    result = evaluate_nmd_escape_rules(row2)
    assert result["nmd_single_exon_rule"] == True
    assert result["nmd_last_exon_rule"] == True
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == True
    assert result["nmd_escape"] == True

    # Example 2: Last exon rule
    #           100 nt       100 nt       100 nt
    #     5' [==========]|[==========]|[=====*====] 3'
    #        0           100          200    250  300
    #                                        *  in the last exon, past the last exon junction at 200
    row = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 250,
        "alt_stop_codon_exons": [3],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [(1, 100), (2, 100), (3, 100)],
        "ref_cds_info": [(1, 100), (2, 100), (3, 100)],
        "alt_cds_info": [(1, 100), (2, 100), (3, 100)],
        "cds_start_in_transcript": 0,
        "total_exon_count": 3,
        "downstream_exon_count": 0,
        "ptc_exon_length": 100,
    }

    result = evaluate_nmd_escape_rules(row)
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == True
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == True

    # Example 3: 50nt from penultimate exon end
    #           100 nt       100 nt       100 nt
    #     5' [==========]|[======*===]|[==========] 3'
    #        0           100     160  200         300
    #                            *--->|  40 nt to the last exon junction, <= 50
    row = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 160,
        "alt_stop_codon_exons": [2],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [(1, 100), (2, 100), (3, 100)],
        "ref_cds_info": [(1, 100), (2, 100), (3, 100)],
        "alt_cds_info": [(1, 100), (2, 100), (3, 100)],
        "cds_start_in_transcript": 0,
        "total_exon_count": 3,
        "downstream_exon_count": 1,
        "ptc_exon_length": 100,
    }

    result = evaluate_nmd_escape_rules(row)
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == True
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == True

    # Example 4: Long exon rule
    #           100 nt              500 nt               100 nt
    #     5' [==========]|[==========*==============]|[==========] 3'
    #        0           100         300             600         700
    #                      <----------------------->  ptc_exon_length = 500, > 407
    row = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 300,
        "alt_stop_codon_exons": [2],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [(1, 100), (2, 500), (3, 100)],
        "ref_cds_info": [(1, 100), (2, 500), (3, 100)],
        "alt_cds_info": [(1, 100), (2, 500), (3, 100)],
        "cds_start_in_transcript": 0,
        "total_exon_count": 3,
        "downstream_exon_count": 1,
        "ptc_exon_length": 500,
    }

    result = evaluate_nmd_escape_rules(row)
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == True
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == True

    # Example 5: Start-proximal rule. The PTC lies at the last exon junction, not upstream of it: no 50 nt rule.
    # The PTC is the first base of exon 2, at 100.
    #           100 nt       100 nt
    #     5' [==========]|[*=========] 3'
    #        0           100         200
    #        <------------->  100 nt from the start codon, < 150
    row = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 100,
        "alt_stop_codon_exons": [2],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [(1, 100), (2, 100)],
        "ref_cds_info": [(1, 100), (2, 100)],
        "alt_cds_info": [(1, 100), (2, 100)],
        "cds_start_in_transcript": 0,
        "total_exon_count": 2,
        "downstream_exon_count": 0,
        "ptc_exon_length": 100,
    }

    result = evaluate_nmd_escape_rules(row)
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == True
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == True
    assert result["nmd_escape"] == True

    # Example 6: multiple escape rules. The PTC lies in the last exon, which is longer than 407 nt, and less than 150 nt
    # from the start codon. It lies past the last exon junction at 100, so the 50 nt rule does not fire.
    #         50 nt   50 nt            500 nt
    #     5' [=====]|[=====]|[==*=======================] 3'
    #        0      50      100 120                     600
    #        <------------------>  120 nt from the start codon, < 150
    #                         <------------------------>  ptc_exon_length = 500, > 407
    row = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 120,
        "alt_stop_codon_exons": [3],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [(1, 50), (2, 50), (3, 500)],
        "ref_cds_info": [(1, 50), (2, 50), (3, 500)],
        "alt_cds_info": [(1, 50), (2, 50), (3, 500)],
        "cds_start_in_transcript": 0,
        "total_exon_count": 3,
        "downstream_exon_count": 0,
        "ptc_exon_length": 500,
    }

    result = evaluate_nmd_escape_rules(row)
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == True
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == True
    assert result["nmd_start_proximal_rule"] == True
    assert result["nmd_escape"] == True

    # Example 7: Not premature → should skip
    row3 = {
        "alt_is_premature": False,
        "alt_first_stop_pos": 150,
        "alt_stop_codon_exons": [2],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [(1, 100), (2, 100)],
        "alt_cds_info": [(1, 100), (2, 100)],
        "total_exon_count": 2,
        "downstream_exon_count": 1,
        "ptc_exon_length": 100,
    }

    result = evaluate_nmd_escape_rules(row3)
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == False

    # Example 8: 50nt rule uses CDS-relative coordinates, not transcript-relative.
    # Transcript has a 200nt 5'UTR in exon 1; CDS spans only part of exon 1 and all of exons 2,3.
    # PTC at CDS pos 160 is in the last 50nt of the penultimate CDS exon (exon 2, CDS-end at 200).
    # Measured in transcript coordinates, without subtracting cds_start_in_transcript, the junction would lie at 400
    # and the rule would not fire.
    #                     300 nt                 100 nt            200 nt
    #     5' [uuuuuuuuuuuuuuuuuuuu==========]|[======*===]|[====================] 3'
    # tx     0                               300     360  400                   600
    # CDS    -200                 0          100     160  200                   400
    #                                                *--->|  40 nt to the last exon junction, <= 50
    row_cds_offset = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 160,
        "alt_stop_codon_exons": [2],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [(1, 300), (2, 100), (3, 200)],
        "ref_cds_info": [(1, 100), (2, 100), (3, 200)],
        "alt_cds_info": [(1, 100), (2, 100), (3, 200)],
        "cds_start_in_transcript": 200,
        "total_exon_count": 3,
        "downstream_exon_count": 1,
        "ptc_exon_length": 100,
    }

    result = evaluate_nmd_escape_rules(row_cds_offset)
    assert result["nmd_50nt_penultimate_rule"] == True
    assert result["nmd_escape"] == True

    # Example 9: PTC sits at CDS pos that would falsely fire without subtracting cds_start_in_transcript.
    # Transcript: exon1 300nt (200 UTR + 100 CDS), exon2 100nt CDS, exon3 200nt CDS.
    # Transcript-relative pen_end would be 400; CDS pos 360 is in [350, 400) → false positive.
    # CDS-relative pen_end is 200; CDS pos 360 is past it → rule must NOT fire.
    #                     300 nt                 100 nt            200 nt
    #     5' [uuuuuuuuuuuuuuuuuuuu==========]|[==========]|[================*===] 3'
    # tx     0                               300          400               560 600
    # CDS    -200                 0          100          200               360 400
    row_false_positive_guard = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 360,
        "alt_stop_codon_exons": [3],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [(1, 300), (2, 100), (3, 200)],
        "ref_cds_info": [(1, 100), (2, 100), (3, 200)],
        "alt_cds_info": [(1, 100), (2, 100), (3, 200)],
        "cds_start_in_transcript": 200,
        "total_exon_count": 3,
        "downstream_exon_count": 0,
        "ptc_exon_length": 200,
    }

    result = evaluate_nmd_escape_rules(row_false_positive_guard)
    assert result["nmd_50nt_penultimate_rule"] == False


def test_nmd_rules_with_utr_only_last_exon():
    # Exon numbers in transcript order. The CDS ends in exon 3; exon 4 holds only 3'UTR.
    # Exon 1: 200 nt (40 5'UTR + 160 CDS), exon 2: 100 nt CDS, exon 3: 60 nt (40 CDS with the stop codon + 20 3'UTR),
    # exon 4: 300 nt 3'UTR. In CDS coordinates, exon 2 ends at 260 and exon 3 ends at 320: the last exon junction.
    #                200 nt            100 nt         60 nt          300 nt
    #     5' [uuuu================]|[==========]|[========uuuu]|[uuuuuuuuuuuuuuu] 3'
    #        -40  0                160          260       300  320              620
    #                                       *-->|  PTC at 230: 30 nt to the end of exon 2
    #                                       *----------------->|  PTC at 230: 90 nt to the last exon junction, > 50
    #                                                 *------->|  PTC at 280: 40 nt, <= 50
    #                                               *--------->|  PTC at 270: 50 nt, <= 50
    #                                              *---------->|  PTC at 269: 51 nt, > 50
    transcript = {
        "alt_is_premature": True,
        "alt_start_codon_pos": 0,
        "has_stop_codon": True,
        "transcript_exon_info": [("1", 200), ("2", 100), ("3", 60), ("4", 300)],
        "ref_cds_info": [(1, 160), (2, 100), (3, 40)],
        "alt_cds_info": [(1, 160), (2, 100), (3, 40)],
        "cds_start_in_transcript": 40,
    }

    def evaluate(stop_pos, stop_exons):
        row = {**transcript, "alt_first_stop_pos": stop_pos, "alt_stop_codon_exons": stop_exons}
        row.update(add_nmd_features(row))
        return evaluate_nmd_escape_rules(row)

    # PTC 30 nt before the end of exon 2, but 90 nt before the last exon junction: no escape
    result = evaluate(230, [2, 3])
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_escape"] == False

    # PTC in exon 3, 40 nt before the last exon junction: escape by the 50 nt rule, not by the last exon rule
    result = evaluate(280, [3, 3])
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == True
    assert result["nmd_escape"] == True

    # The 50 nt rule includes its boundary: a PTC 50 nt before the last exon junction escapes, one 51 nt before it not
    assert evaluate(270, [3])["nmd_50nt_penultimate_rule"] == True
    assert evaluate(269, [3])["nmd_50nt_penultimate_rule"] == False


def test_nmd_features_with_cds_inside_one_exon():
    # Exons of 300/500 nt. Exon 1 holds 100 nt of 5'UTR, the whole CDS (180 nt with the stop codon) and 20 nt of 3'UTR.
    # Exon 2 holds only 3'UTR. In CDS coordinates, the last exon junction lies at 300 - 100 = 200.
    #                            300 nt                           500 nt
    #     5' [uuuuu================================*===uuuu]|[uuuuuuuuuuuuuuu] 3'
    #        -100  0                               165 180  200              700
    #                                              *------->|  ptc_to_intron = 35, the last exon junction: <= 50
    row = {
        "alt_is_premature": True,
        "alt_start_codon_pos": 0,
        "has_stop_codon": True,
        "alt_first_stop_pos": 165,
        "alt_stop_codon_exons": [1, 1],
        "transcript_exon_info": [("1", 300), ("2", 500)],
        "ref_cds_info": [(1, 180)],
        "alt_cds_info": [(1, 180)],
        "cds_start_in_transcript": 100,
        "cds_end_in_transcript": 280,
    }
    row.update(add_nmd_features(row))
    # PTC 35 nt before the last exon junction: escape by the 50 nt rule
    assert row["ptc_to_intron"] == 35
    assert evaluate_nmd_escape_rules(row)["nmd_50nt_penultimate_rule"] == True

    # Exons of 200/300 nt. The CDS (150 nt) ends exactly at the end of exon 1; exon 2 holds only 3'UTR.
    #                       200 nt                      300 nt
    #     5' [uuuuu===========================*==]|[uuuuuuuuuuuuuuu] 3'
    #        -50   0                          135 150              450
    #                                         *-->|  ptc_to_intron = 15, the last exon junction: <= 50
    row = {
        "alt_is_premature": True,
        "alt_start_codon_pos": 0,
        "has_stop_codon": True,
        "alt_first_stop_pos": 135,
        "alt_stop_codon_exons": [1, 1],
        "transcript_exon_info": [("1", 200), ("2", 300)],
        "ref_cds_info": [(1, 150)],
        "alt_cds_info": [(1, 150)],
        "cds_start_in_transcript": 50,
        "cds_end_in_transcript": 200,
    }
    row.update(add_nmd_features(row))
    # PTC 15 nt before the exon junction
    assert row["ptc_to_intron"] == 15
    assert evaluate_nmd_escape_rules(row)["nmd_50nt_penultimate_rule"] == True


def test_calculate_ptc_to_downstream_ej():
    # Case 1: PTC in exon 2, simple transcript
    #           100 nt            200 nt              150 nt
    #     5' [==========]|[===============*====]|[===============] 3'
    #        0           100              250   300              450
    #                                     *---->|  ptc_to_intron = 50
    row1 = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 250,  # PTC position
        "alt_stop_codon_exons": [2],
        "transcript_exon_info": [("1", 100), ("2", 200), ("3", 150)],
        "ref_cds_info": [(1, 100), (2, 200), (3, 150)],
        "alt_cds_info": [(1, 100), (2, 200), (3, 150)],
        "cds_start_in_transcript": 0,
    }
    # End of exon 2: 100 + 200 = 300, distance = 300 - 250 = 50
    assert calculate_ptc_to_downstream_ej(row1) == 50

    # Case 2: PTC in first exon
    #           100 nt            200 nt              150 nt
    #     5' [======*===]|[====================]|[===============] 3'
    #        0      60   100                    300              450
    #               *--->|  ptc_to_intron = 40
    row2 = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 60,
        "alt_stop_codon_exons": [1],
        "transcript_exon_info": [("1", 100), ("2", 200), ("3", 150)],
        "ref_cds_info": [(1, 100), (2, 200), (3, 150)],
        "alt_cds_info": [(1, 100), (2, 200), (3, 150)],
        "cds_start_in_transcript": 0,
    }
    # End of exon 1: 100, distance = 100 - 60 = 40
    assert calculate_ptc_to_downstream_ej(row2) == 40

    # Case 3: PTC in the last exon, which ends at the transcript end
    #           100 nt            200 nt                 200 nt
    #     5' [==========]|[====================]|[=============*======] 3'
    #        0           100                    300            430    500
    #                                                          *----->|  ptc_to_intron = 70
    row3 = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 430,
        "alt_stop_codon_exons": [3],
        "transcript_exon_info": [("1", 100), ("2", 200), ("3", 200)],
        "ref_cds_info": [(1, 100), (2, 200), (3, 200)],
        "alt_cds_info": [(1, 100), (2, 200), (3, 200)],
        "cds_start_in_transcript": 0,
    }
    # The transcript ends at 100 + 200 + 200 = 500, distance = 500 - 430 = 70
    assert calculate_ptc_to_downstream_ej(row3) == 70

    # Case 4: Multiple stop codons, take the smallest exon number
    #           100 nt            200 nt                 200 nt
    #     5' [==========]|[===============*====]|[====================] 3'
    #        0           100              250   300                   500
    #                                     *---->|  ptc_to_intron = 50
    #                                     *-------------------------->|  250 if exon 3 were the PTC exon
    row4 = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 250,
        "alt_stop_codon_exons": [2, 3],
        "transcript_exon_info": [("1", 100), ("2", 200), ("3", 200)],
        "ref_cds_info": [(1, 100), (2, 200), (3, 200)],
        "alt_cds_info": [(1, 100), (2, 200), (3, 200)],
        "cds_start_in_transcript": 0,
    }
    # PTC in exon 2, normal stop codon in exon 3. Smallest exon = 2, end of exon 2: 100 + 200 = 300,
    # distance = 300 - 250 = 50. Exon 3, the last exon, would give 500 - 250 = 250.
    assert calculate_ptc_to_downstream_ej(row4) == 50

    # Case 5: Not premature → should return None
    row5 = {**row1, "alt_is_premature": False}
    assert calculate_ptc_to_downstream_ej(row5) is None

    # Case 6: PTC position missing → should return None
    row6 = {**row1, "alt_first_stop_pos": None}
    assert calculate_ptc_to_downstream_ej(row6) is None

    # Case 7: PTC in the last CDS exon, followed by a UTR-only exon
    #                200 nt            100 nt         60 nt          300 nt
    #     5' [uuuu================]|[==========]|[====*===uuuu]|[uuuuuuuuuuuuuuu] 3'
    #        -40  0                160          260   280 300  320              620
    #                                                 *------->|  ptc_to_intron = 40
    row7 = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 280,
        "alt_stop_codon_exons": [3, 3],
        "transcript_exon_info": [("1", 200), ("2", 100), ("3", 60), ("4", 300)],
        "ref_cds_info": [(1, 160), (2, 100), (3, 40)],
        "alt_cds_info": [(1, 160), (2, 100), (3, 40)],
        "cds_start_in_transcript": 40,
    }
    # The CDS ends at 300, and exon 3 goes on with 20 nt of 3'UTR: the junction is at 320, distance = 320 - 280 = 40
    assert calculate_ptc_to_downstream_ej(row7) == 40

    # Case 8: same transcript, a 1 nt deletion in exon 1 moves the junctions 1 nt upstream in alt CDS coordinates
    #                200 nt                 100 nt              60 nt          300 nt
    # ref 5' [uuuu================]|[====================]|[========uuuu]|[uuuuuuuuuuuuuuu] 3'
    #        -40  0                160                    260       300  320
    # alt 5' [uuuu================]|[=================*==]|[========uuuu]|[uuuuuuuuuuuuuuu] 3'
    #        -40  0                159                250 259       299  319
    #                                                 *-->|  ptc_to_intron = 9
    row8 = {
        **row7,
        "alt_first_stop_pos": 250,
        "alt_stop_codon_exons": [2, 3],
        "alt_cds_info": [(1, 159), (2, 100), (3, 40)],
    }
    # Exon 2 ends at transcript position 300, i.e. at 300 - 40 - 1 = 259 in alt CDS coordinates: distance = 259 - 250 = 9
    assert calculate_ptc_to_downstream_ej(row8) == 9

    # Case 9: transcript exons unknown
    row9 = {**row1, "transcript_exon_info": None}
    assert calculate_ptc_to_downstream_ej(row9) is None

    # Case 10: PTC in the last exon, which goes on with 60 nt of 3'UTR after the CDS
    #           100 nt            200 nt                    260 nt
    #     5' [==========]|[====================]|[=============*======uuuuuu] 3'
    #        0           100                    300            430    500   560
    #                                                          *----------->|  ptc_to_intron = 130
    #                                                          *----->|  70 to the CDS end, not measured
    row10 = {**row3, "transcript_exon_info": [("1", 100), ("2", 200), ("3", 260)]}
    # The transcript ends at 560, distance = 560 - 430 = 130: the 3'UTR that the PTC creates, not the 70 nt to the CDS end
    assert calculate_ptc_to_downstream_ej(row10) == 130

    # Case 11: same transcript, a 1 nt deletion in exon 1 moves the transcript end 1 nt upstream in alt CDS coordinates
    #           100 nt            200 nt                    260 nt
    # ref 5' [==========]|[====================]|[====================uuuuuu] 3'
    #        0           100                    300                   500   560
    # alt 5' [==========]|[====================]|[=============*======uuuuuu] 3'
    #        0           99                     299            429    499   559
    #                                                          *----------->|  ptc_to_intron = 130
    row11 = {**row10, "alt_first_stop_pos": 429, "alt_cds_info": [(1, 99), (2, 200), (3, 200)]}
    # The transcript ends at 560 - 1 = 559, distance = 559 - 429 = 130
    assert calculate_ptc_to_downstream_ej(row11) == 130

    # Case 12: single exon transcript of 300 nt, with 50 nt of 5'UTR and a CDS of 200 nt
    #                     300 nt
    #     5' [uuuuu==========*=========uuuuu] 3'
    #        -50   0         100       200  250
    #                        *------------->|  ptc_to_intron = 150
    row12 = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 100,
        "alt_stop_codon_exons": [1],
        "transcript_exon_info": [("1", 300)],
        "ref_cds_info": [(1, 200)],
        "alt_cds_info": [(1, 200)],
        "cds_start_in_transcript": 50,
    }
    # The transcript ends at 300 - 50 = 250 in CDS coordinates, distance = 250 - 100 = 150
    assert calculate_ptc_to_downstream_ej(row12) == 150


def test_add_likely_misannotated_flag():

    # Baseline: "good" annotation, not misannotated
    row0 = {"cds_in_transcript": True, "ref_start_codon_pos": 0, "ref_valid_stop": True}
    assert add_likely_misannotated_flag(row0) is False

    # CDS not in transcript
    row1 = {"cds_in_transcript": False, "ref_start_codon_pos": 0, "ref_valid_stop": True}
    assert add_likely_misannotated_flag(row1) is True

    # Start position of CDS not at position 0
    row2 = {
        "cds_in_transcript": True,
        "ref_start_codon_pos": 14,  # start codon is not at beginning of CDS
        "ref_valid_stop": True,
    }
    assert add_likely_misannotated_flag(row2) is True

    # No valid stop codon at the end of the CDS sequence
    row3 = {"cds_in_transcript": True, "ref_start_codon_pos": 0, "ref_valid_stop": False}
    assert add_likely_misannotated_flag(row3) is True

    # Multiple conditions that point to a likely misannotation
    row4 = {"cds_in_transcript": True, "ref_start_codon_pos": 13, "ref_valid_stop": False}
    assert add_likely_misannotated_flag(row4) is True

    # Missing information
    row5 = {}
    assert add_likely_misannotated_flag(row5) is True

    # No start codon position
    row5 = {"cds_in_transcript": True, "ref_start_codon_pos": None, "ref_valid_stop": True}
    assert add_likely_misannotated_flag(row5) is True
