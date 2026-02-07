from nmd_scanner.annotation.utils import *
import pandas as pd

def test_calculate_utr_lengths():
    # Example 1: - strand, CDS spans exon 8 to 1
    result1 = calculate_utr_lengths(
        strand="-",
        ref_cds_info=[(8, 30), (7, 105), (6, 173), (5, 70), (4, 123), (3, 174), (2, 97), (1, 98)],
        transcript_exon_info=[('1', 250), ('2', 97), ('3', 174), ('4', 123), ('5', 70),
                                 ('6', 173), ('7', 105), ('8', 5848)],
        ref_cds_start=70925000,
        transcript_start=70920000,
        ref_cds_stop=70929273,
        transcript_end=70950000,
    )
    # Exon 1 has 250 - 98 = 152 nt of 3'UTR
    # Exon 8 has 5848 - 30 = 5818 nt of 5'UTR (on reverse strand)
    assert result1["utr5_length"] == 5818
    assert result1["utr3_length"] == 152

    # Example 2: + strand, CDS starts in exon 3
    result2 = calculate_utr_lengths(
        strand="+",
        ref_cds_info=[(3, 50), (4, 120), (5, 80)],
        transcript_exon_info=[('1', 200), ('2', 150), ('3', 100), ('4', 120), ('5', 80), ('6', 300)],
        ref_cds_start=100500,
        transcript_start=100000,
        ref_cds_stop=102000,
        transcript_end=103000
    )
    # Exons 1 & 2: full UTR = 200 + 150 = 350
    # Exon 3: 100 - 50 = 50 of 5'UTR
    # Exon 6: full UTR = 300
    total_utr5 = 200 + 150 + 50
    total_utr3 = 300
    assert result2["utr5_length"] == total_utr5
    assert result2["utr3_length"] == total_utr3

    # Example 3: single exon on plus strand, CDS fully inside it
    result3 = calculate_utr_lengths(
        strand="+",
        ref_cds_info=[(2, 60)],
        transcript_exon_info=[('2', 150)],
        ref_cds_start=5000,
        transcript_start=4950,
        ref_cds_stop=5060,
        transcript_end=5100
    )
    # UTR5: 5000 - 4950 = 50
    # UTR3: 5100 - 5060 = 40
    assert result3["utr5_length"] == 50
    assert result3["utr3_length"] == 40

    # Example 4: single exon minus strand, CDS fully inside it
    result = calculate_utr_lengths(
        strand="-",
        ref_cds_info=[(1, 60)],
        transcript_exon_info=[('1', 150)],
        ref_cds_start=5060,
        transcript_start=4950,
        ref_cds_stop=5000,
        transcript_end=5100
    )
    assert result["utr5_length"] == 100  # t_end - cds_end = 5100 - 5000
    assert result["utr3_length"] == 110  # cds_start - t_start = 5060 - 4950

    # Example 5: missing information
    # Example 5.1: missing ref_cds_info
    result = calculate_utr_lengths(
        strand="+",
        ref_cds_info=None,
        transcript_exon_info=[('1', 200), ('2', 300)],
        ref_cds_start=None,
        transcript_start=None,
        ref_cds_stop=None,
        transcript_end=None
    )
    assert result["utr5_length"] is None
    assert result["utr3_length"] is None
    # Example 5.2: missing transcript_exon_info
    result = calculate_utr_lengths(
        strand="-",
        ref_cds_info=[(1, 100), (2, 150)],
        transcript_exon_info=None,
        ref_cds_start=None,
        transcript_start=None,
        ref_cds_stop=None,
        transcript_end=None
    )
    assert result["utr5_length"] is None
    assert result["utr3_length"] is None

    # Example 6: not continuous exons, plus strand
    result = calculate_utr_lengths(
        strand="+",
        ref_cds_info=[(1, 60), (3, 80), (5, 100)],
        transcript_exon_info=[('1', 100), ('2', 100), ('3', 100), ('4', 100), ('5', 100)],
        ref_cds_start=5000,
        transcript_start=4800,
        ref_cds_stop=5340,
        transcript_end=5400
    )
    # CDS spans exons 1 (60 used), 3 (80 used), 5 (100 used)
    # 5'UTR = 40 (from exon 1) + 100 (exon 3)
    # 3'UTR = 100 (exon 5)
    assert result["utr5_length"] == 40
    assert result["utr3_length"] == 0

def test_calculate_exon_features():

    # Example from tcga
    result = calculate_exon_features(
        transcript_exon_info=[
            ('1', 250), ('2', 97), ('3', 174), ('4', 123),
            ('5', 70), ('6', 173), ('7', 105), ('8', 5848)
        ],
        alt_stop_codon_exons=[1, 1],
        alt_is_premature=True
    )
    assert result["total_exon_count"] == 8
    assert result["upstream_exon_count"] == 0
    assert result["downstream_exon_count"] == 7

    # Example 1: + strand, PTC in middle exon
    result1 = calculate_exon_features(
        transcript_exon_info=[('1', 100), ('2', 150), ('3', 120), ('4', 110)],
        alt_stop_codon_exons=[2],
        alt_is_premature=True
    )
    assert result1 == {
        "total_exon_count": 4,
        "upstream_exon_count": 1,
        "downstream_exon_count": 2
    }

    # Example 2: - strand, PTC in exon 2 (which is second in reverse)
    result2 = calculate_exon_features(
        transcript_exon_info=[('1', 100), ('2', 150), ('3', 120), ('4', 110)],
        alt_stop_codon_exons=[2],
        alt_is_premature=True
    )
    assert result2 == {
        "total_exon_count": 4,
        "upstream_exon_count": 1,   # reversed: [4,3,2,1] → index 1
        "downstream_exon_count": 2
    }

    # Example 3: PTC in first exon on + strand
    result3 = calculate_exon_features(
        transcript_exon_info=[('1', 100), ('2', 150), ('3', 120)],
        alt_stop_codon_exons=[1],
        alt_is_premature=True
    )
    assert result3 == {
        "total_exon_count": 3,
        "upstream_exon_count": 0,
        "downstream_exon_count": 2
    }

    # Example 4: PTC in last exon on - strand
    result4 = calculate_exon_features(
        transcript_exon_info=[('1', 100), ('2', 150), ('3', 120)],
        alt_stop_codon_exons=[1],
        alt_is_premature=True
    )
    assert result4 == {
        "total_exon_count": 3,
        "upstream_exon_count": 0,
        "downstream_exon_count": 2
    }

    # Example 5: Single exon transcript
    result5 = calculate_exon_features(
        transcript_exon_info=[('1', 500)],
        alt_stop_codon_exons=[1],
        alt_is_premature=True
    )
    assert result5 == {
        "total_exon_count": 1,
        "upstream_exon_count": 0,
        "downstream_exon_count": 0
    }

    # Example 7: PTC exon not present in transcript
    # Example 6: Invalid stop codon exon (doesn't exist)
    result6 = calculate_exon_features(
        transcript_exon_info=[('1', 100), ('2', 200)],
        alt_stop_codon_exons=[99],
        alt_is_premature=True
    )
    assert result6 == {
        "total_exon_count": 2,
        "upstream_exon_count": None,
        "downstream_exon_count": None
    }

    # Example 7: Missing stop codon exons
    result7 = calculate_exon_features(
        transcript_exon_info=[('1', 100), ('2', 200)],
        alt_stop_codon_exons=None,
        alt_is_premature=True
    )
    assert result7 == {
        "total_exon_count": 2,
        "upstream_exon_count": None,
        "downstream_exon_count": None
    }

    # Example 8: is not PTC
    result8 = calculate_exon_features(
        transcript_exon_info=[('1', 100), ('2', 200)],
        alt_stop_codon_exons=[99],
        alt_is_premature=False
    )
    assert result8 == {
        "total_exon_count": 2,
        "upstream_exon_count": None,
        "downstream_exon_count": None
    }

    # Example 9: Transcript with no exons
    result9 = calculate_exon_features(
        transcript_exon_info=[],
        alt_stop_codon_exons=[],
        alt_is_premature=False
    )
    assert result9 == {
        "total_exon_count": None,
        "upstream_exon_count": None,
        "downstream_exon_count": None
    }

def test_calculate_ptc_to_start_distance():
    assert calculate_ptc_to_start_distance(
        alt_is_premature=True,
        alt_start_codon_pos=50,
        alt_first_stop_pos=215
    ) == 165

    # PTC before start codon
    assert calculate_ptc_to_start_distance(
        alt_is_premature=True,
        alt_start_codon_pos=100,
        alt_first_stop_pos=30
    ) is None

    # Same position → distance 0 # can not happen
    assert calculate_ptc_to_start_distance(
        alt_is_premature=True,
        alt_start_codon_pos=120,
        alt_first_stop_pos=120
    ) is None

    # PTC not premature → None
    assert calculate_ptc_to_start_distance(
        alt_is_premature=False,
        alt_start_codon_pos=200,
        alt_first_stop_pos=300
    ) is None

    # Missing alt_first_stop_pos → None
    assert calculate_ptc_to_start_distance(
        alt_is_premature=True,
        alt_start_codon_pos=200,
        alt_first_stop_pos=None
    ) is None

    # Missing alt_start_codon_pos → None
    assert calculate_ptc_to_start_distance(
        alt_is_premature=True,
        alt_start_codon_pos=None,
        alt_first_stop_pos=200
    ) is None

    # All fields missing → None
    assert calculate_ptc_to_start_distance(
        alt_is_premature=None,
        alt_start_codon_pos=None,
        alt_first_stop_pos=None
    ) is None

def test_calculate_ptc_exon_length():
    assert calculate_ptc_exon_length(
        alt_is_premature=True,
        alt_stop_codon_exons=[1, 1],
        transcript_exon_info=[('1', 250), ('2', 97), ('3', 174), ('4', 123), ('5', 70), ('6', 173), ('7', 105), ('8', 5848)]
    ) == 250

    assert calculate_ptc_exon_length(
        alt_is_premature=True,
        alt_stop_codon_exons=[2, 3, 4],
        transcript_exon_info=[('1', 250), ('2', 97), ('3', 174), ('4', 123), ('5', 70), ('6', 173), ('7', 105), ('8', 5848)]
    ) == 97

    assert calculate_ptc_exon_length(
        alt_is_premature=False,
        alt_stop_codon_exons=[2, 3, 4],
        transcript_exon_info=[('1', 250), ('2', 97), ('3', 174), ('4', 123), ('5', 70), ('6', 173), ('7', 105), ('8', 5848)]
    ) is None

    assert calculate_ptc_exon_length(
        alt_is_premature=True,
        alt_stop_codon_exons=[],
        transcript_exon_info=[('1', 250), ('2', 97), ('3', 174), ('4', 123), ('5', 70), ('6', 173), ('7', 105), ('8', 5848)]
    ) is None

def test_calculate_stop_codon_dist():
    # Case 1: PTC upstream of reference stop
    assert calculate_stop_codon_dist(
        ref_first_stop_pos=1000,
        alt_first_stop_pos=800
    ) == 200

    # Case 2: PTC downstream of reference stop (rare, negative distance)
    assert calculate_stop_codon_dist(
        ref_first_stop_pos=800,
        alt_first_stop_pos=1000
    ) == -200

    # Case 3: PTC exactly at reference stop
    assert calculate_stop_codon_dist(
        ref_first_stop_pos=900,
        alt_first_stop_pos=900
    ) == 0

    # Case 4: Missing alt stop codon
    assert calculate_stop_codon_dist(
        ref_first_stop_pos=900,
        alt_first_stop_pos=None
    ) is None

    # Case 5: Missing ref stop codon
    assert calculate_stop_codon_dist(
        ref_first_stop_pos=None,
        alt_first_stop_pos=750
    ) is None

def test_evaluate_nmd_escape_rules():

    # Example 1: Single exon rule
    result = evaluate_nmd_escape_rules(
        alt_is_premature=True,
        alt_first_stop_pos=90,
        alt_stop_codon_exons=[1],
        transcript_exon_info=[(1, 100)],
        alt_start_codon_pos=0,
        total_exon_count=1,
        downstream_exon_count=0,
        ptc_exon_length=100
    )
    assert result["nmd_single_exon_rule"] == True
    assert result["nmd_last_exon_rule"] == True
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == True
    assert result["nmd_escape"] == True

    # Example 2: Last exon rule
    result = evaluate_nmd_escape_rules(
        alt_is_premature=True,
        alt_first_stop_pos=250,
        alt_stop_codon_exons=[3],
        transcript_exon_info=[(1, 100), (2, 100), (3, 100)],
        alt_start_codon_pos=0,
        total_exon_count=3,
        downstream_exon_count=0,
        ptc_exon_length=100
    )
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == True
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == True

    # Example 3: 50nt from penultimate exon end
    result = evaluate_nmd_escape_rules(
        alt_is_premature=True,
        alt_first_stop_pos=160,
        alt_stop_codon_exons=[2],
        transcript_exon_info=[(1, 100), (2, 100), (3, 100)],
        alt_start_codon_pos=0,
        total_exon_count=3,
        downstream_exon_count=1,
        ptc_exon_length=100
    )
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == True
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == True

    # Example 4: Long exon rule
    result = evaluate_nmd_escape_rules(
        alt_is_premature=True,
        alt_first_stop_pos=300,
        alt_stop_codon_exons=[2],
        transcript_exon_info=[(1, 100), (2, 500), (3, 100)],
        alt_start_codon_pos=0,
        total_exon_count=3,
        downstream_exon_count=1,
        ptc_exon_length=500
    )
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == True
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == True

    # Example 5: Start-proximal rule
    result = evaluate_nmd_escape_rules(
        alt_is_premature=True,
        alt_first_stop_pos=100,
        alt_stop_codon_exons=[2],
        transcript_exon_info=[(1, 100), (2, 100)],
        alt_start_codon_pos=0,
        total_exon_count=2,
        downstream_exon_count=0,
        ptc_exon_length=100
    )
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == True
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == True
    assert result["nmd_escape"] == True

    # Example 6: multiple escape rules
    result = evaluate_nmd_escape_rules(
        alt_is_premature=True,
        alt_first_stop_pos=145,
        alt_stop_codon_exons=[3],
        transcript_exon_info=[(1, 100), (2, 100), (3, 500)],
        alt_start_codon_pos=0,
        total_exon_count=3,
        downstream_exon_count=0,
        ptc_exon_length=500
    )
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == True
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == True
    assert result["nmd_start_proximal_rule"] == True
    assert result["nmd_escape"] == True

    # Example 7: Not premature → should skip
    result = evaluate_nmd_escape_rules(
        alt_is_premature=False,
        alt_first_stop_pos=150,
        alt_stop_codon_exons=[2],
        transcript_exon_info=[(1, 100), (2, 100)],
        alt_start_codon_pos=0,
        total_exon_count=2,
        downstream_exon_count=1,
        ptc_exon_length=100
    )

    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == False

def test_calculate_ptc_to_downstream_ej():
    # Case 1: PTC in exon 2, simple transcript
    # End of exon 2: 100 + 200 = 300, distance = 300 - 250 = 50
    assert calculate_ptc_to_downstream_ej(
        alt_is_premature=True,
        alt_stop_codon_exons=[2],
        alt_cds_info=[(1, 100), (2, 200), (3, 150)],
        alt_first_stop_pos=250
    ) == 50

    # Case 2: PTC in first exon
    # End of exon 1: 100, distance = 100 - 60 = 40
    assert calculate_ptc_to_downstream_ej(
        alt_is_premature=True,
        alt_stop_codon_exons=[1],
        alt_cds_info=[(1, 100), (2, 200), (3, 150)],
        alt_first_stop_pos=60
    ) == 40

    # Case 3: PTC in last exon
    # End of exon 3: 100+200+200=500, distance = 500 - 430 = 70
    assert calculate_ptc_to_downstream_ej(
        alt_is_premature=True,
        alt_stop_codon_exons=[3],
        alt_cds_info=[(1, 100), (2, 200), (3, 200)],
        alt_first_stop_pos=430
    ) == 70

    # Case 4: Multiple stop codons, take the smallest exon number
    # Smallest exon = 2, end of exon 2: 100+200=300, distance = 300 - 350 = -50 (PTC past exon end)
    assert calculate_ptc_to_downstream_ej(
        alt_is_premature=True,
        alt_stop_codon_exons=[2, 3],
        alt_cds_info=[(1, 100), (2, 200), (3, 200)],
        alt_first_stop_pos=350
    ) == -50

    # Case 5: Not premature → should return None
    assert calculate_ptc_to_downstream_ej(
        alt_is_premature=False,
        alt_stop_codon_exons=[2],
        alt_cds_info=[(1, 100), (2, 200)],
        alt_first_stop_pos=250
    ) is None

    # Case 6: Missing data → should return None
    assert calculate_ptc_to_downstream_ej(
        alt_is_premature=True,
        alt_stop_codon_exons=[2],
        alt_cds_info=[(1, 100), (2, 200)],
        alt_first_stop_pos=None
    ) is None

def test_add_likely_misannotated_flag():

    # Baseline: "good" annotation, not misannotated
    assert add_likely_misannotated_flag(
        cds_in_transcript=True,
        ref_start_codon_pos=0,
        ref_valid_stop=True
    ) is False

    # CDS not in transcript
    assert add_likely_misannotated_flag(
        cds_in_transcript=False,
        ref_start_codon_pos=0,
        ref_valid_stop=True
    ) is True

    # Start position of CDS not at position 0
    assert add_likely_misannotated_flag(
        cds_in_transcript=True,
        ref_start_codon_pos=14,  # start codon is not at beginning of CDS
        ref_valid_stop=True
    ) is True

    # No valid stop codon at the end of the CDS sequence
    assert add_likely_misannotated_flag(
        cds_in_transcript=True,
        ref_start_codon_pos=0,
        ref_valid_stop=False
    ) is True

    # Multiple conditions that point to a likely misannotation
    assert add_likely_misannotated_flag(
        cds_in_transcript=True,
        ref_start_codon_pos=13,
        ref_valid_stop=False
    ) is True

    # Missing information
    assert add_likely_misannotated_flag(
        cds_in_transcript=None,
        ref_start_codon_pos=None,
        ref_valid_stop=None
    ) is True

    # No start codon position
    assert add_likely_misannotated_flag(
        cds_in_transcript=True,
        ref_start_codon_pos=None,
        ref_valid_stop=True
    ) is True



# def test_calculate_ptc_to_downstream_ej_old():
#     # Case 1: PTC in exon 2, simple transcript
#     row1 = {
#         "alt_is_premature": True,
#         "alt_first_stop_pos": 250,  # PTC position
#         "alt_stop_codon_exons": [2],
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 200),
#             (3, 150)
#         ]
#     }
#     # End of exon 2: 100 + 200 = 300, distance = 300 - 250 = 50
#     assert calculate_ptc_to_downstream_ej(row1) == 50
#
#     # Case 2: PTC in first exon
#     row2 = {
#         "alt_is_premature": True,
#         "alt_first_stop_pos": 60,
#         "alt_stop_codon_exons": [1],
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 200),
#             (3, 150)
#         ]
#     }
#     # End of exon 1: 100, distance = 100 - 60 = 40
#     assert calculate_ptc_to_downstream_ej(row2) == 40
#
#     # Case 3: PTC in last exon
#     row3 = {
#         "alt_is_premature": True,
#         "alt_first_stop_pos": 430,
#         "alt_stop_codon_exons": [3],
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 200),
#             (3, 200)
#         ]
#     }
#     # End of exon 3: 100+200+200=500, distance = 500 - 430 = 70
#     assert calculate_ptc_to_downstream_ej(row3) == 70
#
#     # Case 4: Multiple stop codons, take the smallest exon number
#     row4 = {
#         "alt_is_premature": True,
#         "alt_first_stop_pos": 350,
#         "alt_stop_codon_exons": [2, 3],
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 200),
#             (3, 200)
#         ]
#     }
#     # Smallest exon = 2, end of exon 2: 100+200=300, distance = 300 - 350 = -50 (PTC past exon end)
#     assert calculate_ptc_to_downstream_ej(row4) == -50
#
#     # Case 5: Not premature → should return None
#     row5 = {
#         "alt_is_premature": False,
#         "alt_first_stop_pos": 250,
#         "alt_stop_codon_exons": [2],
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 200)
#         ]
#     }
#     assert calculate_ptc_to_downstream_ej(row5) is None
#
#     # Case 6: Missing data → should return None
#     row6 = {
#         "alt_is_premature": True,
#         "alt_first_stop_pos": None,
#         "alt_stop_codon_exons": [2],
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 200)
#         ]
#     }
#     assert calculate_ptc_to_downstream_ej(row6) is None

# def test_evaluate_nmd_escape_rules_old():
#
#     # Example 1: Single exon rule
#     row2 = {
#         "alt_is_premature": True,
#         "alt_first_stop_pos": 90,
#         "alt_stop_codon_exons": [1],
#         "alt_start_codon_pos": 0,
#         "transcript_exon_info": [
#             (1, 100)
#         ]
#     }
#
#     result = evaluate_nmd_escape_rules(row2)
#     assert result["nmd_single_exon_rule"] == True
#     assert result["nmd_last_exon_rule"] == False
#     assert result["nmd_50nt_penultimate_rule"] == False
#     assert result["nmd_long_exon_rule"] == False
#     assert result["nmd_start_proximal_rule"] == True
#     assert result["nmd_escape"] == True
#
#     # Example 2: Last exon rule
#     row = {
#         "alt_is_premature": True,
#         "alt_first_stop_pos": 250,
#         "alt_stop_codon_exons": [3],
#         "alt_start_codon_pos": 0,
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 100),
#             (3, 100)
#         ]
#     }
#
#     result = evaluate_nmd_escape_rules(row)
#     assert result["nmd_single_exon_rule"] == False
#     assert result["nmd_last_exon_rule"] == True
#     assert result["nmd_50nt_penultimate_rule"] == False
#     assert result["nmd_long_exon_rule"] == False
#     assert result["nmd_start_proximal_rule"] == False
#     assert result["nmd_escape"] == True
#
#     # Example 3: 50nt from penultimate exon end
#     row = {
#         "alt_is_premature": True,
#         "alt_first_stop_pos": 160,
#         "alt_stop_codon_exons": [2],
#         "alt_start_codon_pos": 0,
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 100),
#             (3, 100)
#         ]
#     }
#
#     result = evaluate_nmd_escape_rules(row)
#     assert result["nmd_single_exon_rule"] == False
#     assert result["nmd_last_exon_rule"] == False
#     assert result["nmd_50nt_penultimate_rule"] == True
#     assert result["nmd_long_exon_rule"] == False
#     assert result["nmd_start_proximal_rule"] == False
#     assert result["nmd_escape"] == True
#
#     # Example 4: Long exon rule
#     row = {
#         "alt_is_premature": True,
#         "alt_first_stop_pos": 300,
#         "alt_stop_codon_exons": [2],
#         "alt_start_codon_pos": 0,
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 500),
#             (3, 100)
#         ]
#     }
#
#     result = evaluate_nmd_escape_rules(row)
#     assert result["nmd_single_exon_rule"] == False
#     assert result["nmd_last_exon_rule"] == False
#     assert result["nmd_50nt_penultimate_rule"] == False
#     assert result["nmd_long_exon_rule"] == True
#     assert result["nmd_start_proximal_rule"] == False
#     assert result["nmd_escape"] == True
#
#     # Example 5: Start-proximal rule
#     row = {
#         "alt_is_premature": True,
#         "alt_first_stop_pos": 100,
#         "alt_stop_codon_exons": [2],
#         "alt_start_codon_pos": 0,
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 100)
#         ]
#     }
#
#     result = evaluate_nmd_escape_rules(row)
#     assert result["nmd_single_exon_rule"] == False
#     assert result["nmd_last_exon_rule"] == True
#     assert result["nmd_50nt_penultimate_rule"] == False
#     assert result["nmd_long_exon_rule"] == False
#     assert result["nmd_start_proximal_rule"] == True
#     assert result["nmd_escape"] == True
#
#     # Example 6: multiple escape rules
#     row = {
#         "alt_is_premature": True,
#         "alt_first_stop_pos": 145,
#         "alt_stop_codon_exons": [3],
#         "alt_start_codon_pos": 0,
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 100),
#             (3, 500)
#         ]
#     }
#
#     result = evaluate_nmd_escape_rules(row)
#     assert result["nmd_single_exon_rule"] == False
#     assert result["nmd_last_exon_rule"] == True
#     assert result["nmd_50nt_penultimate_rule"] == False
#     assert result["nmd_long_exon_rule"] == True
#     assert result["nmd_start_proximal_rule"] == True
#     assert result["nmd_escape"] == True
#
#     # Example 7: Not premature → should skip
#     row3 = {
#         "alt_is_premature": False,
#         "alt_first_stop_pos": 150,
#         "alt_stop_codon_exons": [2],
#         "alt_start_codon_pos": 0,
#         "transcript_exon_info": [
#             (1, 100),
#             (2, 100)
#         ]
#     }
#
#     result = evaluate_nmd_escape_rules(row3)
#     assert result["nmd_single_exon_rule"] == False
#     assert result["nmd_last_exon_rule"] == False
#     assert result["nmd_50nt_penultimate_rule"] == False
#     assert result["nmd_long_exon_rule"] == False
#     assert result["nmd_start_proximal_rule"] == False
#     assert result["nmd_escape"] == False