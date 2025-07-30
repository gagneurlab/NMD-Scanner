from nmd_scanner.extra_features import *

def test_calculate_utr_lengths():
    # Example 1: - strand, CDS spans exon 8 to 1
    row1 = {
        "strand": "-",
        "ref_cds_info": [(8, 30), (7, 105), (6, 173), (5, 70), (4, 123), (3, 174), (2, 97), (1, 98)],
        "transcript_exon_info": [('1', 250), ('2', 97), ('3', 174), ('4', 123), ('5', 70),
                                 ('6', 173), ('7', 105), ('8', 5848)],
        "ref_cds_start": 70925000,
        "ref_cds_stop": 70929273,
        "transcript_start": 70920000,
        "transcript_end": 70950000,
    }
    result1 = calculate_utr_lengths(row1)
    # Exon 1 has 250 - 98 = 152 nt of 3'UTR
    # Exon 8 has 5848 - 30 = 5818 nt of 5'UTR (on reverse strand)
    assert result1["utr5_length"] == 5818
    assert result1["utr3_length"] == 152

    # Example 2: + strand, CDS starts in exon 3
    row2 = {
        "strand": "+",
        "ref_cds_info": [(3, 50), (4, 120), (5, 80)],
        "transcript_exon_info": [('1', 200), ('2', 150), ('3', 100), ('4', 120), ('5', 80), ('6', 300)],
        "ref_cds_start": 100500,
        "ref_cds_stop": 102000,
        "transcript_start": 100000,
        "transcript_end": 103000,
    }
    result2 = calculate_utr_lengths(row2)
    # Exons 1 & 2: full UTR = 200 + 150 = 350
    # Exon 3: 100 - 50 = 50 of 5'UTR
    # Exon 6: full UTR = 300
    total_utr5 = 200 + 150 + 50
    total_utr3 = 300
    assert result2["utr5_length"] == total_utr5
    assert result2["utr3_length"] == total_utr3

    # Example 3: single exon on plus strand, CDS fully inside it
    row3 = {
        "strand": "+",
        "ref_cds_info": [(2, 60)],
        "transcript_exon_info": [('2', 150)],
        "ref_cds_start": 5000,
        "ref_cds_stop": 5060,
        "transcript_start": 4950,
        "transcript_end": 5100,
    }
    result3 = calculate_utr_lengths(row3)
    # UTR5: 5000 - 4950 = 50
    # UTR3: 5100 - 5060 = 40
    assert result3["utr5_length"] == 50
    assert result3["utr3_length"] == 40

    # Example 4: single exon minus strand, CDS fully inside it
    row = {
        "strand": "-",
        "ref_cds_info": [(1, 60)],
        "transcript_exon_info": [('1', 150)],
        "ref_cds_start": 5060,
        "ref_cds_stop": 5000,
        "transcript_start": 4950,
        "transcript_end": 5100,
    }
    result = calculate_utr_lengths(row)
    assert result["utr5_length"] == 100  # t_end - cds_end = 5100 - 5000
    assert result["utr3_length"] == 110  # cds_start - t_start = 5060 - 4950

    # Example 5: missing information
    # Example 5.1: missing ref_cds_info
    row = {
        "strand": "+",
        "transcript_exon_info": [('1', 200), ('2', 300)],
    }
    result = calculate_utr_lengths(row)
    assert result["utr5_length"] is None
    assert result["utr3_length"] is None
    # Example 5.2: missing transcript_exon_info
    row = {
        "strand": "-",
        "ref_cds_info": [(1, 100), (2, 150)],
    }
    result = calculate_utr_lengths(row)
    assert result["utr5_length"] is None
    assert result["utr3_length"] is None

    # Example 6: not continuous exons, plus strand
    row = {
        "strand": "+",
        "ref_cds_info": [(1, 60), (3, 80), (5, 100)],
        "transcript_exon_info": [('1', 100), ('2', 100), ('3', 100), ('4', 100), ('5', 100)],
        "ref_cds_start": 5000,
        "ref_cds_stop": 5340,
        "transcript_start": 4800,
        "transcript_end": 5400,
    }
    result = calculate_utr_lengths(row)
    # CDS spans exons 1 (60 used), 3 (80 used), 5 (100 used)
    # 5'UTR = 40 (from exon 1) + 100 (exon 2)
    # 3'UTR = 100 (exon 4)
    assert result["utr5_length"] == 40
    assert result["utr3_length"] == 0

def test_calculate_exon_features():

    # Example from tcga
    row = {
        "alt_is_premature" : True,
        "strand": "+",
        "alt_stop_codon_exons": [1, 1],
        "transcript_exon_info": [
            ('1', 250), ('2', 97), ('3', 174), ('4', 123),
            ('5', 70), ('6', 173), ('7', 105), ('8', 5848)
        ]
    }

    result = calculate_exon_features(row)
    assert result["total_exon_count"] == 8
    assert result["upstream_exon_count"] == 0
    assert result["downstream_exon_count"] == 7

    # Example 1: + strand, PTC in middle exon
    row1 = {
        "alt_is_premature": True,
        "strand": "+",
        "transcript_exon_info": [('1', 100), ('2', 150), ('3', 120), ('4', 110)],
        "alt_stop_codon_exons": [2]
    }
    result1 = calculate_exon_features(row1)
    assert result1 == {
        "total_exon_count": 4,
        "upstream_exon_count": 1,
        "downstream_exon_count": 2
    }

    # Example 2: - strand, PTC in exon 2 (which is second in reverse)
    row2 = {
        "alt_is_premature": True,
        "strand": "-",
        "transcript_exon_info": [('1', 100), ('2', 150), ('3', 120), ('4', 110)],
        "alt_stop_codon_exons": [2]
    }
    result2 = calculate_exon_features(row2)
    assert result2 == {
        "total_exon_count": 4,
        "upstream_exon_count": 2,   # reversed: [4,3,2,1] → index 1
        "downstream_exon_count": 1
    }

    # Example 3: PTC in first exon on + strand
    row3 = {
        "alt_is_premature": True,
        "strand": "+",
        "transcript_exon_info": [('1', 100), ('2', 150), ('3', 120)],
        "alt_stop_codon_exons": [1]
    }
    result3 = calculate_exon_features(row3)
    assert result3 == {
        "total_exon_count": 3,
        "upstream_exon_count": 0,
        "downstream_exon_count": 2
    }

    # Example 4: PTC in last exon on - strand
    row4 = {
        "alt_is_premature": True,
        "strand": "-",
        "transcript_exon_info": [('1', 100), ('2', 150), ('3', 120)],
        "alt_stop_codon_exons": [1]
    }
    result4 = calculate_exon_features(row4)
    assert result4 == {
        "total_exon_count": 3,
        "upstream_exon_count": 2,
        "downstream_exon_count": 0
    }

    # Example 5: Single exon transcript
    row5 = {
        "alt_is_premature": True,
        "strand": "+",
        "transcript_exon_info": [('1', 500)],
        "alt_stop_codon_exons": [1]
    }
    result5 = calculate_exon_features(row5)
    assert result5 == {
        "total_exon_count": 1,
        "upstream_exon_count": 0,
        "downstream_exon_count": 0
    }

    # Example 7: PTC exon not present in transcript
    row6 = {
        "alt_is_premature": True,
        "strand": "-",
        "transcript_exon_info": [('1', 100), ('2', 200)],
        "alt_stop_codon_exons": [99]
    }
    result6 = calculate_exon_features(row6)
    assert result6 == {
        "total_exon_count": 2,
        "upstream_exon_count": None,
        "downstream_exon_count": None
    }

    # Example 7: Missing stop codon exons
    row7 = {
        "alt_is_premature": True,
        "strand": "+",
        "transcript_exon_info": [('1', 100), ('2', 200)],
        "alt_stop_codon_exons": []
    }
    result7 = calculate_exon_features(row7)
    assert result7 == {
        "total_exon_count": None,
        "upstream_exon_count": None,
        "downstream_exon_count": None
    }

    # Example 8: is not PTC
    row8 = {
        "alt_is_premature": False,
        "strand": "-",
        "transcript_exon_info": [('1', 100), ('2', 200)],
        "alt_stop_codon_exons": [99]
    }
    result8 = calculate_exon_features(row8)
    assert result8 == {
        "total_exon_count": None,
        "upstream_exon_count": None,
        "downstream_exon_count": None
    }

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

    # Same position → distance 0
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
        "alt_stop_codon_exons" : [1, 1],
        "transcript_exon_info": [('1', 250), ('2', 97), ('3', 174), ('4', 123), ('5', 70), ('6', 173), ('7', 105), ('8', 5848)],
    }
    assert calculate_ptc_exon_length(row) == 250

    row2 = {
        "alt_is_premature": True,
        "alt_stop_codon_exons": [2, 3, 4],
        "transcript_exon_info": [('1', 250), ('2', 97), ('3', 174), ('4', 123), ('5', 70), ('6', 173), ('7', 105),
                                 ('8', 5848)],
    }
    assert calculate_ptc_exon_length(row2) == 97

    row3 = {
        "alt_is_premature": False,
        "alt_stop_codon_exons": [2, 3, 4],
        "transcript_exon_info": [('1', 250), ('2', 97), ('3', 174), ('4', 123), ('5', 70), ('6', 173), ('7', 105),
                                 ('8', 5848)],
    }
    assert calculate_ptc_exon_length(row3) is None

    row4 = {
        "alt_is_premature": True,
        "alt_stop_codon_exons": [],
        "transcript_exon_info": [('1', 250), ('2', 97), ('3', 174), ('4', 123), ('5', 70), ('6', 173), ('7', 105),
                                 ('8', 5848)],
    }
    assert calculate_ptc_exon_length(row4) is None