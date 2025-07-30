# Import dependencies
from nmd_scanner.rules import *

def test_adjust_last_cds_for_stop_codon():
    df = pd.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx2", "tx2", "tx2"],
        "exon_number": [1, 2, 1, 2, 3],
        "Start": [100, 200, 500, 800, 900],
        "End": [150, 250, 550, 850, 950],
        "Strand": ["+", "+", "-", "-", "-"]
    })

    adjusted = adjust_last_cds_for_stop_codon(df)

    # Get the last exon of tx1 (+ strand)
    tx1_last = adjusted[(adjusted.transcript_id == "tx1") & (adjusted.exon_number == 2)].iloc[0]
    assert tx1_last["Start"] == 200
    assert tx1_last["End"] == 253  # 250 + 3

    # Get the last exon of tx2 (- strand)
    tx2_last = adjusted[(adjusted.transcript_id == "tx2") & (adjusted.exon_number == 3)].iloc[0]
    assert tx2_last["Start"] == 897  # 900 - 3
    assert tx2_last["End"] == 950

def test_apply_variant_edge_aware_with_lengths():
    # need to keep in mind all the cases (variant goes over start or end of exon, indels, SNPs)
    # Maybe can use the test input and output files
        # input = resources/test_files/variants.vcf
        # output = resources/test_output_files/variant_exon_output.tsv

    # Load the test output file with expected results: I cross checked these for correctness
    df_expected = pd.read_csv("resources/test_output_files/variant_exon_output.tsv", sep="\t")

    # Apply the function to each row to get actual results
    df_actual = df_expected.copy()
    actual_cols = df_actual.apply(apply_variant_edge_aware_with_lengths, axis=1)

    # Attach the new columns to compare
    df_actual["Exon_Alt_CDS_seq_actual"] = actual_cols["Exon_Alt_CDS_seq"]
    df_actual["Exon_Alt_CDS_length_actual"] = actual_cols["Exon_Alt_CDS_length"]

    # Run assertions row-by-row to catch mismatches
    for i, row in df_actual.iterrows():
        assert row["Exon_Alt_CDS_seq"] == row["Exon_Alt_CDS_seq_actual"], f"Mismatch in Alt_CDS_seq at row {i}"
        assert row["Exon_Alt_CDS_length"] == row[
            "Exon_Alt_CDS_length_actual"], f"Mismatch in Alt_CDS_length at row {i}"

def test_create_reference_cds_using_file():
    # Load expected output
    expected = pd.read_csv("resources/test_output_files/create_reference_CDS.tsv", sep="\t")

    # Load df3 and cds_df_test from the previous step of your pipeline
    df3 = pd.read_csv("resources/test_output_files/variant_exon_output.tsv", sep="\t")
    cds_df_test = pd.read_csv("resources/test_output_files/cds_df_test.tsv",sep="\t")

    # Run the function
    actual = create_reference_cds(df3, cds_df_test)

    # Sort both for consistent comparison
    expected_sorted = expected.sort_values(["transcript_id", "variant_id"]).reset_index(drop=True)
    actual_sorted = actual.sort_values(["transcript_id", "variant_id"]).reset_index(drop=True)

    # Compare critical columns
    columns_to_check = [
        "transcript_id", "variant_id",
        "ref_cds_seq", "alt_cds_seq",
        "ref_cds_len", "alt_cds_len",
        "strand", "ref", "alt"
    ]

    for col in columns_to_check:
        assert all(expected_sorted[col] == actual_sorted[col]), f"Mismatch in column: {col}"

def test_create_reference_cds():
    # create example dataframe
    cds_df_test = pd.DataFrame({
        "transcript_id": ["tx1"] * 5,
        "exon_number": [2, 3, 4, 5, 6],
        "Chromosome": ["chr1"] * 5,
        "Start": [100, 200, 300, 400, 500],
        "End": [150, 250, 350, 450, 550],
        "Strand": ["+" for _ in range(5)],
        "Exon_CDS_seq": ["AAA", "CCC", "GGG", "TTT", "AAA"]
    })

    # Variant 1: SNP on exon 3
    variant_snp = {
        "transcript_id": "tx1",
        "exon_number": 3,
        "Chromosome": "chr1",
        "Start": 200,
        "End": 250,
        "Strand": "+",
        "ID": "var_snp",
        "Start_variant": 210,
        "End_variant": 211,
        "Ref": "C",
        "Alt": "T",
        "Exon_Alt_CDS_seq": "CCT"  # Same length
    }

    # Variant 2: Insertion on exon 4
    variant_ins = {
        "transcript_id": "tx1",
        "exon_number": 4,
        "Chromosome": "chr1",
        "Start": 300,
        "End": 350,
        "Strand": "+",
        "ID": "var_ins",
        "Start_variant": 325,
        "End_variant": 325,
        "Ref": "-",
        "Alt": "A",
        "Exon_Alt_CDS_seq": "GGGA"  # Inserted A
    }

    # Variant 3: Deletion on exon 5
    variant_del = {
        "transcript_id": "tx1",
        "exon_number": 5,
        "Chromosome": "chr1",
        "Start": 400,
        "End": 450,
        "Strand": "+",
        "ID": "var_del",
        "Start_variant": 440,
        "End_variant": 441,
        "Ref": "T",
        "Alt": "-",
        "Exon_Alt_CDS_seq": "TT"  # One T removed
    }

    # Variant 4: Deletion spanning exon 3 to 4
    spanning_del_3_4 = [
        {
            "transcript_id": "tx1",
            "exon_number": 3,
            "Chromosome": "chr1",
            "Start": 200,
            "End": 250,
            "Strand": "+",
            "ID": "var_spanning",
            "Start_variant": 240,
            "End_variant": 310,
            "Ref": "CCGG",
            "Alt": "-",
            "Exon_Alt_CDS_seq": "C"  # Shortened version
        },
        {
            "transcript_id": "tx1",
            "exon_number": 4,
            "Chromosome": "chr1",
            "Start": 300,
            "End": 350,
            "Strand": "+",
            "ID": "var_spanning",
            "Start_variant": 240,
            "End_variant": 310,
            "Ref": "CCGG",
            "Alt": "-",
            "Exon_Alt_CDS_seq": "G"  # Shortened version
        }
    ]

    # Combine all variants
    df3 = pd.DataFrame([variant_snp, variant_ins, variant_del] + spanning_del_3_4)

    # Run the function
    result = create_reference_cds(df3, cds_df_test)

    # Reference CDS sequence
    ref_seq = "AAACCCGGGTTTAAA"
    ref_len = len(ref_seq)

    for _, row in result.iterrows():
        assert row["ref_cds_seq"] == ref_seq
        assert row["ref_cds_len"] == ref_len

    alt_seqs = {row["variant_id"]: row["alt_cds_seq"] for _, row in result.iterrows()}
    alt_lens = {row["variant_id"]: row["alt_cds_len"] for _, row in result.iterrows()}

    # Variant-specific checks
    assert alt_seqs["var_snp"] == "AAACCTGGGTTTAAA"
    assert alt_lens["var_snp"] == ref_len

    assert alt_seqs["var_ins"] == "AAACCCGGGATTTAAA"
    assert alt_lens["var_ins"] == ref_len + 1

    assert alt_seqs["var_del"] == "AAACCCGGGTTAAA"
    assert alt_lens["var_del"] == ref_len - 1

    assert alt_seqs["var_spanning"] == "AAACGTTTAAA"
    assert alt_lens["var_spanning"] == ref_len - 4

def test_get_transcript_sequence():
    fasta = {"chr1": "AAAAAAAAAACCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"} # ("A" * 10 + "C" * 20 + "G" * 30 + "T" * 40)

    exons_df = pd.DataFrame([
        {"transcript_id": "tx1", "exon_number": 1, "Chromosome": "chr1", "Start": 10, "End": 13, "Strand": "+"},
        {"transcript_id": "tx1", "exon_number": 2, "Chromosome": "chr1", "Start": 20, "End": 24, "Strand": "+"},
        {"transcript_id": "tx1", "exon_number": 3, "Chromosome": "chr1", "Start": 30, "End": 35, "Strand": "+"},

        {"transcript_id": "tx2", "exon_number": 1, "Chromosome": "chr1", "Start": 40, "End": 43, "Strand": "-"},
        {"transcript_id": "tx2", "exon_number": 2, "Chromosome": "chr1", "Start": 50, "End": 53, "Strand": "-"},
        {"transcript_id": "tx2", "exon_number": 3, "Chromosome": "chr1", "Start": 60, "End": 63, "Strand": "-"},
    ])

    # Run function
    transcript_df = get_transcript_sequence(exons_df, fasta)

    # Check tx1
    tx1 = transcript_df[transcript_df["transcript_id"] == "tx1"].iloc[0]
    expected_tx1_seq = "CCCCCCCGGGGG"
    assert tx1["transcript_sequence"] == expected_tx1_seq
    assert tx1["transcript_length"] == len(expected_tx1_seq)
    assert tx1["transcript_exon_info"] == [(1, 3), (2, 4), (3, 5)]

    # Check tx2
    tx2 = transcript_df[transcript_df["transcript_id"] == "tx2"].iloc[0]
    expected_tx2_seq = "AAACCCCCC" # reverse complement of "GGGGGGTTT"
    assert tx2["transcript_sequence"] == expected_tx2_seq
    assert tx2["transcript_length"] == len(expected_tx2_seq)
    assert tx2["transcript_exon_info"] == [(3, 3), (2, 3), (1, 3)]  # reversed for minus strand

def test_get_exon():
    exon_info = [(1, 10), (2, 20), (3, 30)]
    assert get_exon(5, exon_info) == 1
    assert get_exon(25, exon_info) == 2
    assert get_exon(55, exon_info) == 3

def test_analyze_sequence():
    df = pd.DataFrame([{
        "ref_cds_seq": "ATGAAATAG",
        "alt_cds_seq": "ATGAAATAA",
        "ref_cds_info": [(1, 9)],
        "alt_cds_info": [(1, 9)],
    }])
    analyzed = analyze_sequence(df)
    assert analyzed.loc[0, "ref_start_codon_pos"] == 0
    assert analyzed.loc[0, "ref_valid_stop"] == True
    assert analyzed.loc[0, "alt_valid_stop"] == True

def test_start_stop_loss():
    df = pd.DataFrame([{
        "ref_start_codon_pos": 0,
        "alt_start_codon_pos": None,
        "ref_valid_stop": True,
        "alt_valid_stop": False,
        "ref_last_codon": "TAG",
        "alt_last_codon": "GGA"
    }])
    result = start_stop_loss(df)
    assert result["start_loss"].iloc[0] == True
    assert result["stop_loss"].iloc[0] == True

def test_splice_alt_cds_into_transcript():
    row = {
        "ref_cds_seq": "AAAGGGCCC",
        "alt_cds_seq": "AAATTTCCC"
    }
    transcript_seq = "TTTAAAGGGCCCGGG"

    result = splice_alt_cds_into_transcript(row, transcript_seq)
    assert result == "TTTAAATTTCCCGGG"

def test_analyze_transcript():

    df = pd.DataFrame([
        {
            "alt_transcript_seq": "CCCATGAAATAATAGGGG",  # ATG at pos 3, TAA at 9, TAG at 12
            "alt_cds_start": 0,
            "transcript_start": 0,
            "transcript_exon_info": [
                (1, 10, 0, 10),
                (2, 10, 10, 20)
            ],
            "start_loss": True,
            "stop_loss": False
        }
    ])

    result = analyze_transcript(df)

    # Check values
    row = result.loc[0]
    assert row["transcript_start_codon_pos"] == 3
    assert row["transcript_start_codon_exon"] == 1
    assert row["transcript_first_stop_codon"] == "TAA"
    assert row["transcript_first_stop_pos"] == 9
    assert row["transcript_last_codon"] == "GGG"
    assert row["transcript_valid_stop"] is False
    assert row["transcript_num_stop_codons"] == 2
    assert row["transcript_all_stop_codons"] == [(9, "TAA"), (12, "TAG")]
    assert row["transcript_stop_codon_exons"] == [1, 2]

def test_evaluate_nmd_escape_rules():

    # Example 1: Single exon rule
    row2 = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 90,
        "alt_stop_codon_exons": [1],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [
            (1, 100)
        ]
    }

    result = evaluate_nmd_escape_rules(row2)
    assert result["nmd_single_exon_rule"] == True
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == True
    assert result["nmd_escape"] == True

    # Example 2: Last exon rule
    row = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 250,
        "alt_stop_codon_exons": [3],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [
            (1, 100),
            (2, 100),
            (3, 100)
        ]
    }

    result = evaluate_nmd_escape_rules(row)
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == True
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == True

    # Example 3: 50nt from penultimate exon end
    row = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 160,
        "alt_stop_codon_exons": [2],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [
            (1, 100),
            (2, 100),
            (3, 100)
        ]
    }

    result = evaluate_nmd_escape_rules(row)
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == True
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == True

    # Example 4: Long exon rule
    row = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 300,
        "alt_stop_codon_exons": [2],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [
            (1, 100),
            (2, 500),
            (3, 100)
        ]
    }

    result = evaluate_nmd_escape_rules(row)
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == True
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == True

    # Example 5: Start-proximal rule
    row = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 100,
        "alt_stop_codon_exons": [2],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [
            (1, 100),
            (2, 100)
        ]
    }

    result = evaluate_nmd_escape_rules(row)
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == True
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == True
    assert result["nmd_escape"] == True

    # Example 6: multiple escape rules
    row = {
        "alt_is_premature": True,
        "alt_first_stop_pos": 145,
        "alt_stop_codon_exons": [3],
        "alt_start_codon_pos": 0,
        "transcript_exon_info": [
            (1, 100),
            (2, 100),
            (3, 500)
        ]
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
        "transcript_exon_info": [
            (1, 100),
            (2, 100)
        ]
    }

    result = evaluate_nmd_escape_rules(row3)
    assert result["nmd_single_exon_rule"] == False
    assert result["nmd_last_exon_rule"] == False
    assert result["nmd_50nt_penultimate_rule"] == False
    assert result["nmd_long_exon_rule"] == False
    assert result["nmd_start_proximal_rule"] == False
    assert result["nmd_escape"] == False