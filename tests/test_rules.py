import pytest
import pandas as pd
import pysam
import csv
from nmd_scanner.rules import *
from pathlib import Path

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

    print("All variant edge-aware tests passed using the precomputed input/output files.")

def test_create_reference_cds():
    # TODO
    assert None

def test_get_transcript_sequence():
    # TODO
    assert None

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
    assert result["start_loss"].iloc[0] is True
    assert result["stop_loss"].iloc[0] is True

def test_splice_alt_cds_into_transcript():
    row = {
        "ref_cds_seq": "AAAGGGCCC",
        "alt_cds_seq": "AAATTTCCC"
    }
    transcript_seq = "TTTAAAGGGCCCGGG"

    result = splice_alt_cds_into_transcript(row, transcript_seq)
    assert result == "TTTAAATTTCCCGGG"

def test_analyze_transcript():
    # TODO
    assert None

def test_evaluate_nmd_escape_rules():
    # TODO
    assert None