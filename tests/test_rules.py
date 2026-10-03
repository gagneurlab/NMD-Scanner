# Import dependencies
import pandas as pd
import pyranges as pr
import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

from nmd_scanner.rules import (
    analyze_sequence,
    analyze_transcript,
    apply_variant_edge_aware_with_lengths,
    create_reference_cds,
    extract_ptc,
    get_exon,
    get_transcript_sequence,
    merge_stop_codons_into_cds,
    splice_alt_cds_into_transcript,
    start_stop_loss,
)
from nmd_scanner.scan import read_gtf


def _coding_rows(rows):
    """CDS and stop_codon rows of one transcript: (Feature, exon_number, Start, End, Strand)."""
    return pd.DataFrame(
        [
            {"transcript_id": "tx", "Feature": f, "exon_number": str(e), "Start": s, "End": en, "Strand": st}
            for f, e, s, en, st in rows
        ]
    )


def _intervals(df):
    return sorted(zip(df["exon_number"], df["Start"], df["End"], df["Feature"]))


def test_merge_stop_codons_into_cds_gtf():
    # GTF: the stop codon follows the last CDS on + strand, precedes it on - strand
    plus = _coding_rows([("CDS", 1, 100, 150, "+"), ("CDS", 2, 200, 250, "+"), ("stop_codon", 2, 250, 253, "+")])
    assert _intervals(merge_stop_codons_into_cds(plus)) == [(1, 100, 150, "CDS"), (2, 200, 253, "CDS")]

    minus = _coding_rows([("CDS", 2, 500, 550, "-"), ("CDS", 1, 800, 850, "-"), ("stop_codon", 2, 497, 500, "-")])
    assert _intervals(merge_stop_codons_into_cds(minus)) == [(1, 800, 850, "CDS"), (2, 497, 550, "CDS")]


def test_merge_stop_codons_into_cds_without_stop_codon():
    # cds_end_NF: no stop_codon row, so the CDS stays as it is
    rows = _coding_rows([("CDS", 1, 100, 150, "+"), ("CDS", 2, 200, 250, "+")])
    assert _intervals(merge_stop_codons_into_cds(rows)) == [(1, 100, 150, "CDS"), (2, 200, 250, "CDS")]


def test_merge_stop_codons_into_cds_split_stop_codon():
    # stop codon split across an intron: 2 bases at the end of exon 2, 1 base at the start of exon 3,
    # which has no CDS row
    rows = _coding_rows(
        [
            ("CDS", 1, 100, 150, "+"),
            ("CDS", 2, 200, 248, "+"),
            ("stop_codon", 2, 248, 250, "+"),
            ("stop_codon", 3, 300, 301, "+"),
        ]
    )
    assert _intervals(merge_stop_codons_into_cds(rows)) == [
        (1, 100, 150, "CDS"),
        (2, 200, 250, "CDS"),
        (3, 300, 301, "CDS"),
    ]


def test_merge_stop_codons_into_cds_split_stop_codon_minus_strand():
    # stop codon split across an intron on the - strand: 2 bases at the lower end of exon 2, 1 base in exon 3,
    # which lies at lower coordinates and has no CDS row
    rows = _coding_rows(
        [
            ("CDS", 1, 800, 850, "-"),
            ("CDS", 2, 500, 550, "-"),
            ("stop_codon", 2, 498, 500, "-"),
            ("stop_codon", 3, 400, 401, "-"),
        ]
    )
    assert _intervals(merge_stop_codons_into_cds(rows)) == [
        (1, 800, 850, "CDS"),
        (2, 498, 550, "CDS"),
        (3, 400, 401, "CDS"),
    ]


def test_merge_stop_codons_into_cds_gff3():
    # GFF3: the CDS already includes the stop codon, so the union changes nothing
    rows = _coding_rows([("CDS", 1, 100, 150, "+"), ("CDS", 2, 200, 253, "+"), ("stop_codon", 2, 250, 253, "+")])
    assert _intervals(merge_stop_codons_into_cds(rows)) == [(1, 100, 150, "CDS"), (2, 200, 253, "CDS")]


def test_merge_stop_codons_into_cds_has_stop_codon():
    # the flag is per transcript and comes from the stop_codon rows
    with_stop = _coding_rows([("CDS", 1, 100, 150, "+"), ("CDS", 2, 200, 250, "+"), ("stop_codon", 2, 250, 253, "+")])
    without_stop = _coding_rows([("CDS", 1, 100, 150, "+"), ("CDS", 2, 200, 250, "+")]).assign(transcript_id="tx_nf")
    merged = merge_stop_codons_into_cds(pd.concat([with_stop, without_stop], ignore_index=True))
    assert merged.groupby("transcript_id")["has_stop_codon"].agg(set).to_dict() == {"tx": {True}, "tx_nf": {False}}


def test_merge_stop_codons_into_cds_warns_without_stop_codon_rows(caplog):
    rows = _coding_rows([("CDS", 1, 100, 150, "+"), ("CDS", 2, 200, 250, "+")])
    with caplog.at_level("WARNING", logger="nmd_scanner.rules"):
        merge_stop_codons_into_cds(rows)
    assert "No stop_codon rows" in caplog.text
    assert "stop_codon rows together with the CDS rows" in caplog.text

    caplog.clear()
    with_stop = _coding_rows([("CDS", 1, 100, 150, "+"), ("stop_codon", 1, 150, 153, "+")])
    with caplog.at_level("WARNING", logger="nmd_scanner.rules"):
        merge_stop_codons_into_cds(with_stop)
    assert "No stop_codon rows" not in caplog.text


def test_merge_stop_codons_into_cds_rejects_gap():
    rows = _coding_rows([("CDS", 1, 100, 150, "+"), ("stop_codon", 1, 160, 163, "+")])
    with pytest.raises(ValueError, match="tx"):
        merge_stop_codons_into_cds(rows)


def _coding_sequence(coding, fasta):
    coding = coding.sort_values("Start")
    seq = "".join(fasta[c][s:e].seq.upper() for c, s, e in zip(coding["Chromosome"], coding["Start"], coding["End"]))
    return str(Seq(seq).reverse_complement()) if coding["Strand"].iloc[0] == "-" else seq


def test_merge_stop_codons_into_cds_on_real_transcripts():
    gtf_df = read_gtf("resources/chr18.gtf.gz").df
    fasta = Fasta("resources/chr18.fa.gz")
    rows = gtf_df[gtf_df["Feature"].isin(["CDS", "stop_codon"])]

    # ENST00000399496.8: stop codon split across an intron (two stop_codon rows)
    split = rows[rows["transcript_id"] == "ENST00000399496.8"]
    assert (split["Feature"] == "stop_codon").sum() == 2
    seq = _coding_sequence(merge_stop_codons_into_cds(split), fasta)
    assert len(seq) % 3 == 0
    assert seq[-3:] in {"TAA", "TAG", "TGA"}

    # ENST00000454642.3: minus strand, stop codon split across an intron (two stop_codon rows)
    split_minus = rows[rows["transcript_id"] == "ENST00000454642.3"]
    assert (split_minus["Feature"] == "stop_codon").sum() == 2
    assert split_minus["Strand"].iloc[0] == "-"
    seq = _coding_sequence(merge_stop_codons_into_cds(split_minus), fasta)
    assert len(seq) % 3 == 0
    assert seq[-3:] in {"TAA", "TAG", "TGA"}

    # a cds_end_NF transcript has no stop codon: its CDS rows stay as they are
    cds_end_nf = gtf_df.loc[(gtf_df["Feature"] == "transcript") & gtf_df["tag"].str.contains("cds_end_NF", na=False)]
    tx = cds_end_nf["transcript_id"].iloc[0]
    cds = rows[rows["transcript_id"] == tx]
    assert (cds["Feature"] == "stop_codon").sum() == 0
    merged = merge_stop_codons_into_cds(cds)
    assert sorted(zip(merged["Start"], merged["End"])) == sorted(zip(cds["Start"], cds["End"]))


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
        assert row["Exon_Alt_CDS_length"] == row["Exon_Alt_CDS_length_actual"], f"Mismatch in Alt_CDS_length at row {i}"


def test_apply_variant_edge_aware_with_lengths_with_DELs():
    cases = [
        # cds_seq, start, end, var_start, var_end, ref, alt, expected_seq
        ("ATGCGTAC", 100, 108, 101, 107, "N", "<DEL>", "AC"),  # internal deletion
        ("ATGCGTAC", 100, 108, 100, 108, "N", "<DEL>", ""),  # entire CDS deleted
        ("ATGCGTAC", 100, 108, 90, 104, "N", "<DEL>", "GTAC"),  # deletion starts before CDS
        ("ATGCGTAC", 100, 108, 104, 120, "N", "<DEL>", "ATGC"),  # deletion ends after CDS
        ("ATGCGTAC", 100, 108, 200, 210, "N", "<DEL>", None),  # deletion outside CDS
    ]

    for cds_seq, start, end, var_start, var_end, ref, alt, expected_seq in cases:
        row = pd.Series(
            {
                "Exon_CDS_seq": cds_seq,
                "Strand": "+",
                "Start": start,
                "End": end,
                "Start_variant": var_start,
                "End_variant": var_end,
                "Ref": ref,
                "Alt": alt,
            }
        )

        result = apply_variant_edge_aware_with_lengths(row)

        if expected_seq is None:
            assert result["Exon_Alt_CDS_seq"] is None or pd.isna(result["Exon_Alt_CDS_seq"])
            assert result["Exon_Alt_CDS_length"] is None or pd.isna(result["Exon_Alt_CDS_length"])
        else:
            assert result["Exon_Alt_CDS_seq"] == expected_seq
            assert result["Exon_Alt_CDS_length"] == len(expected_seq)


def test_apply_variant_edge_aware_with_lengths_with_DUPs():
    cases = [
        # cds_seq, start, end, var_start, var_end, ref, alt, expected_seq
        ("ATGCGTAC", 100, 108, 101, 107, "N", "<DUP>", "ATGCGTATGCGTAC"),  # internal duplication (TGCGTA duplicated)
        ("ATGCGTAC", 100, 108, 100, 108, "N", "<DUP>", "ATGCGTACATGCGTAC"),  # entire CDS duplicated
        (
            "ATGCGTAC",
            100,
            108,
            90,
            104,
            "N",
            "<DUP>",
            "ATGCATGCGTAC",
        ),  # duplication starts before CDS, overlap "ATGC" duplicated
        (
            "ATGCGTAC",
            100,
            108,
            104,
            120,
            "N",
            "<DUP>",
            "ATGCGTACGTAC",
        ),  # duplication ends after CDS, overlap "GTAC" duplicated
        ("ATGCGTAC", 100, 108, 200, 210, "N", "<DUP>", None),  # duplication outside CDS
    ]

    for cds_seq, start, end, var_start, var_end, ref, alt, expected_seq in cases:
        row = pd.Series(
            {
                "Exon_CDS_seq": cds_seq,
                "Strand": "+",
                "Start": start,
                "End": end,
                "Start_variant": var_start,
                "End_variant": var_end,
                "Ref": ref,
                "Alt": alt,
            }
        )

        result = apply_variant_edge_aware_with_lengths(row)

        if expected_seq is None:
            assert result["Exon_Alt_CDS_seq"] is None or pd.isna(result["Exon_Alt_CDS_seq"])
            assert result["Exon_Alt_CDS_length"] is None or pd.isna(result["Exon_Alt_CDS_length"])
        else:
            assert result["Exon_Alt_CDS_seq"] == expected_seq
            assert result["Exon_Alt_CDS_length"] == len(expected_seq)


def test_create_reference_cds_using_file():
    # Load expected output
    expected = pd.read_csv("resources/test_output_files/create_reference_CDS.tsv", sep="\t")

    # Load df3 and cds_df_test from the previous step of your pipeline
    df3 = pd.read_csv("resources/test_output_files/variant_exon_output.tsv", sep="\t")
    cds_df_test = pd.read_csv("resources/test_output_files/cds_df_adj.tsv", sep="\t")
    # the fixture predates the flag; it gave every transcript a stop codon
    cds_df_test["has_stop_codon"] = True

    # Run the function
    actual = create_reference_cds(df3, cds_df_test)

    # Sort both for consistent comparison
    expected_sorted = expected.sort_values(["transcript_id", "variant_id"]).reset_index(drop=True)
    actual_sorted = actual.sort_values(["transcript_id", "variant_id"]).reset_index(drop=True)

    # Compare critical columns
    columns_to_check = [
        "transcript_id",
        "variant_id",
        "ref_cds_seq",
        "alt_cds_seq",
        "ref_cds_len",
        "alt_cds_len",
        "strand",
        "ref",
        "alt",
    ]

    for col in columns_to_check:
        assert all(expected_sorted[col] == actual_sorted[col]), f"Mismatch in column: {col}"


def test_create_reference_cds():
    # create example dataframe
    cds_df_test = pd.DataFrame(
        {
            "transcript_id": ["tx1"] * 5,
            "exon_number": [2, 3, 4, 5, 6],
            "Chromosome": ["chr1"] * 5,
            "gene_id": ["gene1"] * 5,
            "Start": [100, 200, 300, 400, 500],
            "End": [150, 250, 350, 450, 550],
            "Strand": ["+" for _ in range(5)],
            "Exon_CDS_seq": ["AAA", "CCC", "GGG", "TTT", "AAA"],
            "has_stop_codon": [True] * 5,
        }
    )

    # Variant 1: SNP on exon 3
    variant_snp = {
        "transcript_id": "tx1",
        "exon_number": 3,
        "Chromosome": "chr1",
        "gene_id": "gene1",
        "Start": 200,
        "End": 250,
        "Strand": "+",
        "ID": "var_snp",
        "Start_variant": 210,
        "End_variant": 211,
        "Ref": "C",
        "Alt": "T",
        "Exon_Alt_CDS_seq": "CCT",  # Same length
    }

    # Variant 2: Insertion on exon 4
    variant_ins = {
        "transcript_id": "tx1",
        "exon_number": 4,
        "Chromosome": "chr1",
        "gene_id": ["gene1"],
        "Start": 300,
        "End": 350,
        "Strand": "+",
        "ID": "var_ins",
        "Start_variant": 325,
        "End_variant": 325,
        "Ref": "-",
        "Alt": "A",
        "Exon_Alt_CDS_seq": "GGGA",  # Inserted A
    }

    # Variant 3: Deletion on exon 5
    variant_del = {
        "transcript_id": "tx1",
        "exon_number": 5,
        "Chromosome": "chr1",
        "gene_id": ["gene1"],
        "Start": 400,
        "End": 450,
        "Strand": "+",
        "ID": "var_del",
        "Start_variant": 440,
        "End_variant": 441,
        "Ref": "T",
        "Alt": "-",
        "Exon_Alt_CDS_seq": "TT",  # One T removed
    }

    # Variant 4: Deletion spanning exon 3 to 4
    spanning_del_3_4 = [
        {
            "transcript_id": "tx1",
            "exon_number": 3,
            "Chromosome": "chr1",
            "gene_id": ["gene1"],
            "Start": 200,
            "End": 250,
            "Strand": "+",
            "ID": "var_spanning",
            "Start_variant": 240,
            "End_variant": 310,
            "Ref": "CCGG",
            "Alt": "-",
            "Exon_Alt_CDS_seq": "C",  # Shortened version
        },
        {
            "transcript_id": "tx1",
            "exon_number": 4,
            "Chromosome": "chr1",
            "gene_id": ["gene1"],
            "Start": 300,
            "End": 350,
            "Strand": "+",
            "ID": "var_spanning",
            "Start_variant": 240,
            "End_variant": 310,
            "Ref": "CCGG",
            "Alt": "-",
            "Exon_Alt_CDS_seq": "G",  # Shortened version
        },
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
    fasta = {
        "chr1": "AAAAAAAAAACCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    }  # ("A" * 10 + "C" * 20 + "G" * 30 + "T" * 40)

    exons_df = pd.DataFrame(
        [
            {"transcript_id": "tx1", "exon_number": 1, "Chromosome": "chr1", "Start": 10, "End": 13, "Strand": "+"},
            {"transcript_id": "tx1", "exon_number": 2, "Chromosome": "chr1", "Start": 20, "End": 24, "Strand": "+"},
            {"transcript_id": "tx1", "exon_number": 3, "Chromosome": "chr1", "Start": 30, "End": 35, "Strand": "+"},
            {"transcript_id": "tx2", "exon_number": 1, "Chromosome": "chr1", "Start": 40, "End": 43, "Strand": "-"},
            {"transcript_id": "tx2", "exon_number": 2, "Chromosome": "chr1", "Start": 50, "End": 53, "Strand": "-"},
            {"transcript_id": "tx2", "exon_number": 3, "Chromosome": "chr1", "Start": 60, "End": 63, "Strand": "-"},
        ]
    )

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
    expected_tx2_seq = "AAACCCCCC"  # reverse complement of "GGGGGGTTT"
    assert tx2["transcript_sequence"] == expected_tx2_seq
    assert tx2["transcript_length"] == len(expected_tx2_seq)
    assert tx2["transcript_exon_info"] == [(3, 3), (2, 3), (1, 3)]  # reversed for minus strand


def test_get_exon():
    exon_info = [(1, 10), (2, 20), (3, 30)]
    assert get_exon(5, exon_info) == 1
    assert get_exon(25, exon_info) == 2
    assert get_exon(55, exon_info) == 3


def test_analyze_sequence():
    df = pd.DataFrame(
        [
            {
                "ref_cds_seq": "ATGAAATAG",
                "alt_cds_seq": "ATGAAATAA",
                "ref_cds_info": [(1, 9)],
                "alt_cds_info": [(1, 9)],
                "has_stop_codon": True,
            }
        ]
    )
    analyzed = analyze_sequence(df)
    assert analyzed.loc[0, "ref_start_codon_pos"] == 0
    assert analyzed.loc[0, "ref_valid_stop"] == True
    assert analyzed.loc[0, "alt_valid_stop"] == True


def test_analyze_sequence_without_stop_codon():
    # without an annotated stop codon, the real stop lies downstream of the coding region
    df = pd.DataFrame(
        [
            # stop gained in the last codon: TGG>TAG
            {
                "ref_cds_seq": "ATGAAATGG",
                "alt_cds_seq": "ATGAAATAG",
                "ref_cds_info": [(1, 9)],
                "alt_cds_info": [(1, 9)],
                "has_stop_codon": False,
            },
            # coding region that ends in an incomplete codon: the last 3 bases TAA are out of frame
            {
                "ref_cds_seq": "ATGAAAGTAA",
                "alt_cds_seq": "ATGAAAGTAA",
                "ref_cds_info": [(1, 10)],
                "alt_cds_info": [(1, 10)],
                "has_stop_codon": False,
            },
        ]
    )
    analyzed = analyze_sequence(df)

    # every in-frame stop is premature, including one in the last codon
    assert analyzed.loc[0, "alt_first_stop_pos"] == 6
    assert analyzed.loc[0, "alt_is_premature"] == True
    assert analyzed.loc[0, "ref_is_premature"] == False
    # the last codon is not an annotated stop codon, whatever its bases
    assert analyzed["ref_valid_stop"].tolist() == [False, False]
    assert analyzed["alt_valid_stop"].tolist() == [False, False]


def test_start_stop_loss():
    df = pd.DataFrame(
        [
            {
                "ref_start_codon_pos": 0,
                "alt_start_codon_pos": None,
                "ref_valid_stop": True,
                "alt_valid_stop": False,
                "ref_last_codon": "TAG",
                "alt_last_codon": "GGA",
            },
            # stop codon swap: the last codon changes but still encodes a stop
            {
                "ref_start_codon_pos": 0,
                "alt_start_codon_pos": 0,
                "ref_valid_stop": True,
                "alt_valid_stop": True,
                "ref_last_codon": "TAA",
                "alt_last_codon": "TAG",
            },
        ]
    )
    result = start_stop_loss(df)
    assert result["start_loss"].iloc[0] == True
    assert result["stop_loss"].iloc[0] == True
    assert result["start_loss"].iloc[1] == False
    assert result["stop_loss"].iloc[1] == False


def test_splice_alt_cds_into_transcript():
    row = {"ref_cds_seq": "AAAGGGCCC", "alt_cds_seq": "AAATTTCCC"}
    transcript_seq = "TTTAAAGGGCCCGGG"

    result = splice_alt_cds_into_transcript(row, transcript_seq)
    assert result == "TTTAAATTTCCCGGG"


def test_analyze_transcript():

    df = pd.DataFrame(
        [
            {
                "alt_transcript_seq": "CCCATGAAATAATAGGGG",  # ATG at pos 3, TAA at 9, TAG at 12
                "alt_cds_start": 0,
                "transcript_start": 0,
                "transcript_exon_info": [(1, 10), (2, 10)],
                "start_loss": True,
                "stop_loss": False,
            }
        ]
    )

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


# Synthetic transcript in transcript orientation: a 5' UTR, a 48 bp CDS that ends in the sense codon TGG, the stop
# codon TAA and a 3' UTR. Exon 1 holds the 5' UTR and the first 30 CDS bases, exon 2 the rest.
_UTR5 = "GCCGCCACC"
_CDS = "ATGGCTAGCAAAGGCGAAGAGCTGTTCACCGGCGTGGTGCCCATCTGG"
_STOP = "TAA"
_UTR3 = "GGCTGAATTCCCGGG"
_INTRON = "GTAAGTCCCCCCCCTTTCAG"
_FLANK = "CCCCCCCCCC"


def _extract_ptc_synthetic(tmp_path, strand, has_stop_codon, variants):
    """
    Runs extract_ptc on the synthetic transcript and returns the result indexed by variant_id.

    :param strand: strand of the transcript; on the minus strand, the genome is the reverse complement
    :param has_stop_codon: whether the transcript has a stop_codon row. Without it, the transcript ends with its CDS,
        as one tagged cds_end_NF does.
    :param variants: {variant_id: (position, alt)}: SNVs at a 0-based position in the coding region (CDS plus stop
        codon), with the alt base in transcript orientation
    """
    genome = _FLANK + _UTR5 + _CDS[:30] + _INTRON + _CDS[30:] + _STOP + _UTR3 + _FLANK
    cds_start = len(_FLANK) + len(_UTR5)
    exon2_start = cds_start + 30 + len(_INTRON)
    exon2_end = exon2_start + 18 + (len(_STOP) + len(_UTR3) if has_stop_codon else 0)
    rows = [
        ("exon", 1, len(_FLANK), cds_start + 30),
        ("exon", 2, exon2_start, exon2_end),
        ("CDS", 1, cds_start, cds_start + 30),
        ("CDS", 2, exon2_start, exon2_start + 18),
    ]
    if has_stop_codon:
        rows.append(("stop_codon", 2, exon2_start + 18, exon2_start + 21))

    def genomic(pos):
        # 0-based plus strand position of a position in the coding region
        return cds_start + pos if pos < 30 else exon2_start + pos - 30

    coding = _CDS + _STOP
    snvs = [(variant_id, genomic(pos), coding[pos], alt) for variant_id, (pos, alt) in variants.items()]
    if strand == "-":
        length = len(genome)
        genome = str(Seq(genome).reverse_complement())
        rows = [(f, e, length - end, length - start) for f, e, start, end in rows]
        snvs = [(v, length - 1 - g, str(Seq(r).complement()), str(Seq(a).complement())) for v, g, r, a in snvs]

    (tmp_path / "genome.fa").write_text(f">chrT\n{genome}\n")
    gtf_df = pd.DataFrame(
        [
            {
                "Chromosome": "chrT",
                "Start": start,
                "End": end,
                "Strand": strand,
                "Feature": f,
                "exon_number": str(e),
                "transcript_id": "tx",
                "gene_id": "gene",
            }
            for f, e, start, end in rows
        ]
    )
    vcf = pr.PyRanges(
        pd.DataFrame(
            [{"Chromosome": "chrT", "Start": g, "End": g + 1, "ID": v, "Ref": r, "Alt": a} for v, g, r, a in snvs]
        )
    )
    fasta = Fasta(str(tmp_path / "genome.fa"))
    result = extract_ptc(gtf_df[gtf_df["Feature"] != "exon"], vcf, fasta, gtf_df[gtf_df["Feature"] == "exon"])
    return result.set_index("variant_id")


@pytest.mark.parametrize("strand", ["+", "-"])
def test_extract_ptc_without_stop_codon(tmp_path, strand):
    # cds_end_NF: no stop_codon row, the CDS ends in the sense codon TGG
    result = _extract_ptc_synthetic(tmp_path, strand, False, {"TGG>TAG": (46, "A"), "TGG>TGC": (47, "C")})
    variants = ["TGG>TAG", "TGG>TGC"]

    assert result.loc[variants, "has_stop_codon"].tolist() == [False, False]
    assert result.loc[variants, "ref_valid_stop"].tolist() == [False, False]
    # the real stop lies downstream of the CDS, so a stop gained in the last codon is premature
    assert result.loc["TGG>TAG", "alt_is_premature"] == True
    assert result.loc["TGG>TGC", "alt_is_premature"] == False
    # there is no stop codon to lose
    assert result.loc[variants, "stop_loss"].tolist() == [False, False]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_extract_ptc_stop_codon_change(tmp_path, strand):
    variants = {"TAA>TAG": (50, "G"), "TAA>TGA": (49, "G"), "TAA>CAA": (48, "C"), "TGG>TAG": (46, "A")}
    result = _extract_ptc_synthetic(tmp_path, strand, True, variants)

    # a swap to another stop codon keeps the stop codon at its position: no stop loss and no readthrough
    for swap in ["TAA>TAG", "TAA>TGA"]:
        assert result.loc[swap, "alt_valid_stop"] == True
        assert result.loc[swap, "stop_loss"] == False
        assert result.loc[swap, "alt_is_premature"] == False
        assert pd.isna(result.loc[swap, "transcript_num_stop_codons"])
    # the annotated stop codon no longer encodes a stop
    assert result.loc["TAA>CAA", "stop_loss"] == True
    assert result.loc["TAA>CAA", "alt_is_premature"] == False
    # a stop gained in the last sense codon lies upstream of the annotated stop codon
    assert result.loc["TGG>TAG", "alt_is_premature"] == True
    assert result.loc["TGG>TAG", "stop_loss"] == False
