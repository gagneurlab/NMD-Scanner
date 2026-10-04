# Import dependencies
import pandas as pd
import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

import nmd_scanner.rules
from nmd_scanner.extra_features import add_nmd_features, evaluate_nmd_escape_rules
from nmd_scanner.rules import (
    analyze_sequence,
    analyze_transcript,
    apply_variant_edge_aware_with_lengths,
    cds_range_in_transcript,
    create_reference_cds,
    extract_ptc,
    get_exon,
    get_transcript_sequence,
    join_variants_to_cds,
    splice_alt_cds_into_transcript,
    start_stop_loss,
)


def run_pipeline_on_transcript(tmp_path, strand, exon_seqs, cds_range, variant, stop_codon=True):
    """
    Run extract_ptc, add_nmd_features and evaluate_nmd_escape_rules on one synthetic transcript with one variant.

    The genome holds the exons, separated by introns of 20 nt. On the minus strand, it holds their reverse complement.
    :param exon_seqs: Exon sequences in transcript order (5' to 3')
    :param cds_range: (start, end) of the CDS in transcript coordinates, stop codon excluded
    :param variant: (position, ref, alt) in transcript coordinates and orientation, within one exon
    :param stop_codon: Whether the coding region includes the 3 nt after the CDS as its stop codon. Without them, the
        transcript has no annotated stop codon, as one tagged cds_end_NF.
    :return: The single result row as a dictionary
    """

    flank, intron = "C" * 10, "C" * 20
    chrom = f"chr_{tmp_path.name}"  # unique name, since catch_sequence caches sequences across tests

    # Lay out the transcript 5' to 3' and record each exon as (transcript start, layout start, length)
    layout = flank
    exons = []
    tx_pos = 0
    for exon_seq in exon_seqs:
        exons.append((tx_pos, len(layout), len(exon_seq)))
        layout += exon_seq + intron
        tx_pos += len(exon_seq)
    layout += flank

    # On the minus strand, the layout is the reverse complement of the genome
    if strand == "+":
        genome = layout

        def to_genome(start, end):
            return start, end
    else:
        genome = str(Seq(layout).reverse_complement())

        def to_genome(start, end):
            return len(layout) - end, len(layout) - start

    (tmp_path / "genome.fa").write_text(f">{chrom}\n{genome}\n")
    fasta = Fasta(str(tmp_path / "genome.fa"))

    # Exon rows and the coding regions: CDS rows that include the stop codon, split at the exon boundaries
    coding_end = cds_range[1] + 3 if stop_codon else cds_range[1]
    rows = []
    for number, (tx_start, layout_start, length) in enumerate(exons, start=1):
        rows.append(("exon", number, *to_genome(layout_start, layout_start + length)))
        part_start, part_end = max(cds_range[0], tx_start), min(coding_end, tx_start + length)
        if part_start < part_end:
            start, end = to_genome(layout_start + part_start - tx_start, layout_start + part_end - tx_start)
            rows.append(("CDS", number, start, end))
    annotation = pd.DataFrame(
        [
            {
                "Chromosome": chrom,
                "Start": start,
                "End": end,
                "Strand": strand,
                "Feature": feature,
                "exon_number": str(number),
                "transcript_id": "tx1",
                "gene_id": "gene1",
            }
            for feature, number, start, end in rows
        ]
    )

    # VCF alleles are on the plus strand
    position, ref, alt = variant
    tx_start, layout_start, _ = next(exon for exon in reversed(exons) if exon[0] <= position)
    start, end = to_genome(layout_start + position - tx_start, layout_start + position - tx_start + len(ref))
    if strand == "-":
        ref, alt = str(Seq(ref).reverse_complement()), str(Seq(alt).reverse_complement())
    vcf = pd.DataFrame([{"Chromosome": chrom, "Start": start, "End": end, "ID": "var1", "Ref": ref, "Alt": alt}])

    coding = annotation[annotation["Feature"] == "CDS"].assign(has_stop_codon=stop_codon)
    results = extract_ptc(coding, vcf, fasta, annotation[annotation["Feature"] == "exon"])
    assert len(results) == 1
    row = results.iloc[0].to_dict()
    row.update(add_nmd_features(row))
    row.update(evaluate_nmd_escape_rules(row))
    return row


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


def test_create_reference_cds_carries_has_stop_codon():
    # has_stop_codon is per transcript: tx_stop ends in its stop codon TAA, tx_nf (e.g. cds_end_NF) has none
    cds_df_test = pd.DataFrame(
        {
            "transcript_id": ["tx_stop", "tx_stop", "tx_nf"],
            "exon_number": [1, 2, 1],
            "Chromosome": ["chr1"] * 3,
            "gene_id": ["gene1", "gene1", "gene2"],
            "Start": [100, 200, 500],
            "End": [103, 206, 506],
            "Strand": ["+"] * 3,
            "Exon_CDS_seq": ["ATG", "AAATAA", "ATGAAA"],
            "has_stop_codon": [True, True, False],
        }
    )
    variant = {"Chromosome": "chr1", "Strand": "+", "Ref": "A", "Alt": "C"}
    variants = pd.DataFrame(
        [
            # A>C at the second base of exon 2 of tx_stop
            {
                **variant,
                "transcript_id": "tx_stop",
                "exon_number": 2,
                "gene_id": "gene1",
                "Start": 200,
                "End": 206,
                "ID": "var_stop",
                "Start_variant": 201,
                "End_variant": 202,
                "Exon_Alt_CDS_seq": "ACATAA",
            },
            # A>C at the fifth base of exon 1 of tx_nf
            {
                **variant,
                "transcript_id": "tx_nf",
                "exon_number": 1,
                "gene_id": "gene2",
                "Start": 500,
                "End": 506,
                "ID": "var_nf",
                "Start_variant": 504,
                "End_variant": 505,
                "Exon_Alt_CDS_seq": "ATGACA",
            },
        ]
    )

    result = create_reference_cds(variants, cds_df_test)

    assert dict(zip(result["variant_id"], result["alt_cds_seq"])) == {"var_stop": "ATGACATAA", "var_nf": "ATGACA"}
    assert dict(zip(result["variant_id"], result["has_stop_codon"])) == {"var_stop": True, "var_nf": False}


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


def test_cds_range_in_transcript():
    def cds_range(strand, exons, cds):
        exons_df = pd.DataFrame([{"Start": start, "End": end, "Strand": strand} for start, end in exons])
        cds_df = pd.DataFrame([{"Start": start, "End": end, "Strand": strand} for start, end in cds])
        return cds_range_in_transcript(exons_df, cds_df)

    # Plus strand, exons of 100/300/100 nt, CDS with stop codon inside exon 2 at transcript positions 150 to 330
    assert cds_range("+", [(1000, 1100), (1150, 1450), (1500, 1600)], [(1200, 1380)]) == (150, 330)

    # Minus strand, same transcript: exon 1 is the exon with the largest coordinates
    assert cds_range("-", [(1000, 1100), (1150, 1450), (1500, 1600)], [(1220, 1400)]) == (150, 330)

    # Single exon transcripts
    assert cds_range("+", [(5000, 5150)], [(5050, 5110)]) == (50, 110)
    assert cds_range("-", [(4950, 5100)], [(5000, 5060)]) == (40, 100)

    # Minus strand, exons of 50/60/100/200/80/70/90 nt. Exons 1, 2, 6 and 7 hold only UTR.
    # The CDS starts 30 nt into exon 3 and ends 40 nt into exon 5.
    exons = [(1660, 1710), (1590, 1650), (1480, 1580), (1270, 1470), (1180, 1260), (1100, 1170), (1000, 1090)]
    assert cds_range("-", exons, [(1480, 1550), (1270, 1470), (1220, 1260)]) == (140, 450)

    # Stop codon split across exons: 2 nt at the end of exon 1, 1 nt at the start of exon 2
    exons = [(100, 200), (300, 400)]
    assert cds_range("+", exons, [(150, 200), (300, 301)]) == (50, 101)

    # cds_start_NF: the CDS starts at the first base of the transcript
    assert cds_range("-", [(100, 200), (300, 400)], [(150, 200), (300, 400)]) == (0, 150)

    # cds_end_NF: the CDS has no stop codon and runs to the last base of the transcript
    assert cds_range("+", [(100, 200), (300, 400)], [(150, 200), (300, 400)]) == (50, 200)

    # CDS start outside the exons
    assert cds_range("+", [(100, 200)], [(50, 150)]) is None


def test_extract_ptc_locates_cds_in_transcript(tmp_path_factory):
    # Exons of 100/300/100 nt; the CDS with stop codon lies inside exon 2, at transcript positions 150 to 330
    cds_seq = "ATG" + "CAA" * 58 + "TAA"
    transcript_seq = "C" * 150 + cds_seq + "C" * 170
    exon_seqs = [transcript_seq[:100], transcript_seq[100:400], transcript_seq[400:]]

    for strand in ["+", "-"]:
        tmp_path = tmp_path_factory.mktemp(f"cds_in_exon_2_{'plus' if strand == '+' else 'minus'}")
        # Nonsense variant CAA>TAA at CDS position 174, 76 nt before the end of exon 2
        row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (150, 327), (324, "C", "T"))

        assert row["cds_start_in_transcript"] == 150
        assert row["cds_end_in_transcript"] == 330
        assert row["utr5_length"] == 150
        assert row["utr3_length"] == 170
        assert row["alt_transcript_seq"] == transcript_seq[:324] + "T" + transcript_seq[325:]
        assert row["alt_first_stop_pos"] == 174
        assert row["alt_is_premature"] == True
        assert row["ptc_to_intron"] == 76


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
    row = {
        "ref_cds_seq": "AAAGGGCCC",
        "alt_cds_seq": "AAATTTCCC",
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 12,
    }
    transcript_seq = "TTTAAAGGGCCCGGG"

    result = splice_alt_cds_into_transcript(row, transcript_seq)
    assert result == "TTTAAATTTCCCGGG"

    # The 5'UTR repeats the CDS sequence: splice at the CDS position, not at the first match
    row = {
        "ref_cds_seq": "ATGAAATAA",
        "alt_cds_seq": "ATGTAATAA",
        "cds_start_in_transcript": 11,
        "cds_end_in_transcript": 20,
    }
    result = splice_alt_cds_into_transcript(row, "ATGAAATAACCATGAAATAAGG")
    assert result == "ATGAAATAACCATGTAATAAGG"

    # The transcript does not hold the CDS sequence at the CDS position
    assert splice_alt_cds_into_transcript({**row, "cds_start_in_transcript": 10}, "ATGAAATAACCATGAAATAAGG") is None
    # The CDS position is unknown
    row = {**row, "cds_start_in_transcript": None, "cds_end_in_transcript": None}
    assert splice_alt_cds_into_transcript(row, "ATGAAATAACCATGAAATAAGG") is None


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


def _extract_ptc_synthetic(tmp_path, strand, has_stop_codon, variants, split_stop_codon=False):
    """
    Runs extract_ptc on the exon rows and the coding regions (CDS rows that include the stop codon) of the synthetic
    transcript, and returns the result indexed by variant_id.

    :param strand: strand of the transcript; on the minus strand, the genome is the reverse complement
    :param has_stop_codon: whether the coding region ends in the stop codon. Without it, the transcript ends with its
        CDS, as one tagged cds_end_NF does.
    :param variants: {variant_id: (position, alt)}: SNVs at a 0-based position in the coding region (CDS plus stop
        codon), with the alt base in transcript orientation
    :param split_stop_codon: split the stop codon across an intron. Its first 2 bases end exon 2, its last base
        starts exon 3 and is the only base of the CDS row of exon 3.
    """
    stop = _STOP[:2] + _INTRON + _STOP[2:] if split_stop_codon else _STOP
    genome = _FLANK + _UTR5 + _CDS[:30] + _INTRON + _CDS[30:] + stop + _UTR3 + _FLANK
    cds_start = len(_FLANK) + len(_UTR5)
    exon2_start = cds_start + 30 + len(_INTRON)
    stop_start = exon2_start + 18
    exon3_start = stop_start + 2 + len(_INTRON)
    rows = [("exon", 1, len(_FLANK), cds_start + 30), ("CDS", 1, cds_start, cds_start + 30)]
    if split_stop_codon:
        rows += [
            ("exon", 2, exon2_start, stop_start + 2),
            ("exon", 3, exon3_start, exon3_start + 1 + len(_UTR3)),
            ("CDS", 2, exon2_start, stop_start + 2),
            ("CDS", 3, exon3_start, exon3_start + 1),
        ]
    elif has_stop_codon:
        rows += [("exon", 2, exon2_start, stop_start + 3 + len(_UTR3)), ("CDS", 2, exon2_start, stop_start + 3)]
    else:
        rows += [("exon", 2, exon2_start, stop_start), ("CDS", 2, exon2_start, stop_start)]

    def genomic(pos):
        # 0-based plus strand position of a position in the coding region
        if pos < 30:
            return cds_start + pos
        return exon3_start + pos - 50 if split_stop_codon and pos >= 50 else exon2_start + pos - 30

    coding = _CDS + _STOP
    snvs = [(variant_id, genomic(pos), coding[pos], alt) for variant_id, (pos, alt) in variants.items()]
    if strand == "-":
        length = len(genome)
        genome = str(Seq(genome).reverse_complement())
        rows = [(f, e, length - end, length - start) for f, e, start, end in rows]
        snvs = [(v, length - 1 - g, str(Seq(r).complement()), str(Seq(a).complement())) for v, g, r, a in snvs]

    (tmp_path / "genome.fa").write_text(f">chrT\n{genome}\n")
    annotation = pd.DataFrame(
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
    vcf = pd.DataFrame(
        [{"Chromosome": "chrT", "Start": g, "End": g + 1, "ID": v, "Ref": r, "Alt": a} for v, g, r, a in snvs]
    )
    fasta = Fasta(str(tmp_path / "genome.fa"))
    coding = annotation[annotation["Feature"] == "CDS"].assign(has_stop_codon=has_stop_codon)
    result = extract_ptc(coding, vcf, fasta, annotation[annotation["Feature"] == "exon"])
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


@pytest.mark.parametrize("strand", ["+", "-"])
def test_extract_ptc_split_stop_codon(tmp_path, strand):
    # the stop codon TAA is split across an intron; its last base is the only coding base of exon 3
    variants = {"TAA>TAG": (50, "G"), "TAA>CAA": (48, "C")}
    result = _extract_ptc_synthetic(tmp_path, strand, True, variants, split_stop_codon=True)

    assert result.loc["TAA>TAG", "has_stop_codon"] == True
    assert result.loc["TAA>TAG", "ref_cds_seq"] == _CDS + _STOP
    assert result.loc["TAA>TAG", "ref_cds_info"] == [(1, 30), (2, 20), (3, 1)]
    assert result.loc["TAA>TAG", "cds_in_transcript"] == True
    # the variant in exon 3 changes the last base of the stop codon
    assert result.loc["TAA>TAG", "alt_cds_seq"] == _CDS + "TAG"
    assert result.loc["TAA>TAG", "alt_valid_stop"] == True
    assert result.loc["TAA>TAG", "stop_loss"] == False
    assert result.loc["TAA>CAA", "stop_loss"] == True


def test_extract_ptc_needs_the_coding_regions():
    # CDS and stop_codon rows without has_stop_codon: not the coding regions
    rows = pd.DataFrame(
        {
            "Chromosome": ["chrT", "chrT"],
            "Start": [100, 150],
            "End": [150, 153],
            "Strand": ["+", "+"],
            "Feature": ["CDS", "stop_codon"],
            "exon_number": ["1", "1"],
            "transcript_id": ["tx", "tx"],
        }
    )
    with pytest.raises(ValueError, match="has_stop_codon"):
        extract_ptc(rows, vcf=None, fasta=None, exons_df=None)


# join_variants_to_cds: the CDS x VCF join of extract_ptc, with the semantics of the pyranges 0.x join


def _join_cds():
    """CDS rows in an order that is not sorted by Chromosome, and one on the minus strand."""
    return pd.DataFrame(
        {
            "Chromosome": ["chr2", "chr1", "chr1"],
            "Start": [100, 300, 100],
            "End": [200, 400, 200],
            "Strand": ["+", "-", "+"],
            "transcript_id": ["t_chr2", "t_minus", "t_plus"],
        }
    )


def _join_vcf(rows):
    """Variants from (Chromosome, Start, End, ID) tuples."""
    return pd.DataFrame(
        [{"Chromosome": c, "Start": start, "End": end, "ID": i, "Ref": "N", "Alt": "A"} for c, start, end, i in rows]
    )


def _pairs(joined):
    return list(zip(joined["transcript_id"], joined["ID"]))


def test_join_variants_to_cds_uses_half_open_intervals():
    vcf = _join_vcf(
        [
            ("chr1", 99, 100, "ends_at_cds_start"),
            ("chr1", 200, 201, "starts_at_cds_end"),
            ("chr1", 100, 101, "first_base"),
            ("chr1", 199, 200, "last_base"),
            ("chr1", 95, 105, "deletion_over_cds_start"),
            ("chr1", 195, 210, "deletion_over_cds_end"),
            ("chr1", 250, 260, "between_cds_rows"),
        ]
    )
    joined = join_variants_to_cds(_join_cds(), vcf)
    assert sorted(_pairs(joined)) == [
        ("t_plus", "deletion_over_cds_end"),
        ("t_plus", "deletion_over_cds_start"),
        ("t_plus", "first_base"),
        ("t_plus", "last_base"),
    ]


def test_join_variants_to_cds_matches_the_chromosome_and_ignores_the_strand():
    vcf = _join_vcf(
        [
            ("chr2", 150, 151, "on_chr2"),
            ("chr1", 150, 151, "on_chr1"),
            ("chr1", 350, 351, "in_minus_strand_cds"),
            ("chr3", 150, 151, "on_chr3"),
        ]
    )
    joined = join_variants_to_cds(_join_cds(), vcf)
    assert _pairs(joined) == [("t_chr2", "on_chr2"), ("t_minus", "in_minus_strand_cds"), ("t_plus", "on_chr1")]
    assert joined["Chromosome"].tolist() == ["chr2", "chr1", "chr1"]


def test_join_variants_to_cds_without_overlap_gives_no_rows_and_all_columns():
    joined = join_variants_to_cds(_join_cds(), _join_vcf([("chr1", 10, 20, "upstream"), ("chrX", 150, 151, "other")]))
    assert joined.empty
    assert list(joined.columns) == [
        "Chromosome",
        "Start",
        "End",
        "Strand",
        "transcript_id",
        "Start_variant",
        "End_variant",
        "ID",
        "Ref",
        "Alt",
    ]


def test_join_variants_to_cds_suffixes_the_variant_columns_that_cds_df_has_too():
    cds = _join_cds().iloc[[2]].assign(ID="cds_id")
    joined = join_variants_to_cds(cds, _join_vcf([("chr1", 150, 152, "v1")]))
    assert joined.to_dict("records") == [
        {
            "Chromosome": "chr1",
            "Start": 100,
            "End": 200,
            "Strand": "+",
            "transcript_id": "t_plus",
            "ID": "cds_id",
            "Start_variant": 150,
            "End_variant": 152,
            "ID_variant": "v1",
            "Ref": "N",
            "Alt": "A",
        }
    ]
    pd.testing.assert_index_equal(joined.index, pd.RangeIndex(1))


@pytest.mark.parametrize("reverse", [False, True], ids=["pairs_of_polars_bio", "reversed_pairs"])
def test_join_variants_to_cds_order(monkeypatch, reverse):
    """CDS rows in the order of cds_df; the variants of one CDS row by Start, then End descending, then VCF order."""
    if reverse:
        # polars-bio gives the ties in VCF order already, so only reversed pairs test the last sort key
        overlap = nmd_scanner.rules.pb.overlap
        monkeypatch.setattr(nmd_scanner.rules.pb, "overlap", lambda *args, **kwargs: overlap(*args, **kwargs).reverse())
    vcf = _join_vcf(
        [
            ("chr1", 150, 151, "last"),
            ("chr1", 120, 121, "snv_a"),
            ("chr1", 120, 125, "deletion"),
            ("chr1", 120, 121, "snv_b"),
            ("chr1", 90, 110, "first"),
            ("chr1", 350, 351, "minus"),
            ("chr1", 120, 121, "snv_c"),
        ]
    )
    joined = join_variants_to_cds(_join_cds(), vcf)
    assert _pairs(joined) == [
        ("t_minus", "minus"),
        ("t_plus", "first"),
        ("t_plus", "deletion"),
        ("t_plus", "snv_a"),
        ("t_plus", "snv_b"),
        ("t_plus", "snv_c"),
        ("t_plus", "last"),
    ]


def test_join_variants_to_cds_gives_repeated_cds_text_as_category():
    cds = pd.DataFrame(
        {
            "Chromosome": ["chr1"] * 4,
            "Start": [100, 300, 500, 700],
            "End": [200, 400, 600, 800],
            "transcript_id": ["t1", "t2", "t3", "t4"],
        }
    )
    joined = join_variants_to_cds(cds, _join_vcf([("chr1", 150, 151, "v1"), ("chr1", 350, 351, "v2")]))
    assert _pairs(joined) == [("t1", "v1"), ("t2", "v2")]
    assert isinstance(joined["Chromosome"].dtype, pd.CategoricalDtype)
    assert pd.api.types.is_string_dtype(joined["transcript_id"].dtype)
    assert pd.api.types.is_string_dtype(joined["ID"].dtype)


@pytest.mark.parametrize("side", ["cds_df", "vcf"])
def test_join_variants_to_cds_rejects_an_end_above_the_limit_of_polars_bio(side):
    cds, vcf = _join_cds(), _join_vcf([("chr1", 150, 151, "v1")])
    if side == "cds_df":
        cds.loc[0, "End"] = 2**31
    else:
        vcf.loc[0, "End"] = 2**31
    with pytest.raises(ValueError, match=f"Cannot join {side}: it has an End above 2147483647"):
        join_variants_to_cds(cds, vcf)
