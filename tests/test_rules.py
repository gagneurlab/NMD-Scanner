"""
The transcript drawings show a transcript 5' to 3' and are not to scale. `u` is UTR, `=` is CDS, `[...]` is an exon, `|`
between two exons is an exon junction, `*` is the PTC, and `v` marks an indel. The numbers under a transcript are
positions in CDS coordinates, as alt_first_stop_pos: 0 is the first base of the start codon. A row labelled tx gives
transcript coordinates instead. The numbers under an alt transcript are in alt CDS coordinates. `*--->|` is the distance
from the PTC to an exon junction or to the transcript end, and `<--->` is a length. "Technical Notes.md" defines the
features and the NMD escape rules with figures in the same style.
"""

# Import dependencies
import logging

import pandas as pd
import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

import nmd_scanner.rules
from nmd_scanner.extra_features import add_nmd_features, evaluate_nmd_escape_rules
from nmd_scanner.rules import (
    analyze_sequence,
    analyze_transcript,
    annotated_stop_in_alt,
    cds_range_in_transcript,
    create_reference_cds,
    drop_missing_alleles,
    drop_symbolic_alleles,
    extract_ptc,
    get_exon,
    join_variant_windows,
    join_variants_to_cds,
    splice_alt_cds_into_transcript,
    start_stop_loss,
)


def run_pipeline_on_transcript(
    tmp_path, strand, exon_seqs, cds_range, variant, stop_codon=True, start_codon=True, exon_rows=True
):
    """
    Run extract_ptc, add_nmd_features and evaluate_nmd_escape_rules on one synthetic transcript with one variant.

    The genome holds the exons, separated by introns of 20 nt. On the minus strand, it holds their reverse complement.
    :param exon_seqs: Exon sequences in transcript order (5' to 3')
    :param cds_range: (start, end) of the CDS in transcript coordinates, stop codon excluded
    :param variant: (position, ref, alt) in transcript coordinates and orientation, within one exon
    :param stop_codon: Whether the coding region includes the 3 nt after the CDS as its stop codon. Without them, the
        transcript has no annotated stop codon, as one tagged cds_end_NF.
    :param start_codon: Whether the first 3 nt of the CDS are an annotated start codon. Without one, the transcript is
        like one tagged cds_start_NF.
    :param exon_rows: Whether the annotation has the exon rows. Without them, it has only the coding regions.
    :return: The single result row as a dictionary
    """

    flank = "C" * 10
    intron = "C" * 20
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
        part_start = max(cds_range[0], tx_start)
        part_end = min(coding_end, tx_start + length)
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
                # GFF3 phase: the CDS starts with a complete codon
                "Frame": "0",
            }
            for feature, number, start, end in rows
        ]
    )

    # VCF alleles are on the plus strand
    position, ref, alt = variant
    tx_start, layout_start, _ = next(exon for exon in reversed(exons) if exon[0] <= position)
    start, end = to_genome(layout_start + position - tx_start, layout_start + position - tx_start + len(ref))
    if strand == "-":
        ref = str(Seq(ref).reverse_complement())
        alt = str(Seq(alt).reverse_complement())
    vcf = pd.DataFrame([{"Chromosome": chrom, "Start": start, "End": end, "ID": "var1", "Ref": ref, "Alt": alt}])

    coding = annotation[annotation["Feature"] == "CDS"].assign(has_start_codon=start_codon, has_stop_codon=stop_codon)
    exons_df = annotation[annotation["Feature"] == "exon"]
    results = extract_ptc(coding, vcf, fasta, exons_df if exon_rows else exons_df.iloc[:0])
    assert len(results) == 1
    row = results.iloc[0].to_dict()
    row.update(add_nmd_features(row))
    row.update(evaluate_nmd_escape_rules(row))
    return row


def test_create_reference_cds_using_file():
    # Load expected output
    expected = pd.read_csv("resources/test_output_files/create_reference_CDS.tsv", sep="\t")

    df3 = pd.read_csv("resources/test_output_files/variant_exon_output.tsv", sep="\t")
    cds_df_test = pd.read_csv("resources/test_output_files/cds_df_adj.tsv", sep="\t")
    # the fixture has no has_start_codon and has_stop_codon columns, and its expected output gives every transcript a
    # start and a stop codon
    cds_df_test["has_start_codon"] = True
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
            "has_start_codon": [True] * 5,
            "has_stop_codon": [True] * 5,
            "Frame": ["0"] * 5,
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
            "has_start_codon": [True] * 3,
            "has_stop_codon": [True, True, False],
            "Frame": ["0"] * 3,
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


def _values(row, columns):
    """The value of each column of a result row; None means null."""
    return {
        column: None if pd.api.types.is_scalar(row[column]) and pd.isna(row[column]) else row[column]
        for column in columns
    }


@pytest.mark.parametrize("strand", ["+", "-"])
def test_extract_ptc_of_a_transcript_without_exon_rows(tmp_path, strand):
    """
    The only transcript of the variant has no exon rows, only its coding regions. Its row has null transcript columns,
    and the flags come from the alt CDS. The missense variant AGC>AGA at `x` changes no stop codon.

           11 nt        10 nt
    5' [uuuu=======]|[====x===uu] 3'   coding regions only, no exon rows
    tx 0    4        11   15  19 21

    The drawing is in transcript orientation, also on the minus strand.
    """
    row = run_pipeline_on_transcript(
        tmp_path, strand, ["GACCATGGATG", "TAAGCTAAGC"], (4, 16), (15, "C", "A"), exon_rows=False
    )

    expected = {
        "ref_cds_seq": "ATGGATGTAAGCTAA",
        "alt_cds_seq": "ATGGATGTAAGATAA",
        "cds_in_transcript": False,
        "transcript_start": None,
        "transcript_end": None,
        "transcript_seq": None,
        "transcript_length": None,
        "cds_start_in_transcript": None,
        "cds_end_in_transcript": None,
        "transcript_exon_info": None,
        "alt_transcript_seq": None,
        "alt_transcript_length": None,
        "alt_cds_start_in_transcript": None,
        "alt_is_premature": False,
        "start_loss": False,
        "stop_loss": False,
        "utr5_length": None,
        "utr3_length": None,
        "total_exon_count": None,
        "annotated_stop_distance": 0,
        "likely_misannotated": True,
        "nmd_escape": False,
    }
    assert _values(row, expected) == expected


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
                "has_start_codon": True,
                "has_stop_codon": True,
                "cds_frame": 0,
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
                "has_start_codon": True,
                "has_stop_codon": False,
                "cds_frame": 0,
            },
            # coding region that ends in an incomplete codon: the last 3 bases TAA are out of frame
            {
                "ref_cds_seq": "ATGAAAGTAA",
                "alt_cds_seq": "ATGAAAGTAA",
                "ref_cds_info": [(1, 10)],
                "alt_cds_info": [(1, 10)],
                "has_start_codon": True,
                "has_stop_codon": False,
                "cds_frame": 0,
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


def test_start_loss_judges_the_annotated_start_codon():
    """
    Only a variant that changes the annotated start codon, the first 3 nt of the CDS, causes a start loss. The start
    codon can be a non-ATG codon such as CTG. A CDS without an annotated start codon, e.g. of a cds_start_NF
    transcript, has no start codon to lose.

    ref CDS  CTG AAA CCC TAA    annotated start codon CTG
    alt CDS  CCG AAA CCC TAA    CTG>CCG: start loss, the first row
             0   3   6   9
    """
    df = pd.DataFrame(
        {
            "has_start_codon": [True, True, True, True, False],
            "ref_cds_seq": ["CTGAAACCCTAA"] * 2 + ["ATGAAACCCTAA"] * 3,
            "alt_cds_seq": [
                "CCGAAACCCTAA",  # CTG>CCG
                "CTGAGACCCTAA",  # missense AAA>AGA after the start codon CTG
                "ACGAAACCCTAA",  # ATG>ACG
                "ATGCAAACCCTAA",  # insertion after the start codon
                "GTGAAACCCTAA",  # A>G at the first base, without an annotated start codon
            ],
            "ref_cds_info": [[(1, 12)]] * 5,
            "alt_cds_info": [[(1, 12)]] * 3 + [[(1, 13)]] + [[(1, 12)]],
            "has_stop_codon": [True] * 5,
            "cds_frame": [0] * 5,
        }
    )

    result = start_stop_loss(analyze_sequence(df))

    assert result["start_loss"].tolist() == [True, False, True, False, False]


def test_start_codon_pos_is_the_annotated_start_codon():
    """
    The start codon position is the annotated start codon at CDS position 0, also if it is CTG and an in-frame ATG
    follows. Without an annotated start codon, the true start lies upstream of the CDS, and an in-frame ATG is an
    internal Met.

    ref CDS  CTG AAA ATG CCC TAA    annotated start codon CTG: ref_start_codon_pos 0, not 6
             0       6
    """
    df = pd.DataFrame(
        {
            # the last 2 CDS have no annotated start codon (e.g. cds_start_NF)
            "has_start_codon": [True, True, False, False],
            "ref_cds_seq": ["CTGAAAATGCCCTAA", "ATGAAAATGCCCTAA", "ATGAAAATGCCCTAA", "CTGAAAATGCCCTAA"],
            "alt_cds_seq": [
                "CTGAAAATGCCATAA",  # CCC>CCA
                "ACGAAAATGCCCTAA",  # ATG>ACG
                "ATGAAAATGCCATAA",  # CCC>CCA
                "CTGAAAATGCCATAA",  # CCC>CCA
            ],
            "ref_cds_info": [[(1, 15)]] * 4,
            "alt_cds_info": [[(1, 15)]] * 4,
            "has_stop_codon": [True] * 4,
            "cds_frame": [0] * 4,
        }
    )

    result = analyze_sequence(df)

    assert result["ref_start_codon_pos"].tolist() == [0, 0, None, None]
    assert result["ref_start_codon_exon"].tolist() == [1, 1, None, None]
    # ATG>ACG changes the annotated start codon, so the alt CDS has none
    assert result["alt_start_codon_pos"].tolist() == [0, None, None, None]


def test_splice_alt_cds_into_transcript():
    # Single exon transcripts; the variant lies inside the CDS
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
    unknown = {**row, "cds_start_in_transcript": None, "cds_end_in_transcript": None}
    assert splice_alt_cds_into_transcript(unknown, "ATGAAATAACCATGAAATAAGG") is None

    # A deletion of TAAG at t17 reaches 1 nt past the stop codon into the 3'UTR, which loses that nt too
    deletion = {**row, "alt_cds_seq": "ATGAAA", "utr3_change": ("G", "")}
    assert splice_alt_cds_into_transcript(deletion, "ATGAAATAACCATGAAATAAGG") == "ATGAAATAACCATGAAAG"
    # The transcript does not hold the ref 3'UTR bases
    assert splice_alt_cds_into_transcript({**deletion, "utr3_change": ("C", "")}, "ATGAAATAACCATGAAATAAGG") is None

    # An indel at the start codon edge changes the 5'UTR: here CC before the start codon becomes GGG
    utr5 = {**row, "utr5_change": ("CC", "GGG")}
    assert splice_alt_cds_into_transcript(utr5, "ATGAAATAACCATGAAATAAGG") == "ATGAAATAAGGGATGTAATAAGG"
    # The transcript does not hold the ref 5'UTR bases
    assert splice_alt_cds_into_transcript({**utr5, "utr5_change": ("AC", "")}, "ATGAAATAACCATGAAATAAGG") is None


def test_analyze_transcript_without_cds_start_in_transcript():
    # Without the CDS position in the transcript, the alt CDS position is unknown too, and the scan does not run
    df = pd.DataFrame(
        [
            {
                "alt_transcript_seq": "CCCATGAAATAATAGGGG",
                "cds_start_in_transcript": None,
                "alt_cds_start_in_transcript": None,
                "transcript_exon_info": [(1, 10), (2, 10)],
                "alt_transcript_exon_info": [(1, 10), (2, 10)],
                "start_loss": True,
                "stop_loss": False,
            }
        ]
    )

    row = analyze_transcript(df).loc[0]

    assert row["transcript_start_codon_pos"] is None
    assert row["transcript_num_stop_codons"] is None
    assert row["transcript_all_stop_codons"] is None


def test_analyze_transcript_reads_from_the_alt_cds_start():
    """
    The variant shortens the 5' UTR by 1 nt and changes the stop codon TAA to CAA. The classification and the scan read
    the alt transcript from the alt CDS start at t2, and find the TGA at t14: a stop loss. From the ref CDS start at
    t3, they would read the TGA at t3 in the wrong frame, as a PTC.

    ref tx  CCC ATG AAA TAA CTG TGA GG
            0   3       9   12
    alt tx  CC ATG AAA CAA CTG TGA GG
            0  2       8       14
    """
    row = {
        "transcript_seq": "CCCATGAAATAACTGTGAGG",
        "alt_transcript_seq": "CCATGAAACAACTGTGAGG",
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 12,
        "alt_cds_start_in_transcript": 2,
        "cds_frame": 0,
        "has_start_codon": True,
        "has_stop_codon": True,
        "ref_cds_seq": "ATGAAATAA",
        "alt_cds_seq": "ATGAAACAA",
        "transcript_exon_info": [(1, 20)],
        "alt_transcript_exon_info": [(1, 19)],
        "alt_is_premature": False,
        "start_loss": False,
        "stop_loss": True,
    }

    result = analyze_transcript(pd.DataFrame([row])).loc[0]

    assert result["alt_is_premature"] == False
    assert result["stop_loss"] == True
    assert result["transcript_start_codon_pos"] == 2
    assert result["transcript_first_stop_pos"] == 14
    assert result["transcript_all_stop_codons"] == [(14, "TGA")]


# Transcript parts for the scan tests: a 5'UTR of 13 nt with an ATG at transcript position 2, and a 3'UTR from position
# 34 with a TAG at 40 in frame with the CDS and a TAA at 45 out of frame.
_SCAN_UTR5 = "CCATGCCGCCGCC"
_SCAN_UTR3 = "GCTGCTTAGCCTAA" + "C" * 31


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize("exon_starts", [[25], [5, 25]], ids=["no_5utr_intron", "5utr_intron"])
def test_stop_loss_scan_starts_at_the_cds_start(tmp_path, strand, exon_starts):
    """
    After a stop loss, the scan reads on from the CDS start in the transcript, in the frame of the CDS.

    The CDS starts at transcript position 13. With the intron in the 5' UTR, the genomic distance from the transcript
    start to the CDS start is 33 nt. On the minus strand, it is the 45 nt of the 3' UTR. Both give the wrong frame. `x`
    is the lost stop codon TGA>CGA, and `s` is the first stop codon in the frame of the CDS, the TAG at 40.

       5 nt    20 nt            54 nt
    5' [uu]|[uuu======]|[==xxxuuusssuuuuuuuuu] 3'
    tx 0    5   13      25 31 34 40          79

    Without the 5' UTR intron, exons 1 and 2 are one exon of 25 nt. The drawing is in transcript orientation, also on
    the minus strand.
    """
    transcript_seq = _SCAN_UTR5 + "ATG" + "GCT" * 5 + "TGA" + _SCAN_UTR3
    bounds = [0] + exon_starts + [len(transcript_seq)]
    exon_seqs = [transcript_seq[start:end] for start, end in zip(bounds, bounds[1:])]

    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 31), (31, "T", "C"))

    assert row["stop_loss"] == True
    assert row["cds_start_in_transcript"] == 13
    assert row["alt_cds_start_in_transcript"] == 13
    assert row["transcript_start_codon_pos"] == 13
    assert row["transcript_first_stop_codon"] == "TAG"
    assert row["transcript_first_stop_pos"] == 40
    assert row["transcript_all_stop_codons"] == [(40, "TAG")]


@pytest.mark.parametrize(
    ("ref_seq", "alt_seq", "stop", "expected"),
    [
        # SNV upstream of the stop codon at 6, and SNV in it
        ("ATGGCCTAAGG", "ATGACCTAAGG", 6, (6, 6)),
        ("ATGGCCTAAGG", "ATGGCCTAGGG", 6, (6, 6)),
        # 1 nt deletion upstream
        ("ATGGCCTAAGG", "ATGCCTAAGG", 6, (5, 5)),
        # GCC inserted right before the stop codon: the stop codon moves by 3
        ("ATGGCCTAAGG", "ATGGCCGCCTAAGG", 6, (9, 9)),
        # TAA inserted right before the stop codon TAA: placed 3'-most, the insertion follows the stop codon
        ("ATGGCCTAAGG", "ATGGCCTAATAAGG", 6, (9, 6)),
        # TCC inserted right before the stop codon TAA: placed 3'-most, the insertion follows its T
        ("ATGGCCTAAGG", "ATGGCCTCCTAAGG", 6, (9, 6)),
        # the last sense codon TCC deleted: placed 3'-most, the deletion is CCT and takes the T of the stop codon. The
        # replacement ends before the stop codon in the alt: the stop codon maps to the position anchored on the 3'UTR
        ("ATGTCCTAAGG", "ATGTAAGG", 6, (3, 3)),
        # G inserted inside the stop codon: TAA>TGAA
        ("ATGGCCTAAGG", "ATGGCCTGAAGG", 6, (6, 6)),
        # deletion of CCT, from upstream into the stop codon: the stop codon maps to the position anchored on the 3'UTR
        ("ATGGCCTAAGG", "ATGGAAGG", 6, (3, 3)),
        # CT replaced by AGC: the stop codon maps into the replacement, as far in as in the reference
        ("ATGGCCTAAGG", "ATGGCAGCAAGG", 6, (6, 6)),
    ],
    ids=[
        "snv_upstream",
        "snv_in_stop",
        "deletion_upstream",
        "insertion_before_stop",
        "stop_inserted_before_stop",
        "t_insertion_before_stop",
        "t_deletion_before_stop",
        "insertion_in_stop",
        "deletion_into_stop",
        "delins_into_stop",
    ],
)
def test_annotated_stop_in_alt(ref_seq, alt_seq, stop, expected):
    assert annotated_stop_in_alt(ref_seq, alt_seq, stop) == expected


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
                "Frame": "0",
            }
            for f, e, start, end in rows
        ]
    )
    vcf = pd.DataFrame(
        [{"Chromosome": "chrT", "Start": g, "End": g + 1, "ID": v, "Ref": r, "Alt": a} for v, g, r, a in snvs]
    )
    fasta = Fasta(str(tmp_path / "genome.fa"))
    coding = annotation[annotation["Feature"] == "CDS"].assign(has_start_codon=True, has_stop_codon=has_stop_codon)
    result = extract_ptc(coding, vcf, fasta, annotation[annotation["Feature"] == "exon"])
    return result.set_index("variant_id")


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
def test_extract_ptc_gives_a_row_per_vcf_record(tmp_path, strand):
    """
    Two VCF records, var1 and var2, have the same CHROM, POS, REF and ALT: TGG>TAG at `x`, CDS position 46. The PTC
    TAG at 45 lies upstream of the stop codon `s` at 48. Each record gets its own row, with its ID as variant_id.

            39 nt            36 nt
    5' [uuu==========]|[======x=sssuuuu] 3'
       -9  0           30       48

    The drawing is in transcript orientation, also on the minus strand.
    """
    result = _extract_ptc_synthetic(tmp_path, strand, True, {"var1": (46, "A"), "var2": (46, "A")})

    assert list(result.index) == ["var1", "var2"]
    assert result["alt_is_premature"].tolist() == [True, True]
    assert result["alt_first_stop_pos"].tolist() == [45, 45]
    # Apart from variant_id, the two rows are the same
    assert _values(result.loc["var1"], result.columns) == _values(result.loc["var2"], result.columns)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_extract_ptc_skips_records_with_alt_dot_or_star(tmp_path, strand, caplog):
    """
    Three VCF records change TGG at `x`, CDS position 46: var1 has ALT ".", var2 has ALT "*", and var3 is TGG>TAG.
    ALT "." means no alternate allele, and "*" stands for the bases that an overlapping deletion removes. So var1 and
    var2 change no base. They get no row, and a warning counts them. var3 gives the PTC TAG at 45, upstream of the
    stop codon `s` at 48.

            39 nt            36 nt
    5' [uuu==========]|[======x=sssuuuu] 3'
       -9  0           30       48

    The drawing is in transcript orientation, also on the minus strand.
    """
    with caplog.at_level(logging.WARNING):
        result = _extract_ptc_synthetic(
            tmp_path, strand, True, {"var1": (46, "."), "var2": (46, "*"), "var3": (46, "A")}
        )

    assert list(result.index) == ["var3"]
    assert result.loc["var3", "alt_is_premature"] == True
    assert result.loc["var3", "alt_first_stop_pos"] == 45
    warnings = [record.getMessage() for record in caplog.records if record.levelno == logging.WARNING]
    assert [message for message in warnings if "ALT" in message] == [
        'Skipping 2 variant(s) with ALT "." or "*". Such a record changes no base: "*" stands for the bases that an '
        "overlapping deletion removes, and that deletion comes from its own record."
    ]


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


# drop_symbolic_alleles: the records that extract_ptc skips


@pytest.mark.parametrize(
    "alt",
    ["<DEL>", "<DUP>", "<INS>", "<INV>", "<CNV>", "<DUP:TANDEM>", "<INS:ME:ALU>", "<*>"]
    + ["G]chr2:100]", "]chr2:100]G", "G[chr2:100[", "[chr2:100[G", "G.", ".G", "GTA."],
)
def test_drop_symbolic_alleles_drops_symbolic_alleles_and_breakends(alt):
    vcf = pd.DataFrame({"Alt": ["T", alt, "GT"]}, index=[10, 11, 12])

    assert drop_symbolic_alleles(vcf)["Alt"].tolist() == ["T", "GT"]
    assert drop_symbolic_alleles(vcf).index.tolist() == [10, 12]


@pytest.mark.parametrize("alt", ["T", "GT", "ACGTN", "acgt"])
def test_drop_symbolic_alleles_keeps_sequence_alleles(alt, caplog):
    vcf = pd.DataFrame({"Alt": [alt]})

    pd.testing.assert_frame_equal(drop_symbolic_alleles(vcf), vcf)
    assert caplog.text == ""


# drop_missing_alleles: the records with ALT "." or "*" that extract_ptc skips


@pytest.mark.parametrize("alt", [".", "*"])
def test_drop_missing_alleles_drops_alt_dot_and_star(alt):
    vcf = pd.DataFrame({"Alt": ["T", alt, "GT"]}, index=[10, 11, 12])

    assert drop_missing_alleles(vcf)["Alt"].tolist() == ["T", "GT"]
    assert drop_missing_alleles(vcf).index.tolist() == [10, 12]


@pytest.mark.parametrize("alt", ["T", "GT", "ACGTN", "<*>", ".G", "G."])
def test_drop_missing_alleles_keeps_other_alleles(alt, caplog):
    vcf = pd.DataFrame({"Alt": [alt]})

    pd.testing.assert_frame_equal(drop_missing_alleles(vcf), vcf)
    assert caplog.text == ""


# join_variants_to_cds: the CDS x VCF join of extract_ptc


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
    cds = _join_cds()
    vcf = _join_vcf([("chr1", 150, 151, "v1")])
    if side == "cds_df":
        cds.loc[0, "End"] = 2**31
    else:
        vcf.loc[0, "End"] = 2**31
    with pytest.raises(ValueError, match=f"Cannot join {side}: it has an End above 2147483647"):
        join_variants_to_cds(cds, vcf)


def test_join_variant_windows_reaches_3_bases_past_the_coding_row():
    """
    A variant joins a coding row if its window comes within 3 bases of it: 2 bases for the splice dinucleotide, and 1
    more for an insertion right next to it. `w` marks a window base that joins, `x` one that does not.

    genomic  96 97                100                200    202 203
             x  w   ..........    [==================]    ..  w  x

    The joined rows keep the coordinates of the coding row in Start and End, and those of the VCF record in
    Start_variant and End_variant.
    """
    cds = pd.DataFrame({"Chromosome": ["chr1"], "Start": [100], "End": [200], "transcript_id": ["tx"]})
    windows = [(96, 97, "x_before"), (97, 98, "w_before"), (202, 203, "w_after"), (203, 203, "x_insertion_after")]
    variants = pd.DataFrame(
        [
            {
                "Chromosome": "chr1",
                "Start": start - 1,
                "End": end + 1,
                "ID": variant_id,
                "Window_Start": start,
                "Window_End": end,
            }
            for start, end, variant_id in windows
        ]
    )

    joined = join_variant_windows(cds, variants)

    assert list(joined["ID"]) == ["w_before", "w_after"]
    assert list(joined["Start"]) == [100, 100] and list(joined["End"]) == [200, 200]
    assert list(joined["Start_variant"]) == [96, 201] and list(joined["End_variant"]) == [99, 204]
