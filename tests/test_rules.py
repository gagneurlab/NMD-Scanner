"""
The transcript drawings show a transcript 5' to 3' and are not to scale. `u` is UTR, `=` is CDS, `[...]` is an exon, `|`
between two exons is an exon junction, `*` is the PTC, and `v` marks an indel. The numbers under a transcript are
positions in CDS coordinates, as alt_first_stop_pos: 0 is the first base of the start codon. A row labelled tx gives
transcript coordinates instead. The numbers under an alt transcript are in alt CDS coordinates. `*--->|` is the distance
from the PTC to an exon junction or to the transcript end, and `<--->` is a length. "Technical Notes.md" defines the
features and the NMD escape rules with figures in the same style.
"""

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
    annotated_stop_in_alt,
    cds_range_in_transcript,
    create_reference_cds,
    drop_symbolic_alleles,
    extract_ptc,
    get_exon,
    get_transcript_sequence,
    join_variant_windows,
    join_variants_to_cds,
    splice_alt_cds_into_transcript,
    start_stop_loss,
)


def run_pipeline_on_transcript(tmp_path, strand, exon_seqs, cds_range, variant, stop_codon=True, start_codon=True):
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
        ref, alt = str(Seq(ref).reverse_complement()), str(Seq(alt).reverse_complement())
    vcf = pd.DataFrame([{"Chromosome": chrom, "Start": start, "End": end, "ID": "var1", "Ref": ref, "Alt": alt}])

    coding = annotation[annotation["Feature"] == "CDS"].assign(has_start_codon=start_codon, has_stop_codon=stop_codon)
    results = extract_ptc(coding, vcf, fasta, annotation[annotation["Feature"] == "exon"])
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
    #           100 nt                          300 nt                          100 nt
    #     5' [uuuuuuuuuu]|[uuuuu==============================*=====uuuuuuu]|[uuuuuuuuuu] 3'
    # tx     0           100    150                           324   330     400         500
    # CDS    -150        -50    0                             174   180     250         350
    #                                                         *------------>|  ptc_to_intron = 76
    #        <------------------>  utr5_length = 150
    #                                                               <------------------->  utr3_length = 170
    # Both strands: the drawing is in transcript orientation, 5' to 3'. On the minus strand, the genomic
    # coordinates run the other way.
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


# A transcript of 400 nt: 40 nt of 5'UTR, a CDS of 252 nt with its stop codon, and 108 nt of 3'UTR. The CDS repeats CAA,
# whose other two frames hold no stop codon. Two codon pairs put a TAA into one of these frames: CTA ACA at CDS position
# 180 into the frame after a 1 nt deletion, CCT AAC at CDS position 198 into the frame after a 1 nt insertion.
_PTC_CDS = "ATG" + "CAA" * 59 + "CTAACA" + "CAA" * 4 + "CCTAAC" + "CAA" * 15 + "TAA"
_PTC_TRANSCRIPT = "C" * 40 + _PTC_CDS + "C" * 108
_PTC_CDS_RANGE = (40, 289)

# The inputs of the NMD efficiency model best_model.pkl (see scripts/train_new.ipynb). It cannot score a row in which
# one of them is null.
MODEL_INPUTS = [
    "start_loss",
    "stop_loss",
    "total_exon_count",
    "ptc_less_than_150nt_to_start",
    "nmd_long_exon_rule",
    "nmd_start_proximal_rule",
    "nmd_single_exon_rule",
    "nmd_escape",
    "downstream_exon_count",
    "nmd_last_exon_rule",
    "ptc_to_start_codon",
    "stop_codon_distance",
    "ptc_exon_length",
    "ptc_to_intron",
    "upstream_exon_count",
    "nmd_50nt_penultimate_rule",
    "utr5_length",
    "utr3_length",
    "transcript_length",
]


def _ptc_exons(*exon_ends):
    """Split _PTC_TRANSCRIPT into exons that end at the given transcript positions."""
    exon_starts = [0, *exon_ends[:-1]]
    return [_PTC_TRANSCRIPT[start:end] for start, end in zip(exon_starts, exon_ends)]


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    ("exon_ends", "variant", "ptc_pos", "ptc_to_intron"),
    [
        # Nonsense variant CAA>TAA at CDS position 210, in the last exon: 400 - 40 - 210 = 150
        # Both strands: the drawing is in transcript orientation, 5' to 3'. On the minus strand, the genomic
        # coordinates run the other way.
        #           100 nt       100 nt            200 nt
        #     5' [uuuu======]|[==========]|[=====*===uuuuuuuuuuu] 3'
        #        -40  0      60           160    210 252        360
        #                                        *------------->|  ptc_to_intron = 150
        pytest.param((100, 200, 400), (250, "C", "T"), 210, 150, id="last_exon"),
        #                          400 nt
        #     5' [uuuu=====================*===uuuuuuuuuuu] 3'
        #        -40  0                    210 252        360
        #                                  *------------->|  ptc_to_intron = 150
        pytest.param((400,), (250, "C", "T"), 210, 150, id="single_exon"),
        # A 1 nt deletion in exon 2 makes the TAA at CDS position 181 the PTC, at 180 in alt CDS coordinates.
        # The alt transcript has 399 nt: 399 - 40 - 180 = 179
        #           100 nt       100 nt            200 nt
        #                            v  CA>C deletes the A at CDS position 121
        # ref 5' [uuuu======]|[==========]|[=========uuuuuuuuuuu] 3'
        #        -40  0      60           160        252        360
        # alt 5' [uuuu======]|[==========]|[==*======uuuuuuuuuuu] 3'
        #        -40  0      60           159 180    251        359
        #                                     *---------------->|  ptc_to_intron = 179
        pytest.param((100, 200, 400), (160, "CA", "C"), 180, 179, id="last_exon_after_deletion"),
        # A 1 nt insertion in exon 2 makes the TAA at CDS position 200 the PTC, at 201 in alt CDS coordinates.
        # The alt transcript has 401 nt: 401 - 40 - 201 = 160
        #           100 nt       100 nt            200 nt
        #                            v  C>CA inserts an A after CDS position 120
        # ref 5' [uuuu======]|[==========]|[=========uuuuuuuuuuu] 3'
        #        -40  0      60           160        252        360
        # alt 5' [uuuu======]|[==========]|[====*====uuuuuuuuuuu] 3'
        #        -40  0      60           161   201  253        361
        #                                       *-------------->|  ptc_to_intron = 160
        pytest.param((100, 200, 400), (160, "C", "CA"), 201, 160, id="last_exon_after_insertion"),
    ],
)
def test_ptc_to_intron_of_a_ptc_in_the_last_exon(tmp_path, strand, exon_ends, variant, ptc_pos, ptc_to_intron):
    # The last exon ends at the transcript end, so ptc_to_intron is the length of the 3'UTR that the PTC creates
    row = run_pipeline_on_transcript(tmp_path, strand, _ptc_exons(*exon_ends), _PTC_CDS_RANGE, variant)

    assert row["alt_first_stop_pos"] == ptc_pos
    assert row["alt_is_premature"] == True
    assert row["nmd_last_exon_rule"] == True
    assert row["ptc_to_intron"] == ptc_to_intron
    assert row["ptc_to_intron"] == row["alt_transcript_length"] - row["cds_start_in_transcript"] - ptc_pos


@pytest.mark.parametrize(
    ("exon_ends", "variant", "ptc_to_intron"),
    [
        # Nonsense variant CAA>TAA at CDS position 90, in exon 2 of 3: 200 - 40 - 90 = 70
        #           100 nt       100 nt            200 nt
        #     5' [uuuu======]|[===*======]|[=========uuuuuuuuuuu] 3'
        #        -40  0      60   90      160        252        360
        #                         *------>|  ptc_to_intron = 70
        pytest.param((100, 200, 400), (130, "C", "T"), 70, id="internal_exon"),
        # CAA>TAA at CDS position 210, in exon 3, which holds the end of the CDS. Exon 4 holds only 3'UTR: 320 - 40 - 210
        #           100 nt       100 nt        120 nt        80 nt
        #     5' [uuuu======]|[==========]|[=====*===uuu]|[uuuuuuuu] 3'
        #        -40  0      60           160    210 252 280       360
        #                                        *------>|  ptc_to_intron = 70
        pytest.param((100, 200, 320, 400), (250, "C", "T"), 70, id="last_cds_exon_before_utr_exon"),
        # last_exon and single_exon: see the drawings in test_ptc_to_intron_of_a_ptc_in_the_last_exon
        pytest.param((100, 200, 400), (250, "C", "T"), 150, id="last_exon"),
        pytest.param((400,), (250, "C", "T"), 150, id="single_exon"),
    ],
)
def test_ptc_rows_have_all_model_inputs(tmp_path, exon_ends, variant, ptc_to_intron):
    row = run_pipeline_on_transcript(tmp_path, "+", _ptc_exons(*exon_ends), _PTC_CDS_RANGE, variant)

    assert row["alt_is_premature"] == True
    assert row["ptc_to_intron"] == ptc_to_intron
    assert [name for name in MODEL_INPUTS if pd.isna(row[name])] == []


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


def test_analyze_sequence_reads_codons_in_the_cds_frame():
    """
    A cds_start_NF CDS with phase 1: its first base A belongs to no codon, and the codons TGC AAA CCC TAA follow.

    ref CDS  A TGC AAA CCC TAA
             0 1   4   7   10
    """
    df = pd.DataFrame(
        {
            "ref_cds_seq": ["ATGCAAACCCTAA"] * 3,
            "alt_cds_seq": [
                "GTGCAAACCCTAA",  # A>G at CDS position 0
                "ATGTAAACCCTAA",  # C>T at CDS position 3: TGC>TGT, and TAA out of frame
                "ATGCTAACCCTAA",  # A>T at CDS position 4: AAA>TAA in frame
            ],
            "ref_cds_info": [[(1, 13)]] * 3,
            "alt_cds_info": [[(1, 13)]] * 3,
            "has_start_codon": [False] * 3,
            "has_stop_codon": [True] * 3,
            "cds_frame": [1] * 3,
        }
    )

    result = start_stop_loss(analyze_sequence(df))

    # Without an annotated start codon, the out-of-frame ATG at CDS position 0 is no start codon, so changing it is no
    # start loss
    assert result["ref_start_codon_pos"].tolist() == [None] * 3
    assert result["start_loss"].tolist() == [False] * 3
    # Out-of-frame stop codon: the first in-frame stop codon is still the one at the CDS end
    assert result.loc[1, "alt_is_premature"] == False
    assert result.loc[1, "alt_first_stop_pos"] == 10
    # In-frame stop codon
    assert result.loc[2, "alt_is_premature"] == True
    assert result.loc[2, "alt_first_stop_pos"] == 4


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


def test_start_loss_with_a_deleted_stop_codon_reads_from_the_next_atg_into_the_3utr():
    """
    A deletion removes the start codon and the stop codon of ATG AAA CCC TAA. Only the first base A of the CDS is left.
    After the start loss, the scan takes the next ATG, at t4, and reads its frame on into the 3' UTR, to the TGA at t10.
    This ATG lies in the former 3' UTR, downstream of the deleted stop codon. So its ORF does not overlap the CDS, and
    the row is neither a PTC nor a stop loss.

    ref tx  CC ATG AAA CCC TAA GATGCCCTGACC
            0  2           11
    alt tx  CC A G ATG CCC TGA CC
            0  2   4       10
    """
    row = {
        "transcript_seq": "CC" + "ATGAAACCCTAA" + "GATGCCCTGACC",
        "alt_transcript_seq": "CC" + "A" + "GATGCCCTGACC",
        "cds_start_in_transcript": 2,
        "cds_end_in_transcript": 14,
        "alt_cds_start_in_transcript": 2,
        "cds_frame": 0,
        "has_stop_codon": True,
        "ref_cds_seq": "ATGAAACCCTAA",
        "alt_cds_seq": "A",
        "transcript_exon_info": [(1, 15)],
        "alt_is_premature": False,
        "start_loss": True,
        "stop_loss": True,
    }

    result = analyze_transcript(pd.DataFrame([row])).loc[0]

    assert result["alt_is_premature"] == False
    assert result["stop_loss"] == False
    assert result["transcript_start_codon_pos"] == 4
    assert result["transcript_first_stop_codon"] == "TGA"
    assert result["transcript_first_stop_pos"] == 10


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


def test_analyze_transcript():

    df = pd.DataFrame(
        [
            {
                "alt_transcript_seq": "CCCATGAAATAATAGGGG",  # ATG at pos 3, TAA at 9, TAG at 12
                "transcript_seq": "CCCATGAAATAATAGGGG",
                "cds_start_in_transcript": 0,
                "alt_cds_start_in_transcript": 0,
                "cds_frame": 0,
                "cds_end_in_transcript": 12,
                "has_stop_codon": True,
                "ref_cds_seq": "CCCATGAAATAA",
                "alt_cds_seq": "CCCATGAAATAA",
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


def test_analyze_transcript_without_cds_start_in_transcript():
    # Without the CDS position in the transcript, the alt CDS position is unknown too, and the scan does not run
    df = pd.DataFrame(
        [
            {
                "alt_transcript_seq": "CCCATGAAATAATAGGGG",
                "cds_start_in_transcript": None,
                "alt_cds_start_in_transcript": None,
                "transcript_exon_info": [(1, 10), (2, 10)],
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


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize("exon_starts", [[25], [5, 25]], ids=["no_5utr_intron", "5utr_intron"])
def test_start_loss_scan_starts_at_the_cds_start(tmp_path, strand, exon_starts):
    """
    After a start loss, the scan takes the first ATG at or after the CDS start in the transcript.

    The scan starts at transcript position 13, so it skips the ATG at 2 in the 5' UTR. `x` is the lost start codon
    ATG>ATA, and `a` is the in-frame ATG at 19 that the scan finds. `s` is the first stop codon in its frame, the TGA at
    position 31. With the intron in the 5' UTR or on the minus strand, the genomic distance is 33 or 45 nt. A scan from
    there finds no ATG.

       5 nt          20 nt              54 nt
    5' [uu]|[uuuxxx===aaa===]|[======sssuuuuuuuuu] 3'
    tx 0    5   13    19      25     31          79

    Without the 5' UTR intron, exons 1 and 2 are one exon of 25 nt. The drawing is in transcript orientation, also on
    the minus strand.
    """
    transcript_seq = _SCAN_UTR5 + "ATG" + "GCT" + "ATG" + "GCT" * 3 + "TGA" + _SCAN_UTR3
    bounds = [0] + exon_starts + [len(transcript_seq)]
    exon_seqs = [transcript_seq[start:end] for start, end in zip(bounds, bounds[1:])]

    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 31), (15, "G", "A"))

    assert row["start_loss"] == True
    assert row["cds_start_in_transcript"] == 13
    assert row["alt_cds_start_in_transcript"] == 13
    assert row["transcript_start_codon_pos"] == 19
    assert row["transcript_first_stop_codon"] == "TGA"
    assert row["transcript_first_stop_pos"] == 31


# Transcript for the stop codon classification tests, exons split at t25. A 5'UTR of 13 nt, the CDS from t13 with the
# stop codon TAA at t40, and a 3'UTR from t43. In the 3'UTR, a TGA at t47 lies in the frame shifted by -1 nt,
# a TAG at t52 in the frame of the CDS, and the frame shifted by +1 nt has no stop codon.
# In the frame shifted by -1 nt, the codons CTG ACC at t25 read as TGA.
#
#             25 nt                                   50 nt
#    5' [uuuuuuuuuuuuu============]|[===============sssuuuuaaauubbbuuuuuuuuuuuuuuuuuuuu] 3'
#    tx  0            13             25             40 43  47   52                     75
#
# `s` is the stop codon TAA, `a` the TGA in the frame shifted by -1 nt, and `b` the TAG in the frame of the CDS. The
# drawing is to scale and in transcript orientation, also on the minus strand.
_STOP_UTR5 = "CCGCCGCCACCGC"
_STOP_CDS = "ATG" + "GCC" * 3 + "CTG" + "ACC" + "GCC" * 3
_STOP_UTR3 = "CCCC" + "TGA" + "CC" + "TAG" + "C" * 20
_STOP_TRANSCRIPT = _STOP_UTR5 + _STOP_CDS + "TAA" + _STOP_UTR3


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


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize(
    ("variant", "is_premature", "stop_loss", "alt_first_stop_pos", "new_stop_pos"),
    [
        # 1 nt deletion at t18: the shifted frame reads TGA at alt t25, inside the CDS
        ((17, "CC", "C"), True, False, 12, None),
        # 1 nt deletion at t33: the shifted frame has no stop in the CDS and reads TGA at alt t46, in the 3'UTR
        ((32, "CC", "C"), False, True, None, 46),
        # 1 nt insertion after t32: the shifted frame has no stop before the transcript end (nonstop)
        ((32, "C", "CG"), False, True, None, None),
        # TAA>CAA: reading on in frame, the next stop is the TAG at t52
        ((40, "T", "C"), False, True, None, 52),
        # TAA>TGAA: TGA at the annotated position of the stop codon
        ((40, "T", "TG"), False, False, 27, None),
        # GCC inserted right before the stop codon: the stop codon moves to alt t43
        ((39, "C", "CGCC"), False, False, 30, None),
        # GCC deleted right before the stop codon: the stop codon moves to alt t37
        ((36, "CGCC", "C"), False, False, 24, None),
        # TAA inserted right before the stop codon TAA: the same alt transcript as TAA inserted after it
        ((39, "C", "CTAA"), False, False, 27, None),
        # TCCTAG inserted right before the stop codon: the TAG at alt t43 lies upstream of the stop codon at alt t46
        ((39, "C", "CTCCTAG"), True, False, 30, None),
    ],
    ids=[
        "frameshift_ptc_in_cds",
        "frameshift_stop_in_utr3",
        "frameshift_nonstop",
        "stop_codon_snv",
        "insertion_in_stop_codon",
        "inframe_insertion_before_stop_codon",
        "inframe_deletion_before_stop_codon",
        "stop_codon_inserted_before_stop_codon",
        "stop_codon_gained_before_stop_codon",
    ],
)
def test_first_stop_codon_classification(
    tmp_path, strand, variant, is_premature, stop_loss, alt_first_stop_pos, new_stop_pos
):
    exon_seqs = [_STOP_TRANSCRIPT[:25], _STOP_TRANSCRIPT[25:]]
    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 40), variant)

    assert row["alt_is_premature"] == is_premature
    assert row["stop_loss"] == stop_loss
    assert row["alt_first_stop_pos"] == alt_first_stop_pos
    if stop_loss:
        # A stop loss reports the next stop codon in frame in the readthrough columns, for a frameshift as for an SNV
        assert row["transcript_first_stop_pos"] == new_stop_pos
        assert row["transcript_num_stop_codons"] == (0 if new_stop_pos is None else 1)
    else:
        assert row["transcript_first_stop_pos"] is None
    if not (is_premature or stop_loss):
        # The first stop codon is the annotated one
        assert row["stop_codon_distance"] == 0


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize(
    ("variant", "distance"),
    [
        # 1 nt deletion at t18: the stop codon moves to alt t39, and the shifted frame reads TGA at alt t25
        ((17, "CC", "C"), 14),
        # TCCTAG inserted right before the stop codon: TAG at alt t43, the stop codon at alt t46
        ((39, "C", "CTCCTAG"), 3),
        # TAA>TGAA: TGA at the annotated position of the stop codon
        ((40, "T", "TG"), 0),
        # TAA>CAA: reading on in frame, the next stop is the TAG at t52
        ((40, "T", "C"), -12),
        # 1 nt deletion at t33: the stop codon moves to alt t39, and the shifted frame reads TGA at alt t46
        ((32, "CC", "C"), -7),
        # 1 nt insertion after t32: no stop codon up to the transcript end (nonstop)
        ((32, "C", "CG"), None),
    ],
    ids=[
        "frameshift_ptc_in_cds",
        "stop_codon_gained_before_stop_codon",
        "insertion_in_stop_codon",
        "stop_codon_snv",
        "frameshift_stop_in_utr3",
        "frameshift_nonstop",
    ],
)
def test_stop_codon_distance(tmp_path, strand, variant, distance):
    """
    The distance runs from the first in-frame stop codon of the alt transcript to the annotated stop codon, both in
    alt transcript coordinates. It is negative for a stop loss: the new stop codon lies downstream.

    The 1 nt deletion at t18 (`v`) moves the stop codon `s` to alt t39. The shifted frame reads the PTC `*`, a TGA at
    alt t25:

                          v
    5' [uuuuuuuuuuuuu===========]|[=***===========sssuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu] 3'
    tx  0            13            24             39                                 74
                                    *------------>  stop_codon_distance = 14

    TAA>CAA (`x`) loses the stop codon at t40. Reading on in frame, the next stop codon `s` is the TAG at t52:

    5' [uuuuuuuuuuuuu============]|[===============xxxuuuuuuuuusssuuuuuuuuuuuuuuuuuuuu] 3'
    tx  0            13             25             40          52                     75
                                                   <---------->  stop_codon_distance = -12

    Both drawings are to scale and in transcript orientation, also on the minus strand.
    """
    exon_seqs = [_STOP_TRANSCRIPT[:25], _STOP_TRANSCRIPT[25:]]
    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 40), variant)

    assert row["stop_codon_distance"] == distance


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize(
    ("last_codon", "variant"),
    [
        ("GCC", (37, "G", "GCCT")),
        ("GCC", (39, "C", "CTCC")),
        ("GCC", (40, "T", "TCCT")),
        ("GCC", (39, "C", "CTGGCCC")),
        ("TCC", (36, "CTCC", "C")),
        ("TCC", (37, "TCCT", "T")),
    ],
    ids=[
        "tcc_inserted_after_t37",
        "tcc_inserted_after_t39",
        "tcc_inserted_after_t40",
        "tggccc_inserted",
        "tcc_deleted_after_t36",
        "tcc_deleted_after_t37",
    ],
)
def test_inframe_indel_before_stop_codon_starting_with_t(tmp_path, strand, last_codon, variant):
    # The CDS ends in last_codon, followed by the stop codon TAA at t40. An in-frame indel right before the stop codon
    # starts with T, as TAA does. So it can also be placed after the T of the stop codon: TCC inserted before TAA
    # (after t39) gives the same alt transcript as CCT inserted after T (after t40) or CCT after t37. Deleting the last
    # sense codon TCC equals deleting CCT from t38. The alt transcript reads TAA at the shifted position of the stop
    # codon: neither a PTC nor a stop loss.
    transcript_seq = _STOP_UTR5 + _STOP_CDS[:-3] + last_codon + "TAA" + _STOP_UTR3
    exon_seqs = [transcript_seq[:25], transcript_seq[25:]]
    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 40), variant)

    assert row["alt_is_premature"] == False
    assert row["stop_loss"] == False
    assert row["nmd_escape"] == False
    assert row["stop_codon_distance"] == 0


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize("variant", [(36, "C", "T"), (20, "C", "A")], ids=["synonymous", "missense"])
def test_annotated_stop_codon_without_stop_in_reference(tmp_path, strand, variant):
    # The stop_codon rows lie on TCA, which is no stop codon: the annotation does not match the genome. The reference
    # transcript reads its first in-frame stop codon at t52, in the 3'UTR. The row keeps the flags from the CDS: no
    # stop codon to lose, and no PTC.
    transcript_seq = _STOP_UTR5 + _STOP_CDS + "TCA" + _STOP_UTR3
    exon_seqs = [transcript_seq[:25], transcript_seq[25:]]
    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 40), variant)

    assert row["ref_valid_stop"] == False
    assert row["alt_is_premature"] == False
    assert row["stop_loss"] == False
    assert row["stop_codon_distance"] is None


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize(
    ("variant", "stop_loss"),
    [((40, "T", "C"), True), ((36, "C", "T"), False)],
    ids=["stop_codon_snv", "synonymous"],
)
def test_internal_stop_codon_in_reference(tmp_path, strand, variant, stop_loss):
    # The CDS holds an in-frame TGA at t28, as a selenoprotein holds a selenocysteine codon. The pipeline does not
    # read selenocysteine annotations, so it cannot tell TGA from a misannotated stop codon. The reference transcript
    # does not read through to the annotated stop codon, so the row keeps the flags from the CDS: the TGA is the first
    # stop codon of the alt CDS, and TAA>CAA loses the annotated stop codon.
    transcript_seq = _STOP_UTR5 + _STOP_CDS[:15] + "TGA" + _STOP_CDS[18:] + "TAA" + _STOP_UTR3
    exon_seqs = [transcript_seq[:25], transcript_seq[25:]]
    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 40), variant)

    assert row["alt_is_premature"] == True
    assert row["alt_first_stop_pos"] == 15
    assert row["stop_loss"] == stop_loss
    assert row["stop_codon_distance"] == 12


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize(
    "variant",
    [(40, "TAACC", "T"), (41, "AACCC", "C")],
    ids=["anchor_in_cds", "anchor_in_utr3"],
)
def test_deletion_from_stop_codon_into_utr3(tmp_path, strand, variant):
    # Both variants delete t41 to t44: the last 2 nt of the stop codon and the first 2 nt of the 3'UTR.
    # On the plus strand, a VCF anchors the deletion at t40, in the CDS; on the minus strand at t45, in the 3'UTR.
    exon_seqs = [_STOP_TRANSCRIPT[:25], _STOP_TRANSCRIPT[25:]]
    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 40), variant)

    # The alt transcript lacks the 3'UTR part of the deletion too. Reading on in frame, TCC at t40 is followed by TGA.
    assert row["alt_transcript_seq"] == _STOP_TRANSCRIPT[:41] + _STOP_TRANSCRIPT[45:]
    assert row["alt_is_premature"] == False
    assert row["stop_loss"] == True
    assert row["transcript_first_stop_codon"] == "TGA"
    assert row["transcript_first_stop_pos"] == 43


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize(
    "variant",
    [(38, "CCTAAC", "C"), (39, "CTAACC", "C")],
    ids=["from_cds", "into_utr3"],
)
def test_deletion_across_stop_codon_placements(tmp_path, strand, variant):
    # Here the 3'UTR starts with CCTAG. Deleting t39 to t43 (C, the stop codon TAA, C) gives the same alt transcript
    # as deleting t40 to t44 (the stop codon TAA, CC). The CDS keeps its last codon GCC, and the TAG from the 3'UTR
    # lands at t40, the annotated position of the stop codon. Both placements are neither a PTC nor a stop loss.
    transcript_seq = _STOP_UTR5 + _STOP_CDS + "TAA" + "CCTAG" + "C" * 20
    exon_seqs = [transcript_seq[:25], transcript_seq[25:]]
    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 40), variant)

    assert row["alt_transcript_seq"] == transcript_seq[:39] + transcript_seq[44:]
    assert row["alt_is_premature"] == False
    assert row["stop_loss"] == False
    assert row["stop_codon_distance"] == 0


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize(
    "variant",
    [(40, "TA", "T"), (41, "AA", "A"), (42, "AA", "A")],
    ids=["ref_t40_t41", "ref_t41_t42", "ref_t42_t43"],
)
def test_deletion_in_stop_codon_run(tmp_path, strand, variant):
    # The CDS ends in TGG TAA, and the 3'UTR starts with ACTG. Each variant deletes one A of the run AAA at t41 to t43,
    # which leaves TAA at t40. On the plus strand, ref_t40_t41 is the left-normalized TA>T. On the minus strand,
    # ref_t42_t43 is anchored on the 3'UTR base t43. The alt CDS can end in GTA, but the alt transcript still reads TAA
    # at the annotated position: neither a PTC nor a stop loss.
    transcript_seq = _STOP_UTR5 + _STOP_CDS[:-3] + "TGG" + "TAA" + "ACTG" + "C" * 20
    exon_seqs = [transcript_seq[:25], transcript_seq[25:]]
    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 40), variant)

    assert row["alt_transcript_seq"] == transcript_seq[:43] + transcript_seq[44:]
    assert row["alt_is_premature"] == False
    assert row["stop_loss"] == False
    assert row["stop_codon_distance"] == 0


# Deletions and delins that remove the start of the stop codon, on the transcript of the stop codon tests. The last
# codon, the stop codon and the 3'UTR start are given; the variant is (position, length of ref, alt) in transcript
# coordinates, a deletion with alt "" and its first deleted base as position. The stop codon is lost if every
# representation of the variant removes it. A delins is not split into an SNV plus an indel.
_STOP_START_REMOVED = [
    # CCTA>G: TCC TAA CCCC gives TGA CCCC, the 3'UTR unchanged: the stop codon is still there
    ("TCC", "TAA", "CCCC", (38, 4, "G"), None),
    ("TCT", "TAA", "CCCC", (38, 4, "G"), None),
    ("TCC", "TAA", "TAACC", (38, 4, "G"), None),
    # the last sense codon and the stop codon deleted, with and without a stop codon repeat
    ("TCC", "TAA", "TAACC", (37, 6, ""), -3),
    ("GTA", "TAG", "TAGCAT", (37, 6, ""), -3),
    ("GTA", "TAG", "TAGCAT", (36, 6, ""), -3),
    ("GCT", "TAA", "CTAAGCC", (36, 7, ""), -4),
    # ATAG>T: a synonymous change plus a stop codon deletion in a stop codon repeat. A delins is not split
    ("GTA", "TAG", "TAGCAT", (39, 4, "T"), -3),
]
_STOP_START_REMOVED_IDS = [
    "cctag_to_g",
    "tct_taa_to_tga",
    "cctag_to_g_before_stop_repeat",
    "last_codon_and_stop_deleted_in_repeat",
    "last_codon_and_stop_deleted_in_stop_repeat",
    "last_codon_and_stop_deleted_shifted_left",
    "seven_nt_deleted",
    "atag_to_t_in_stop_repeat",
]


def _run_stop_start_removed(tmp_path, strand, last_codon, stop_codon, utr3, variant):
    transcript_seq = _STOP_UTR5 + _STOP_CDS[:-3] + last_codon + stop_codon + utr3 + "C" * 20
    exon_seqs = [transcript_seq[:25], transcript_seq[25:]]
    position, ref_length, alt = variant
    if alt == "":
        # VCF anchors a deletion on the base before it
        position, ref_length, alt = position - 1, ref_length + 1, transcript_seq[position - 1]
    ref = transcript_seq[position : position + ref_length]
    return run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 40), (position, ref, alt))


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize(
    ("last_codon", "stop_codon", "utr3", "variant", "distance"), _STOP_START_REMOVED, ids=_STOP_START_REMOVED_IDS
)
def test_variant_removing_the_start_of_the_stop_codon(
    tmp_path, strand, last_codon, stop_codon, utr3, variant, distance
):
    row = _run_stop_start_removed(tmp_path, strand, last_codon, stop_codon, utr3, variant)

    # Without a distance, the stop codon stays at the annotated position, followed by the unchanged 3'UTR: neither a
    # PTC nor a stop loss.
    assert row["alt_is_premature"] == False
    assert row["stop_loss"] == (distance is not None)
    assert row["nmd_escape"] == False


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize(
    ("last_codon", "stop_codon", "utr3", "variant", "distance"), _STOP_START_REMOVED, ids=_STOP_START_REMOVED_IDS
)
def test_stop_codon_distance_of_variant_removing_the_start_of_the_stop_codon(
    tmp_path, strand, last_codon, stop_codon, utr3, variant, distance
):
    row = _run_stop_start_removed(tmp_path, strand, last_codon, stop_codon, utr3, variant)

    assert row["stop_codon_distance"] == (0 if distance is None else distance)


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
@pytest.mark.parametrize(
    ("variant", "is_premature"),
    [((17, "CC", "C"), True), ((32, "CC", "C"), False), ((32, "C", "CG"), False)],
    ids=["frameshift_ptc_in_cds", "frameshift_stop_past_cds", "frameshift_nonstop"],
)
def test_stop_codon_after_frameshift_without_annotated_stop(tmp_path, strand, variant, is_premature):
    # No stop_codon rows: the CDS ends at t40, and the TAA there is not annotated.
    # A stop inside the CDS is premature. A stop past its end, or none, is neither premature nor a stop loss.
    exon_seqs = [_STOP_TRANSCRIPT[:25], _STOP_TRANSCRIPT[25:]]
    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (13, 40), variant, stop_codon=False)

    assert row["alt_is_premature"] == is_premature
    assert row["stop_loss"] == False


@pytest.mark.parametrize("strand", ["+", "-"], ids=["plus", "minus"])
def test_stop_loss_with_stop_codon_out_of_frame(tmp_path, strand):
    # The CDS starts at t14, 1 nt after the ATG, without an annotated start codon and with phase 0, as a misannotated
    # cds_start_NF CDS can. Read from t14, the annotated stop codon at t40 is out of frame, and TGA at t26 is in frame.
    # TAA>CAA keeps the stop loss that the CDS shows.
    exon_seqs = [_STOP_TRANSCRIPT[:25], _STOP_TRANSCRIPT[25:]]
    row = run_pipeline_on_transcript(tmp_path, strand, exon_seqs, (14, 40), (40, "T", "C"), start_codon=False)

    assert row["stop_loss"] == True


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
    cds, vcf = _join_cds(), _join_vcf([("chr1", 150, 151, "v1")])
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
