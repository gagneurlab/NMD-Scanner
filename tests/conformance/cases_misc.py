"""
Conformance cases of what gives a row and what gives none (REF mismatch, intergenic and UTR variants, a
chromosome without a coding region, a VCF without a record), several transcripts or variants in one run, variant_id,
the location of the CDS in the transcript, reassign_exons, an Ensembl GFF3, and the values of unknown_reason. A
variant on a chromosome that the FASTA lacks is an error, and so is a CDS row outside the exon rows of its transcript.

The drawings show the transcript 5' to 3', also for the minus strand. The runner renders their layout block from the
case. The line ref is the transcript and the line alt is the transcript with the change. The two lines are aligned,
and `-` fills a gap. `[...]` is an exon and `|` is an exon junction. Upper case is the coding region (the CDS with its
stop codon), in codons of the annotated frame; lower case is UTR. An intron or a flank that the change touches is
drawn in lower case outside the brackets. `..N..` leaves out N bases. `^` marks the change as ref>alt. A further
variant in the VCF gets a `^` line of its own, with its ID; the alt line does not hold it. A mark line puts a
character under bases, e.g. `~~~` under a copy of the CDS, and `<-- label -->` spans a length. A ruler gives tx
positions. Positions are 0-based.
"""

import re

from .runner import (
    IDS,
    NO_PTC_FEATURES,
    NO_RULE,
    NOT_SCANNED,
    SAME_EXONS,
    UNKNOWN_ALT,
    Case,
    Change,
    Layout,
    Mark,
    NoRow,
    Raises,
    Ruler,
    Span,
    Transcript,
    names,
    per_strand,
)

# 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
#     <--5'UTR 8-->                        <--3'UTR 10->
#     tx 0          tx 14         tx 23                 tx 39
# The chromosome is 99 nt: flank 10, exon 1 at 10-24, intron 24-44, exon 2 at 44-53, intron 53-73, exon 3 at 73-89,
# flank 89-99. Both flanks have 10 nt, so the minus strand gives the same transcript and CDS bounds.
TX_SEQ = "GTCAGACCATGGCC" + "AGGCTGGGC" + "TCCTAAGCAGCCAGGC"
CDS_SEQ = "ATGGCC" + "AGGCTGGGC" + "TCCTAA"
THREE_EXONS = Layout(
    Transcript(("gtcagaccATGGCC", "AGGCTGGGC", "TCCTAAgcagccaggc")),
    {
        **IDS,
        "ref_cds_start": per_strand(18, 20),
        "ref_cds_stop": per_strand(79, 81),
        "ref_cds_seq": CDS_SEQ,
        "ref_cds_len": 21,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 6), (2, 9), (3, 6)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 18,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(18, "TAA")],
        "ref_stop_codon_exons": [3],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 89,
        "transcript_seq": TX_SEQ,
        "transcript_length": 39,
        "cds_start_in_transcript": 8,
        "cds_end_in_transcript": 29,
        "transcript_exon_info": [(1, 14), (2, 9), (3, 16)],
        "utr3_length": 10,
        "utr5_length": 8,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# CTG>CAG at CDS 10, in the middle of exon 2: a missense SNV
MISSENSE = Change("AGGC[T>A]GGGC")
MISSENSE_DRAWING = """
ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
alt 5' [gtcagacc ATG GCC]|[AGG CAG GGC]|[TCC TAA gcagccaggc] 3'
                                ^ T>A
T>A at CDS 10: CTG>CAG
"""
# The row of MISSENSE: the alt transcript is known, and no stop codon changes
MISSENSE_ROW = {
    "variant_id": "var1",
    "ref": per_strand("T", "A"),
    "alt": per_strand("A", "T"),
    "start_variant": per_strand(48, 50),
    "end_variant": per_strand(49, 51),
    "alt_cds_start": per_strand(18, 20),
    "alt_cds_stop": per_strand(79, 81),
    "alt_cds_seq": "ATGGCC" + "AGGCAGGGC" + "TCCTAA",
    "alt_cds_len": 21,
    "alt_cds_info": [(1, 6), (2, 9), (3, 6)],
    "alt_start_codon_pos": 0,
    "alt_start_codon_exon": 1,
    "alt_last_codon": "TAA",
    "alt_valid_stop": True,
    "alt_first_stop_codon": "TAA",
    "alt_first_stop_pos": 18,
    "alt_num_stop_codons": 1,
    "alt_all_stop_codons": [(18, "TAA")],
    "alt_stop_codon_exons": [3],
    "alt_is_premature": False,
    "start_loss": False,
    "stop_loss": False,
    "alt_transcript_seq": "GTCAGACCATGGCC" + "AGGCAGGGC" + "TCCTAAGCAGCCAGGC",
    "alt_transcript_length": 39,
    "alt_cds_start_in_transcript": 8,
    **NOT_SCANNED,
    "unknown_reason": None,
    **NO_PTC_FEATURES,
    "annotated_stop_distance": 0,
    **NO_RULE,
}


def _is_exon(row, number):
    """Whether a GFF3 row is the exon row with the exon_number attribute number."""
    return row[2] == "exon" and re.search(rf"exon_number={number}(;|$)", row[8]) is not None


def _move_3prime_end(row, strand, length):
    """Move the 3' end of a GFF3 row by length nt upstream, in transcript orientation."""
    if strand == "+":
        row[4] = str(int(row[4]) - length)
    else:
        row[3] = str(int(row[3]) + length)
    return row


def chromosome_2(rows, strand):
    """Put every GFF3 row on chr2. The FASTA and the VCF stay on chr1."""
    return [["chr2", *row[1:]] for row in rows]


def second_transcript_with_a_shorter_3utr(rows, strand):
    """
    Add the transcript tx2 to the gene g1: a copy of tx1 whose exon 3 ends 4 nt earlier, so its 3' UTR has 6 nt
    instead of 10.
    """
    copies = []
    for row in rows[1:]:
        copy = [*row[:8], row[8].replace("tx1", "tx2")]
        if copy[2] == "transcript" or _is_exon(copy, 3):
            copy = _move_3prime_end(copy, strand, 4)
        copies.append(copy)
    return rows + copies


def copy_of_the_gene_on_chromosome_2(rows, strand):
    """
    Add a copy of every GFF3 row on chr2, as the gene g2 with the transcript tx2. The FASTA and the VCF stay on chr1.
    """
    copies = [["chr2", *row[1:8], row[8].replace("tx1", "tx2").replace("=g1", "=g2")] for row in rows]
    return rows + copies


def exon_numbers_1_and_3_swapped(rows, strand):
    """Swap the exon_number attributes 1 and 3 on every row, so that they no longer follow the transcript order."""
    swap = {"1": "3", "3": "1"}
    return [
        [*row[:8], re.sub(r"exon_number=(\d+)", lambda m: "exon_number=" + swap.get(m[1], m[1]), row[8])]
        for row in rows
    ]


def exon_1_ends_2_nt_before_its_cds_row(rows, strand):
    """Move the 3' end of the exon 1 row 2 nt upstream. The CDS row of exon 1 then reaches 2 nt past the exon."""
    return [_move_3prime_end(row, strand, 2) if _is_exon(row, 1) else row for row in rows]


# MI-13: the 5' UTR holds the CDS sequence ATG GCC TAA too
#   5' [g atggcctaa cc ATG GCC TAA gc] 3'
#        <-copy--->    <--CDS--->
#   tx 0             tx 12      tx 21
# The chromosome is 43 nt: flank 10, the exon at 10-33, flank 33-43.
REPEATED_CDS = Layout(
    Transcript(("gatggcctaaccATGGCCTAAgc",)),
    {
        **IDS,
        "ref_cds_start": per_strand(22, 12),
        "ref_cds_stop": per_strand(31, 21),
        "ref_cds_seq": "ATGGCCTAA",
        "ref_cds_len": 9,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 9)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 6,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(6, "TAA")],
        "ref_stop_codon_exons": [1],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 33,
        "transcript_seq": "GATGGCCTAACCATGGCCTAAGC",
        "transcript_length": 23,
        "cds_start_in_transcript": 12,
        "cds_end_in_transcript": 21,
        "transcript_exon_info": [(1, 23)],
        "utr3_length": 2,
        "utr5_length": 12,
        "total_exon_count": 1,
        "likely_misannotated": False,
    },
)


# MI-19: Ensembl GFF3. The CDS rows of exon 2 and exon 3 have phase 2, because exon 1 holds 7 CDS nt and exon 2 9.
# has_start_codon and has_stop_codon come from the FASTA: ATG at CDS 0 with phase 0, and TAA as the last 3 CDS nt.
#   5' [gtcagacc ATG GCC A]|[GG CTG GGC T]|[CC TAA gcagccaggc] 3'
ENSEMBL = Layout(
    Transcript(("gtcagaccATGGCCA", "GGCTGGGCT", "CCTAAgcagccaggc"), flavor="ensembl"),
    {
        **THREE_EXONS.ref,
        "ref_cds_info": [(1, 7), (2, 9), (3, 5)],
        "transcript_exon_info": [(1, 15), (2, 9), (3, 15)],
    },
)

CASES = [
    # MI-02, MI-18: a VCF whose only record has a REF that the genome does not have
    Case(
        "only_variant_with_a_ref_mismatch_gives_no_row",
        """
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GTC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
                              ^ C>T
        C>T at CDS 4: the VCF has REF A, the genome C.
        """,
        THREE_EXONS,
        Change("CCATGG[C>T]CGTAAG", vcf_ref="A"),
        NoRow("REF mismatch"),
    ),
    # MI-01: the variant with the REF mismatch gives no row, the other variant keeps its row
    Case(
        "variant_with_a_ref_mismatch_gives_no_row_and_the_other_variant_keeps_its_row",
        """
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC]|[AGG CAG GGC]|[TCC TAA gcagccaggc] 3'
                                        ^ T>A
                              ^ C>T var2
        var1 T>A at CDS 10: the row.
        var2 C>T at CDS 4: the VCF has REF A, the genome C: no row.
        """,
        THREE_EXONS,
        MISSENSE,
        {**MISSENSE_ROW, "alt_transcript_exon_info": SAME_EXONS, "nmd_model_status": "no_ptc"},
        more_changes=(Change("CCATGG[C>T]CGTAAG", vcf_ref="A", vcf_id="var2"),),
    ),
    # MI-03, MI-18
    Case(
        "intergenic_snv_gives_no_row",
        """
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc]cccccccccc 3'
        alt 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc]ccccaccccc 3'
                                                                        ^ C>A
        C>A in the 3' flank, the 5th base after the transcript end.
        """,
        THREE_EXONS,
        Change("CAGGCCCCC[C>A]CCCCC"),
        NoRow("touches no coding region"),
    ),
    # MI-04
    Case(
        "snv_in_the_5utr_away_from_the_start_codon_gives_no_row",
        """
        tx       1       8
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gacagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
                 ^ T>A
        T>A at tx 1, 7 nt before the start codon.
        """,
        THREE_EXONS,
        Change("CG[T>A]CAGACC"),
        NoRow("touches no coding region"),
        ruler=Ruler((1, 8)),
    ),
    # MI-04
    Case(
        "snv_in_the_3utr_away_from_the_stop_codon_gives_no_row",
        """
        tx                                               29    35
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagcctggc] 3'
                                                               ^ A>T
        A>T at tx 35, 6 nt after the stop codon.
        """,
        THREE_EXONS,
        Change("CAGCC[A>T]GGC"),
        NoRow("touches no coding region"),
        ruler=Ruler((29, 35)),
    ),
    # MI-09: the variant lies next to the coding region but changes only the UTR
    Case(
        "snv_right_before_the_start_codon_gives_no_row",
        """
        tx             7
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagaca ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
                       ^ C>A
        C>A at tx 7, the last base of the 5' UTR.
        """,
        THREE_EXONS,
        Change("AGAC[C>A]ATGGCC"),
        NoRow("touches no coding region"),
        ruler=Ruler((7,)),
    ),
    # MI-09
    Case(
        "snv_right_after_the_stop_codon_gives_no_row",
        """
        tx                                               29
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA ccagccaggc] 3'
                                                         ^ G>C
        G>C at tx 29, the first base of the 3' UTR.
        """,
        THREE_EXONS,
        Change("CTAA[G>C]CAGCC"),
        NoRow("touches no coding region"),
        ruler=Ruler((29,)),
    ),
    # MI-05
    Case(
        "snv_on_a_chromosome_without_a_coding_region_gives_no_row",
        """
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC]|[AGG CAG GGC]|[TCC TAA gcagccaggc] 3'
                                        ^ T>A
        The GFF3 has the transcript on chr2. The VCF and the FASTA are on chr1, which has no coding region.
        The variant T>A lies on chr1 where chr2 has CDS 10.
        """,
        Layout(Transcript(THREE_EXONS.transcript.exons, edit_gff3=chromosome_2), THREE_EXONS.ref),
        MISSENSE,
        NoRow("touches no coding region"),
    ),
    # MI-08, MI-18
    Case(
        "vcf_without_a_record_gives_no_row",
        """
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        The VCF has its header and no record.
        """,
        THREE_EXONS,
        None,
        NoRow("no record"),
    ),
    # MI-06: one row per transcript; each row has the transcript columns of its own transcript
    Case(
        "variant_in_the_cds_of_two_transcripts_gives_a_row_for_each",
        """
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC]|[AGG CAG GGC]|[TCC TAA gcagccaggc] 3'
                                        ^ T>A
                                                         <--------> tx1: utr3_length = 10
                                                         <----> tx2: utr3_length = 6
        The block shows tx1. tx2 is a copy of tx1 whose exon 3 ends 4 nt earlier, so its 3' UTR is gcagcc.
        T>A at CDS 10: CTG>CAG in both.
        """,
        Layout(
            Transcript(THREE_EXONS.transcript.exons, edit_gff3=second_transcript_with_a_shorter_3utr), THREE_EXONS.ref
        ),
        MISSENSE,
        {**MISSENSE_ROW, "alt_transcript_exon_info": SAME_EXONS, "nmd_model_status": "no_ptc"},
        more_rows=(
            {
                "transcript_id": "tx2",
                "transcript_start": per_strand(10, 14),
                "transcript_end": per_strand(85, 89),
                "transcript_seq": "GTCAGACCATGGCC" + "AGGCTGGGC" + "TCCTAAGCAGCC",
                "transcript_length": 35,
                "transcript_exon_info": [(1, 14), (2, 9), (3, 12)],
                "utr3_length": 6,
                "alt_transcript_seq": "GTCAGACCATGGCC" + "AGGCAGGGC" + "TCCTAAGCAGCC",
                "alt_transcript_length": 35,
            },
        ),
        marks=(Span("ref", 29, 39, "tx1: utr3_length = 10"), Span("ref", 29, 35, "tx2: utr3_length = 6")),
    ),
    # Pins the row rule of "Output columns" in a run with a touched and an untouched transcript: a row is a variant
    # in a transcript whose coding region it touches. So tx1 gives its row and tx2 gives none. In the closest cases,
    # every transcript is touched (MI-06) or none is (MI-05). The FASTA lacks chr2. That is no error, because chr2
    # has no variant.
    Case(
        "variant_in_tx1_with_an_untouched_copy_of_the_gene_on_chr2_gives_one_row",
        """
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC]|[AGG CAG GGC]|[TCC TAA gcagccaggc] 3'
                                        ^ T>A
        The block shows tx1 on chr1. The GFF3 also has every row of tx1 and its gene on chr2, as the gene g2 with
        the transcript tx2. The VCF and the FASTA have chr1 only.
        T>A at CDS 10: CTG>CAG in tx1. No variant touches the coding region of tx2, so tx2 gives no row.
        """,
        Layout(Transcript(THREE_EXONS.transcript.exons, edit_gff3=copy_of_the_gene_on_chromosome_2), THREE_EXONS.ref),
        MISSENSE,
        {**MISSENSE_ROW, "alt_transcript_exon_info": SAME_EXONS, "nmd_model_status": "no_ptc"},
    ),
    # MI-07: one row per variant; each row applies only its own variant
    Case(
        "two_variants_in_one_transcript_give_a_row_each_with_only_their_own_change",
        """
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC]|[AGG CAG GGC]|[TCC TAA gcagccaggc] 3'
                                        ^ T>A
                                            ^ G>A var2
        var1 T>A at CDS 10: CTG>CAG. var2 G>A at CDS 13: GGC>GAC.
        The row of each variant applies only its own change.
        """,
        THREE_EXONS,
        MISSENSE,
        {**MISSENSE_ROW, "alt_transcript_exon_info": SAME_EXONS, "nmd_model_status": "no_ptc"},
        more_changes=(Change("CTGG[G>A]CGTAAG", vcf_id="var2"),),
        more_rows=(
            {
                "variant_id": "var2",
                "ref": per_strand("G", "C"),
                "alt": per_strand("A", "T"),
                "start_variant": per_strand(51, 47),
                "end_variant": per_strand(52, 48),
                "alt_cds_seq": "ATGGCC" + "AGGCTGGAC" + "TCCTAA",
                "alt_transcript_seq": "GTCAGACCATGGCC" + "AGGCTGGAC" + "TCCTAAGCAGCCAGGC",
            },
        ),
    ),
    # MI-10, UR-01: a record without ID
    Case(
        "variant_without_an_id_gets_the_variant_id_dot",
        MISSENSE_DRAWING + "The VCF record has the ID '.'.\n",
        THREE_EXONS,
        Change("AGGC[T>A]GGGC", vcf_id="."),
        {
            **MISSENSE_ROW,
            "variant_id": ".",
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
    ),
    # MI-10
    Case(
        "numeric_variant_id_keeps_its_leading_zeros",
        MISSENSE_DRAWING + "The VCF record has the ID 007.\n",
        THREE_EXONS,
        Change("AGGC[T>A]GGGC", vcf_id="007"),
        {
            **MISSENSE_ROW,
            "variant_id": "007",
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
    ),
    # MI-10
    Case(
        "variant_id_na_stays_text",
        MISSENSE_DRAWING + "The VCF record has the ID NA.\n",
        THREE_EXONS,
        Change("AGGC[T>A]GGGC", vcf_id="NA"),
        {
            **MISSENSE_ROW,
            "variant_id": "NA",
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
    ),
    # MI-13: the CDS is located by its coordinates, at tx 12, not at the copy at tx 1
    Case(
        "cds_sequence_repeated_in_the_5utr_is_located_by_its_coordinates",
        """
        tx       1           12
        ref 5' [gatggcctaacc ATG GCC TAA gc] 3'
        alt 5' [gatggcctaacc ATG GAC TAA gc] 3'
                                  ^ C>A
                 ~~~~~~~~~ copy of the CDS
        C>A at CDS 4: GCC>GAC. The CDS starts at tx 12, by its coordinates. The copy at tx 1 in the 5' UTR stays.
        """,
        REPEATED_CDS,
        Change("CCATGG[C>A]C"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(26, 16),
            "end_variant": per_strand(27, 17),
            "alt_cds_start": per_strand(22, 12),
            "alt_cds_stop": per_strand(31, 21),
            "alt_cds_seq": "ATGGACTAA",
            "alt_cds_len": 9,
            "alt_cds_info": [(1, 9)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 6,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(6, "TAA")],
            "alt_stop_codon_exons": [1],
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GATGGCCTAACCATGGACTAAGC",
            "alt_transcript_length": 23,
            "alt_cds_start_in_transcript": 12,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        marks=(Mark("ref", 1, 10, "~", "copy of the CDS"),),
        ruler=Ruler((1, 12)),
    ),
    # MI-15: reassign_exons numbers the exons by their position, 5' to 3'. The row is that of the correct numbers.
    Case(
        "reassign_exons_numbers_the_exons_by_their_position",
        MISSENSE_DRAWING
        + "The GFF3 has exon_number 3 on exon 1 and 1 on exon 3. annotate() runs with reassign_exons.\n",
        Layout(Transcript(THREE_EXONS.transcript.exons, edit_gff3=exon_numbers_1_and_3_swapped), THREE_EXONS.ref),
        MISSENSE,
        {**MISSENSE_ROW, "alt_transcript_exon_info": SAME_EXONS, "nmd_model_status": "no_ptc"},
        reassign_exons=True,
    ),
    # MI-16: README.md, section Arguments: the chromosome names must match in the VCF, the GFF3 and the FASTA
    Case(
        "variant_on_a_chromosome_that_the_fasta_lacks_is_an_error_that_names_the_chromosome",
        MISSENSE_DRAWING
        + "The GFF3 and the VCF are on chr1, the FASTA has chr2 only. annotate() rejects the input with a ValueError\n"
        + "that names the chromosome chr1.\n",
        Layout(Transcript(THREE_EXONS.transcript.exons, fasta_contig="chr2"), THREE_EXONS.ref),
        MISSENSE,
        Raises(ValueError, match=r"\bchromosome\(s\) with variants and CDS rows: chr1\b"),
    ),
    # MI-17, PF-25 (CDS rows partly outside the exons), SC-18
    Case(
        "cds_row_that_reaches_past_its_exon_end_is_an_error_that_names_the_transcript_and_the_cds_row",
        """
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC]|[AGG CAG GGC]|[TCC TAA gcagccaggc] 3'
                                        ^ T>A
                              xx past the end of the exon 1 row
        The block shows the CDS rows. The exon 1 row ends 2 nt before the CDS row of exon 1, so the CDS row
        reaches 2 nt into the intron. T>A at CDS 10.
        annotate() rejects the GFF3 with a ValueError that names the transcript tx1 and the CDS row.
        """,
        Layout(
            Transcript(THREE_EXONS.transcript.exons, edit_gff3=exon_1_ends_2_nt_before_its_cds_row), THREE_EXONS.ref
        ),
        MISSENSE,
        Raises(
            ValueError,
            match=per_strand(
                names("CDS row of transcript tx1 at chr1:19-24"), names("CDS row of transcript tx1 at chr1:76-81")
            ),
        ),
        marks=(Mark("ref", 12, 14, "x", "past the end of the exon 1 row"),),
    ),
    # MI-19: an Ensembl GFF3 with phase 2 CDS rows
    Case(
        "ensembl_gff3_with_phase_2_cds_rows",
        """
        ref 5' [gtcagacc ATG GCC A]|[GG CTG GGC T]|[CC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC A]|[GG CAG GGC T]|[CC TAA gcagccaggc] 3'
                                         ^ T>A
        The CDS rows of exon 1, 2 and 3 have phase 0, 2 and 2.
        T>A at CDS 10: CTG>CAG.
        """,
        ENSEMBL,
        Change("GGC[T>A]GGGCT"),
        {
            **MISSENSE_ROW,
            "alt_cds_info": [(1, 7), (2, 9), (3, 5)],
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
    ),
    # UR-02: no placement keeps the donor GT after exon 1
    Case(
        "snv_in_the_donor_gt_gives_splice_site_destroyed",
        """
        ref 5' [gtcagacc ATG GCC]gtaagtcccccccctttcag[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC]ataagtcccccccctttcag[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
                                 ^ G>A
        G>A at intron 1 +1, the G of the donor GT.
        """,
        THREE_EXONS,
        Change("CC[G>A]TAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(24, 74),
            "end_variant": per_strand(25, 75),
            **UNKNOWN_ALT,
            "unknown_reason": "splice_site_destroyed",
        },
    ),
    # UR-04: AG inserted at CAG|AG. One placement keeps the exon start at the old AG, another puts the inserted AG
    # into the exon.
    Case(
        "ag_inserted_at_cag_ag_gives_exon_boundary_ambiguous",
        """
        ref 5' [gtcagacc ATG GCC]gtaagtcccccccctttcag[--AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC]gtaagtcccccccctttcag[AGAGG CTG GGC]|[TCC TAA gcagccaggc] 3'
                                                      ^^ ->AG
        AG inserted between the intron 1 end tttcag and the exon 2 start AGG. One placement keeps the exon
        start at the old AG: ...tttcagag[AGG... Another puts the inserted AG into the exon: ...tttcag[AGAGG...
        The alt line shows the second one.
        """,
        THREE_EXONS,
        Change("TTTCAG[>AG]AGGCTG"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "T"),
            "alt": per_strand("GAG", "TCT"),
            "start_variant": per_strand(43, 54),
            "end_variant": per_strand(44, 55),
            **UNKNOWN_ALT,
            "unknown_reason": "exon_boundary_ambiguous",
        },
        equivalent=(Change("TTTC[>AG]AGAGGC"), Change("CAGAG[>AG]GCTG")),
    ),
    # UR-06: an MNV has one placement, so it is never exon_boundary_ambiguous, also where its ALT forms an AG one
    # base into the exon
    Case(
        "mnv_at_cag_ag_whose_alt_forms_another_ag_is_not_ambiguous",
        """
        ref 5' [gtcagacc ATG GCC]|[AGG CTG GGC]|[TCC TAA gcagccaggc] 3'
        alt 5' [gtcagacc ATG GCC]|[GAG CTG GGC]|[TCC TAA gcagccaggc] 3'
                                   ^^ AG>GA
        AG>GA at CDS 6 and 7: AGG>GAG, at the start of exon 2. Intron 1 ends in tttcag. The alt GAG forms
        another AG one base into the exon.
        """,
        THREE_EXONS,
        Change("TTTCAG[AG>GA]GCTG"),
        {
            **MISSENSE_ROW,
            "ref": per_strand("AG", "CT"),
            "alt": per_strand("GA", "TC"),
            "start_variant": per_strand(44, 53),
            "end_variant": per_strand(46, 55),
            "alt_cds_seq": "ATGGCC" + "GAGCTGGGC" + "TCCTAA",
            "alt_transcript_seq": "GTCAGACCATGGCC" + "GAGCTGGGC" + "TCCTAAGCAGCCAGGC",
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("TTTCAG[AGG>GAG]CTG"),),
    ),
]
