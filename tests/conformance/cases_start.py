"""
Conformance cases of the start codon and of the reading frame.

The start codon cases pin start_loss, the scan of the alt transcript for the next ATG after a start loss, the
classification of the rescued ORF, non-ATG start codons, and the start codon of an Ensembl GFF3, which comes from the
FASTA. The frame cases pin cds_frame, the GFF3 phase of the 5'-most CDS row, and the codon scans that start at the
first complete codon. A 5'-most CDS row with the phase "." is an error.

The drawings show the transcript 5' to 3', also for the minus strand. The runner renders their layout block from the
case. The line ref is the transcript and the line alt is the transcript with the change. The two lines are aligned,
and `-` fills a gap. `[...]` is an exon and `|` is an exon junction. Upper case is the coding region (the CDS with its
stop codon), in codons of the annotated frame; lower case is UTR. `..N..` leaves out N bases. `^` marks the change as
ref>alt. A mark line puts a character under bases: `***` under the PTC, `aaa` under the ATG of the scan and `sss`
under a stop codon. `<-- label -->` spans a length. A ruler gives tx, alt tx or CDS positions, as labelled.
"""

from .runner import (
    IDS,
    NO_PTC_FEATURES,
    NO_RULE,
    NOT_SCANNED,
    SAME_EXONS,
    Case,
    Change,
    Layout,
    Mark,
    Raises,
    Ruler,
    Span,
    Transcript,
    per_strand,
)

ALT_CODON_COLUMNS = (
    "alt_start_codon_pos",
    "alt_start_codon_exon",
    "alt_last_codon",
    "alt_valid_stop",
    "alt_first_stop_codon",
    "alt_first_stop_pos",
    "alt_num_stop_codons",
    "alt_all_stop_codons",
    "alt_stop_codon_exons",
)
SCAN_COLUMNS = tuple(NOT_SCANNED)
PTC_FEATURE_COLUMNS = tuple(NO_PTC_FEATURES)
RULE_COLUMNS = tuple(NO_RULE)


def variant(ref, alt, start, end):
    """The columns of the VCF record var1: ref, alt, start_variant and end_variant."""
    return {"variant_id": "var1", "ref": ref, "alt": alt, "start_variant": start, "end_variant": end}


def columns(names, *values):
    """The expected values of the columns names, in their order."""
    return dict(zip(names, values, strict=True))


# A CDS of 15 nt over two exons: exon 1 holds the 3 nt 5' UTR and the first 5 CDS bases, exon 2 the rest of the CDS
# and the 5 nt 3' UTR. On the minus strand, the layout of 63 nt is reverse complemented.
ATG_START = Layout(
    Transcript(("gggATGAA", "ACCCGACTAAggggg")),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 15),
        "ref_cds_stop": per_strand(48, 50),
        "ref_cds_seq": "ATGAAACCCGACTAA",
        "ref_cds_len": 15,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 5), (2, 10)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 12,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(12, "TAA")],
        "ref_stop_codon_exons": [2],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 53,
        "transcript_seq": "GGGATGAAACCCGACTAAGGGGG",
        "transcript_length": 23,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 18,
        "transcript_exon_info": [(1, 8), (2, 15)],
        "utr3_length": 5,
        "utr5_length": 3,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)

# As ATG_START, with the non-ATG start codon CTG, which the start_codon rows mark
CTG_START = Layout(
    Transcript(("gggCTGAA", "ACCCGACTAAggggg")),
    {
        **ATG_START.ref,
        "ref_cds_seq": "CTGAAACCCGACTAA",
        "transcript_seq": "GGGCTGAAACCCGACTAAGGGGG",
    },
)

# The 9 scan columns after a start loss whose scan finds no ATG: the scan reads no stop codon
NO_ATG_FOUND = {
    "transcript_start_codon_pos": None,
    "transcript_start_codon_exon": None,
    "transcript_last_codon": "GGG",
    "transcript_valid_stop": False,
    "transcript_first_stop_codon": None,
    "transcript_first_stop_pos": None,
    "transcript_num_stop_codons": 0,
    "transcript_all_stop_codons": [],
    "transcript_stop_codon_exons": [],
}

CASES = [
    # SC-01, SC-11, SC-22 (has_start_codon from the start_codon rows)
    Case(
        "atg_to_acg_is_a_start_loss_and_without_a_next_atg_neither_a_ptc_nor_a_stop_loss",
        """
        tx      0   3        8         15  18   23
        ref 5' [ggg ATG AA]|[A CCC GAC TAA ggggg] 3'
        alt 5' [ggg ACG AA]|[A CCC GAC TAA ggggg] 3'
                     ^ T>C
        ATG>ACG changes the start codon. From tx 3 on, the alt transcript has no ATG, so the scan reads no stop codon.
        No ORF overlaps the CDS: stop_codon_distance = null.
        """,
        ATG_START,
        Change("gggA[T>C]GAA"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "start_variant": per_strand(14, 48),
            "end_variant": per_strand(15, 49),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(48, 50),
            "alt_cds_seq": "ACGAAACCCGACTAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 5), (2, 10)],
            "alt_start_codon_pos": None,
            "alt_start_codon_exon": None,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 12,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(12, "TAA")],
            "alt_stop_codon_exons": [2],
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGACGAAACCCGACTAAGGGGG",
            "alt_transcript_length": 23,
            "alt_cds_start_in_transcript": 3,
            **NO_ATG_FOUND,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        ruler=Ruler((0, 3, 8, 15, 18, 23)),
    ),
    # SC-02
    Case(
        "ctg_to_ccg_is_a_start_loss_of_a_non_atg_start_codon",
        """
        tx      0   3        8         15  18   23
        ref 5' [ggg CTG AA]|[A CCC GAC TAA ggggg] 3'
        alt 5' [ggg CCG AA]|[A CCC GAC TAA ggggg] 3'
                     ^ T>C
        The start_codon rows mark CTG as the start codon. CTG>CCG changes it. The scan finds no ATG.
        """,
        CTG_START,
        Change("gggC[T>C]GAA"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "start_variant": per_strand(14, 48),
            "end_variant": per_strand(15, 49),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(48, 50),
            "alt_cds_seq": "CCGAAACCCGACTAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 5), (2, 10)],
            "alt_start_codon_pos": None,
            "alt_start_codon_exon": None,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 12,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(12, "TAA")],
            "alt_stop_codon_exons": [2],
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGCCGAAACCCGACTAAGGGGG",
            "alt_transcript_length": 23,
            "alt_cds_start_in_transcript": 3,
            **NO_ATG_FOUND,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        ruler=Ruler((0, 3, 8, 15, 18, 23)),
    ),
    # SC-03
    Case(
        "missense_after_a_ctg_start_codon_keeps_the_start_codon",
        """
        tx      0   3        8  10     15  18   23
        ref 5' [ggg CTG AA]|[A CCC GAC TAA ggggg] 3'
        alt 5' [ggg CTG AA]|[A CGC GAC TAA ggggg] 3'
                                ^ C>G
        CCC>CGC is a missense. The alt CDS starts with the annotated start codon CTG too: alt_start_codon_pos = 0.
        """,
        CTG_START,
        Change("AC[C>G]CGAC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("G", "C"),
            "start_variant": per_strand(40, 22),
            "end_variant": per_strand(41, 23),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(48, 50),
            "alt_cds_seq": "CTGAAACGCGACTAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 5), (2, 10)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 12,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(12, "TAA")],
            "alt_stop_codon_exons": [2],
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GGGCTGAAACGCGACTAAGGGGG",
            "alt_transcript_length": 23,
            "alt_cds_start_in_transcript": 3,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        ruler=Ruler((0, 3, 8, 10, 15, 18, 23)),
    ),
]

# SC-04: the CDS starts with the annotated start codon CTG and has an in-frame ATG at CDS 30
INTERNAL_MET_CDS = "CTG" + "AAA" * 9 + "ATG" + "AAA" * 42 + "TGGGACTAA"
INTERNAL_MET = Layout(
    Transcript(("ggg" + INTERNAL_MET_CDS[:100], INTERNAL_MET_CDS[100:] + "ggggg")),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 15),
        "ref_cds_stop": per_strand(201, 203),
        "ref_cds_seq": INTERNAL_MET_CDS,
        "ref_cds_len": 168,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 100), (2, 68)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 165,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(165, "TAA")],
        "ref_stop_codon_exons": [2],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 206,
        "transcript_seq": "GGG" + INTERNAL_MET_CDS + "GGGGG",
        "transcript_length": 176,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 171,
        "transcript_exon_info": [(1, 103), (2, 73)],
        "utr3_length": 5,
        "utr5_length": 3,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)

# SC-20: the start codon CTG, an in-frame ATG at CDS 6, and a TAG in the 3' UTR in the frame of the CDS
CTG_THEN_ATG = Layout(
    Transcript(("ggCTGAA", "AATGCCCTAAgggtagcc")),
    {
        **IDS,
        "ref_cds_start": per_strand(12, 18),
        "ref_cds_stop": per_strand(47, 53),
        "ref_cds_seq": "CTGAAAATGCCCTAA",
        "ref_cds_len": 15,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 5), (2, 10)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 12,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(12, "TAA")],
        "ref_stop_codon_exons": [2],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 55,
        "transcript_seq": "GGCTGAAAATGCCCTAAGGGTAGCC",
        "transcript_length": 25,
        "cds_start_in_transcript": 2,
        "cds_end_in_transcript": 17,
        "transcript_exon_info": [(1, 7), (2, 18)],
        "utr3_length": 8,
        "utr5_length": 2,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)


# The REF_COLUMNS that the layouts with the 5' UTR ggg, an annotated ATG, one in-frame stop codon, the annotated TAA, and
# the 3' UTR ggggg share
ATG_TAA_BASE = {
    **IDS,
    "ref_cds_start": per_strand(13, 15),
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "cds_in_transcript": True,
    "ref_start_codon_pos": 0,
    "ref_start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_num_stop_codons": 1,
    "ref_is_premature": False,
    "transcript_start": 10,
    "cds_start_in_transcript": 3,
    "utr3_length": 5,
    "utr5_length": 3,
    "likely_misannotated": False,
}

# SC-16: an ATG at CDS 6, in frame, and a CAG at CDS 9
MET_RESCUE_CDS = "ATGGCCATGCAG" + "AAA" * 36 + "GACTAA"
MET_RESCUE = Layout(
    Transcript(("ggg" + MET_RESCUE_CDS[:60], MET_RESCUE_CDS[60:120], MET_RESCUE_CDS[120:] + "ggggg")),
    {
        **ATG_TAA_BASE,
        "ref_cds_stop": per_strand(179, 181),
        "ref_cds_seq": MET_RESCUE_CDS,
        "ref_cds_len": 126,
        "ref_cds_info": [(1, 60), (2, 60), (3, 6)],
        "ref_first_stop_pos": 123,
        "ref_all_stop_codons": [(123, "TAA")],
        "ref_stop_codon_exons": [3],
        "transcript_end": 184,
        "transcript_seq": "GGG" + MET_RESCUE_CDS + "GGGGG",
        "transcript_length": 134,
        "cds_end_in_transcript": 129,
        "transcript_exon_info": [(1, 63), (2, 60), (3, 11)],
        "total_exon_count": 3,
    },
)

# SC-08: an ATG at CDS 4, out of frame, whose frame reads CTA AAA as TAA
OUT_OF_FRAME_RESCUE_CDS = "ATGGATGCCCTA" + "AAA" * 50 + "GACTAA"
OUT_OF_FRAME_RESCUE = Layout(
    Transcript(
        ("ggg" + OUT_OF_FRAME_RESCUE_CDS[:60], OUT_OF_FRAME_RESCUE_CDS[60:120], OUT_OF_FRAME_RESCUE_CDS[120:] + "ggggg")
    ),
    {
        **ATG_TAA_BASE,
        "ref_cds_stop": per_strand(221, 223),
        "ref_cds_seq": OUT_OF_FRAME_RESCUE_CDS,
        "ref_cds_len": 168,
        "ref_cds_info": [(1, 60), (2, 60), (3, 48)],
        "ref_first_stop_pos": 165,
        "ref_all_stop_codons": [(165, "TAA")],
        "ref_stop_codon_exons": [3],
        "transcript_end": 226,
        "transcript_seq": "GGG" + OUT_OF_FRAME_RESCUE_CDS + "GGGGG",
        "transcript_length": 176,
        "cds_end_in_transcript": 171,
        "transcript_exon_info": [(1, 63), (2, 60), (3, 53)],
        "total_exon_count": 3,
    },
)

# SC-05, SC-09: GTA AGC after the start codon, and an ATG at CDS 9 in frame
FRAMESHIFT_RESCUE = Layout(
    Transcript(("gggATGGT", "AAGCATGGCCAAAGACTAAggggg")),
    {
        **ATG_TAA_BASE,
        "ref_cds_stop": per_strand(57, 59),
        "ref_cds_seq": "ATGGTAAGCATGGCCAAAGACTAA",
        "ref_cds_len": 24,
        "ref_cds_info": [(1, 5), (2, 19)],
        "ref_first_stop_pos": 21,
        "ref_all_stop_codons": [(21, "TAA")],
        "ref_stop_codon_exons": [2],
        "transcript_end": 62,
        "transcript_seq": "GGGATGGTAAGCATGGCCAAAGACTAAGGGGG",
        "transcript_length": 32,
        "cds_end_in_transcript": 27,
        "transcript_exon_info": [(1, 8), (2, 24)],
        "total_exon_count": 2,
    },
)

# The REF_COLUMNS of the two exon layouts ggg ATG CC | C AAA GAC TAA, with a 3' UTR of their own
ATG_CCC_AAA_GAC = {
    **ATG_TAA_BASE,
    "ref_cds_seq": "ATGCCCAAAGACTAA",
    "ref_cds_len": 15,
    "ref_cds_info": [(1, 5), (2, 10)],
    "ref_first_stop_pos": 12,
    "ref_all_stop_codons": [(12, "TAA")],
    "ref_stop_codon_exons": [2],
    "cds_end_in_transcript": 18,
    "total_exon_count": 2,
}

# SC-12: an ATG in the 3' UTR, 5 nt after the first base of the stop codon
ATG_IN_THE_3UTR = Layout(
    Transcript(("gggATGCC", "CAAAGACTAAggatgccctgagg")),
    {
        **ATG_CCC_AAA_GAC,
        "ref_cds_start": per_strand(13, 23),
        "ref_cds_stop": per_strand(48, 58),
        "transcript_end": 61,
        "transcript_seq": "GGGATGCCCAAAGACTAAGGATGCCCTGAGG",
        "transcript_length": 31,
        "transcript_exon_info": [(1, 8), (2, 23)],
        "utr3_length": 13,
    },
)

# SC-13: an ATG that starts at the last base of the stop codon: TAA tg
ATG_AT_THE_LAST_STOP_CODON_BASE = Layout(
    Transcript(("gggATGCC", "CAAAGACTAAtgccctgagg")),
    {
        **ATG_CCC_AAA_GAC,
        "ref_cds_start": per_strand(13, 20),
        "ref_cds_stop": per_strand(48, 55),
        "transcript_end": 58,
        "transcript_seq": "GGGATGCCCAAAGACTAATGCCCTGAGG",
        "transcript_length": 28,
        "transcript_exon_info": [(1, 8), (2, 20)],
        "utr3_length": 10,
    },
)

# SC-13: an ATG that starts 1 nt before the stop codon TGA: CCA TGA
ATG_BEFORE_THE_STOP_CODON = Layout(
    Transcript(("gggATGCC", "CAAACCATGAcctaggg")),
    {
        **ATG_CCC_AAA_GAC,
        "ref_cds_start": per_strand(13, 17),
        "ref_cds_stop": per_strand(48, 52),
        "ref_cds_seq": "ATGCCCAAACCATGA",
        "ref_last_codon": "TGA",
        "ref_first_stop_codon": "TGA",
        "ref_all_stop_codons": [(12, "TGA")],
        "transcript_end": 55,
        "transcript_seq": "GGGATGCCCAAACCATGACCTAGGG",
        "transcript_length": 25,
        "transcript_exon_info": [(1, 8), (2, 17)],
        "utr3_length": 7,
    },
)

# SC-10: an ATG at CDS 4, out of frame, whose frame reads on past the stop codon to a TGA in the 3' UTR
RESCUE_PAST_THE_STOP = Layout(
    Transcript(("gggATGGATGC", "CCCCAAAGACTAAgtgacc")),
    {
        **ATG_TAA_BASE,
        "ref_cds_start": per_strand(13, 16),
        "ref_cds_stop": per_strand(54, 57),
        "ref_cds_seq": "ATGGATGCCCCCAAAGACTAA",
        "ref_cds_len": 21,
        "ref_cds_info": [(1, 8), (2, 13)],
        "ref_first_stop_pos": 18,
        "ref_all_stop_codons": [(18, "TAA")],
        "ref_stop_codon_exons": [2],
        "transcript_end": 60,
        "transcript_seq": "GGGATGGATGCCCCCAAAGACTAAGTGACC",
        "transcript_length": 30,
        "cds_end_in_transcript": 24,
        "transcript_exon_info": [(1, 11), (2, 19)],
        "utr3_length": 6,
        "total_exon_count": 2,
    },
)

CASES += [
    # SC-04
    Case(
        "ptc_after_a_ctg_start_codon_is_measured_from_the_ctg_not_from_the_internal_met",
        """
        CDS         0                  30                                           159     165
        ref 5' [ggg CTG AAA ..21.. AAA ATG AAA ..60.. AAA A]|[AA AAA ..48.. AAA AAA TGG GAC TAA ggggg] 3'
        alt 5' [ggg CTG AAA ..21.. AAA ATG AAA ..60.. AAA A]|[AA AAA ..48.. AAA AAA TAG GAC TAA ggggg] 3'
                                                                                     ^ G>A
                                                                                    *** PTC
                    <------------ ptc_to_start_codon = 159, not < 150 ------------>
                                                                                    <---------------> ptc_to_intron = 14
        exon 1: 103 nt, exon 2: 73 nt
        TGG>TAG is a PTC in the last exon. The ATG at CDS 30 is an internal Met: from it, the PTC would be 129 nt away.
        """,
        INTERNAL_MET,
        Change("AAT[G>A]GGAC"),
        {
            **variant(per_strand("G", "C"), per_strand("A", "T"), per_strand(193, 22), per_strand(194, 23)),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(201, 203),
            "alt_cds_seq": "CTG" + "AAA" * 9 + "ATG" + "AAA" * 42 + "TAGGACTAA",
            "alt_cds_len": 168,
            "alt_cds_info": [(1, 100), (2, 68)],
            **columns(ALT_CODON_COLUMNS, 0, 1, "TAA", True, "TAG", 159, 2, [(159, "TAG"), (165, "TAA")], [2, 2]),
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GGGCTG" + "AAA" * 9 + "ATG" + "AAA" * 42 + "TAGGACTAAGGGGG",
            "alt_transcript_length": 176,
            "alt_cds_start_in_transcript": 3,
            **NOT_SCANNED,
            "unknown_reason": None,
            **columns(PTC_FEATURE_COLUMNS, 1, 0, 159, False, 73, 14),
            "stop_codon_distance": 6,
            **columns(RULE_COLUMNS, True, False, False, False, False, True),
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 162, 165, "*", "PTC"),
            Span("alt", 3, 162, "ptc_to_start_codon = 159, not < 150"),
            Span("alt", 162, 176, "ptc_to_intron = 14"),
        ),
        ruler=Ruler((0, 30, 159, 165), "CDS"),
    ),
    # SC-20
    Case(
        "stop_loss_scan_starts_at_the_annotated_ctg_not_at_the_in_frame_atg",
        """
        tx      0  2        7 8       14     20   25
        ref 5' [gg CTG AA]|[A ATG CCC TAA gggtagcc] 3'
        alt 5' [gg CTG AA]|[A ATG CCC CAA gggtagcc] 3'
                                      ^ T>C
        TAA>CAA loses the stop codon. The scan reads on in the frame of the CDS to the TAG at tx 20. Its start codon
        is the annotated CTG at tx 2, not the in-frame ATG at tx 8. stop_codon_distance = 14 - 20 = -6.
        """,
        CTG_THEN_ATG,
        Change("CCC[T>C]AAggg"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(44, 20), per_strand(45, 21)),
            "alt_cds_start": per_strand(12, 18),
            "alt_cds_stop": per_strand(47, 53),
            "alt_cds_seq": "CTGAAAATGCCCCAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 5), (2, 10)],
            **columns(ALT_CODON_COLUMNS, 0, 1, "CAA", False, None, None, 0, [], []),
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": "GGCTGAAAATGCCCCAAGGGTAGCC",
            "alt_transcript_length": 25,
            "alt_cds_start_in_transcript": 2,
            **columns(SCAN_COLUMNS, 2, 1, "GCC", False, "TAG", 20, 1, [(20, "TAG")], [2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": -6,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        ruler=Ruler((0, 2, 7, 8, 14, 20, 25)),
    ),
    # SC-16
    Case(
        "start_loss_with_a_ptc_3_nt_after_the_atg_of_the_scan_escapes_by_the_start_proximal_rule",
        """
        tx      0   3       9   12                           63                   123 126      134
        ref 5' [ggg ATG GCC ATG CAG AAA AAA ..36.. AAA AAA]|[AAA AAA ..51.. AAA]|[GAC TAA ggggg] 3'
        alt 5' [ggg ACG GCC ATG TAG AAA AAA ..36.. AAA AAA]|[AAA AAA ..51.. AAA]|[GAC TAA ggggg] 3'
                     ^^^^^^^^^^^^ TGGCCATGC>CGGCCATGT
                            aaa ATG of the scan
                                *** PTC
                            <-> ptc_to_start_codon = 3
                                <-- ptc_to_intron = 51 -->
        exon 1: 63 nt, exon 2: 60 nt, exon 3: 11 nt
        One MNV changes ATG>ACG and CAG>TAG. The scan finds the ATG at tx 9, and its first stop codon is the TAG at
        tx 12: a PTC, 114 nt upstream of the stop codon at tx 126. From the CDS start, the PTC would be 9 nt away.
        """,
        MET_RESCUE,
        Change("gggA[TGGCCATGC>CGGCCATGT]AGAAA"),
        {
            **variant(
                per_strand("TGGCCATGC", "GCATGGCCA"),
                per_strand("CGGCCATGT", "ACATGGCCG"),
                per_strand(14, 171),
                per_strand(23, 180),
            ),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(179, 181),
            "alt_cds_seq": "ACGGCCATGTAG" + "AAA" * 36 + "GACTAA",
            "alt_cds_len": 126,
            "alt_cds_info": [(1, 60), (2, 60), (3, 6)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAG", 9, 2, [(9, "TAG"), (123, "TAA")], [1, 3]),
            "alt_is_premature": True,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGACGGCCATGTAG" + "AAA" * 36 + "GACTAAGGGGG",
            "alt_transcript_length": 134,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 9, 1, "GGG", False, "TAG", 12, 2, [(12, "TAG"), (126, "TAA")], [1, 3]),
            "unknown_reason": None,
            **columns(PTC_FEATURE_COLUMNS, 0, 2, 3, True, 63, 51),
            "stop_codon_distance": 114,
            **columns(RULE_COLUMNS, False, False, False, True, False, True),
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 9, 12, "a", "ATG of the scan"),
            Mark("alt", 12, 15, "*", "PTC"),
            Span("alt", 9, 12, "ptc_to_start_codon = 3"),
            Span("alt", 12, 63, "ptc_to_intron = 51"),
        ),
        ruler=Ruler((0, 3, 9, 12, 63, 123, 126, 134)),
    ),
    # SC-08, PF-08, PF-11, PF-14, PF-21
    Case(
        "start_loss_classifies_the_first_stop_codon_of_an_out_of_frame_rescued_orf_as_a_ptc",
        """
        tx      0   3    7       13                          63                   123                    168      176
        ref 5' [ggg ATG GAT GCC CTA AAA AAA ..36.. AAA AAA]|[AAA AAA ..51.. AAA]|[AAA AAA ..33.. AAA GAC TAA ggggg] 3'
        alt 5' [ggg ACG GAT GCC CTA AAA AAA ..36.. AAA AAA]|[AAA AAA ..51.. AAA]|[AAA AAA ..33.. AAA GAC TAA ggggg] 3'
                     ^ T>C
                         aaaa ATG of the scan
                                 **** PTC
                         <------> ptc_to_start_codon = 6
                                 <--------------------- stop_codon_distance = 155 --------------------->
                                 <-----------------------> ptc_to_intron = 50
        exon 1: 63 nt, exon 2: 60 nt, exon 3: 53 nt
        ATG>ACG is a start loss. The scan finds the ATG at tx 7, out of frame with the stop codon. Its ORF ATG CCC TAA
        ends at tx 13, upstream of the stop codon at tx 168: a PTC in exon 1.
        """,
        OUT_OF_FRAME_RESCUE,
        Change("gggA[T>C]GGATG"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(14, 221), per_strand(15, 222)),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(221, 223),
            "alt_cds_seq": "ACGGATGCCCTA" + "AAA" * 50 + "GACTAA",
            "alt_cds_len": 168,
            "alt_cds_info": [(1, 60), (2, 60), (3, 48)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 165, 1, [(165, "TAA")], [3]),
            "alt_is_premature": True,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGACGGATGCCCTA" + "AAA" * 50 + "GACTAAGGGGG",
            "alt_transcript_length": 176,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 7, 1, "GGG", False, "TAA", 13, 1, [(13, "TAA")], [1]),
            "unknown_reason": None,
            **columns(PTC_FEATURE_COLUMNS, 0, 2, 6, True, 63, 50),
            "stop_codon_distance": 155,
            **columns(RULE_COLUMNS, False, False, False, True, False, True),
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 7, 10, "a", "ATG of the scan"),
            Mark("alt", 13, 16, "*", "PTC"),
            Span("alt", 7, 13, "ptc_to_start_codon = 6"),
            Span("alt", 13, 168, "stop_codon_distance = 155"),
            Span("alt", 13, 63, "ptc_to_intron = 50"),
        ),
        ruler=Ruler((0, 3, 7, 13, 63, 123, 168, 176)),
    ),
    # SC-05, SC-09
    Case(
        "deletion_in_the_atg_is_a_start_loss_and_the_next_atg_in_frame_reads_to_the_annotated_stop_codon",
        """
        ref    5' [ggg ATG GT]|[A AGC ATG GCC AAA GAC TAA ggggg] 3'
        alt    5' [ggg A-G GT]|[A AGC ATG GCC AAA GAC TAA ggggg] 3'
                        ^ T>-
        alt tx     0   3    6   7     11              23       31
                            sssssss TAA in the frame of the CDS start
                                      aaa ATG of the scan
                                                      sss the annotated stop codon
        No placement of the deletion keeps an ATG at the CDS start. The scan finds the ATG at tx 11, which is the Met
        at CDS 9 of the ref CDS. Its ORF ends at the annotated stop codon at tx 23: stop_codon_distance = 0.
        """,
        FRAMESHIFT_RESCUE,
        Change("gggA[T>]GGT"),
        {
            **variant(per_strand("AT", "CA"), per_strand("A", "C"), per_strand(13, 56), per_strand(15, 58)),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(57, 59),
            "alt_cds_seq": "AGGTAAGCATGGCCAAAGACTAA",
            "alt_cds_len": 23,
            "alt_cds_info": [(1, 4), (2, 19)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 3, 1, [(3, "TAA")], [1]),
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGAGGTAAGCATGGCCAAAGACTAAGGGGG",
            "alt_transcript_length": 31,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 11, 2, "GGG", False, "TAA", 23, 1, [(23, "TAA")], [2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": [(1, 7), (2, 24)],
        },
        equivalent=(Change("gggA[TG>G]GT"),),
        marks=(
            Ruler((0, 3, 6, 7, 11, 23, 31), "tx", "alt"),
            Mark("alt", 6, 9, "s", "TAA in the frame of the CDS start"),
            Mark("alt", 11, 14, "a", "ATG of the scan"),
            Mark("alt", 23, 26, "s", "the annotated stop codon"),
        ),
    ),
    # SC-12
    Case(
        "start_loss_with_the_next_atg_in_the_3utr_is_neither_a_ptc_nor_a_stop_loss",
        """
        tx      0   3        8         15    20    26   31
        ref 5' [ggg ATG CC]|[C AAA GAC TAA ggatgccctgagg] 3'
        alt 5' [ggg ACG CC]|[C AAA GAC TAA ggatgccctgagg] 3'
                     ^ T>C
                                       sss the annotated stop codon
                                             aaa ATG of the scan
                                                   sss the stop codon of the scan
        The scan finds the ATG at tx 20, downstream of the first base of the stop codon at tx 15. Its ORF does not
        overlap the CDS: neither flag, stop_codon_distance = null. The scan columns are set.
        """,
        ATG_IN_THE_3UTR,
        Change("gggA[T>C]GCC"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(14, 56), per_strand(15, 57)),
            "alt_cds_start": per_strand(13, 23),
            "alt_cds_stop": per_strand(48, 58),
            "alt_cds_seq": "ACGCCCAAAGACTAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 5), (2, 10)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 12, 1, [(12, "TAA")], [2]),
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGACGCCCAAAGACTAAGGATGCCCTGAGG",
            "alt_transcript_length": 31,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 20, 2, "AGG", False, "TGA", 26, 1, [(26, "TGA")], [2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 15, 18, "s", "the annotated stop codon"),
            Mark("alt", 20, 23, "a", "ATG of the scan"),
            Mark("alt", 26, 29, "s", "the stop codon of the scan"),
        ),
        ruler=Ruler((0, 3, 8, 15, 20, 26, 31)),
    ),
    # SC-10
    Case(
        "start_loss_whose_rescued_orf_ends_downstream_of_the_annotated_stop_codon_is_a_stop_loss",
        """
        tx      0   3    7       11            21   25   30
        ref 5' [ggg ATG GAT GC]|[C CCC AAA GAC TAA gtgacc] 3'
        alt 5' [ggg ACG GAT GC]|[C CCC AAA GAC TAA gtgacc] 3'
                     ^ T>C
                         aaaa ATG of the scan
                                               sss the annotated stop codon
                                                    sss the stop codon of the scan
                                               <---> stop_codon_distance = 21 - 25 = -4
        The scan finds the ATG at tx 7, out of frame with the stop codon at tx 21. Its frame reads past it to the TGA
        at tx 25 in the 3' UTR: a stop loss.
        """,
        RESCUE_PAST_THE_STOP,
        Change("gggA[T>C]GGATG"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(14, 55), per_strand(15, 56)),
            "alt_cds_start": per_strand(13, 16),
            "alt_cds_stop": per_strand(54, 57),
            "alt_cds_seq": "ACGGATGCCCCCAAAGACTAA",
            "alt_cds_len": 21,
            "alt_cds_info": [(1, 8), (2, 13)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 18, 1, [(18, "TAA")], [2]),
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": True,
            "alt_transcript_seq": "GGGACGGATGCCCCCAAAGACTAAGTGACC",
            "alt_transcript_length": 30,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 7, 1, "ACC", False, "TGA", 25, 1, [(25, "TGA")], [2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": -4,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 7, 10, "a", "ATG of the scan"),
            Mark("alt", 21, 24, "s", "the annotated stop codon"),
            Mark("alt", 25, 28, "s", "the stop codon of the scan"),
            Span("alt", 21, 25, "stop_codon_distance = 21 - 25 = -4"),
        ),
        ruler=Ruler((0, 3, 7, 11, 21, 25, 30)),
    ),
    # SC-13, the side of an ATG upstream of the first base of the stop codon
    Case(
        "start_loss_with_the_next_atg_1_nt_before_the_stop_codon_is_classified",
        """
        tx      0   3        8         15    20   25
        ref 5' [ggg ATG CC]|[C AAA CCA TGA cctaggg] 3'
        alt 5' [ggg ACG CC]|[C AAA CCA TGA cctaggg] 3'
                     ^ T>C
                                     aaaa ATG of the scan
                                       sss the annotated stop codon
                                             sss the stop codon of the scan
        The scan finds the ATG at tx 14, 1 nt before the first base of the stop codon TGA at tx 15. So its ORF is
        classified: its frame reads to the TAG at tx 20, a stop loss. stop_codon_distance = 15 - 20 = -5.
        """,
        ATG_BEFORE_THE_STOP_CODON,
        Change("gggA[T>C]GCC"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(14, 50), per_strand(15, 51)),
            "alt_cds_start": per_strand(13, 17),
            "alt_cds_stop": per_strand(48, 52),
            "alt_cds_seq": "ACGCCCAAACCATGA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 5), (2, 10)],
            **columns(ALT_CODON_COLUMNS, None, None, "TGA", True, "TGA", 12, 1, [(12, "TGA")], [2]),
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": True,
            "alt_transcript_seq": "GGGACGCCCAAACCATGACCTAGGG",
            "alt_transcript_length": 25,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 14, 2, "GGG", False, "TAG", 20, 1, [(20, "TAG")], [2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": -5,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 14, 17, "a", "ATG of the scan"),
            Mark("alt", 15, 18, "s", "the annotated stop codon"),
            Mark("alt", 20, 23, "s", "the stop codon of the scan"),
        ),
        ruler=Ruler((0, 3, 8, 15, 20, 25)),
    ),
    # SC-13, the side of an ATG downstream of the first base of the stop codon
    Case(
        "start_loss_with_the_next_atg_at_the_last_base_of_the_stop_codon_has_no_orf",
        """
        tx      0   3        8           17     23   28
        ref 5' [ggg ATG CC]|[C AAA GAC TAA tgccctgagg] 3'
        alt 5' [ggg ACG CC]|[C AAA GAC TAA tgccctgagg] 3'
                     ^ T>C
                                       sss the annotated stop codon
                                         aaaa ATG of the scan
                                                sss the stop codon of the scan
        The scan finds the ATG at tx 17, 2 nt after the first base of the stop codon at tx 15. An ATG 1 nt after it
        would need TAT, which is no stop codon. So this is the nearest ATG on this side: its ORF does not overlap the
        CDS, neither flag, stop_codon_distance = null.
        """,
        ATG_AT_THE_LAST_STOP_CODON_BASE,
        Change("gggA[T>C]GCC"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(14, 53), per_strand(15, 54)),
            "alt_cds_start": per_strand(13, 20),
            "alt_cds_stop": per_strand(48, 55),
            "alt_cds_seq": "ACGCCCAAAGACTAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 5), (2, 10)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 12, 1, [(12, "TAA")], [2]),
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGACGCCCAAAGACTAATGCCCTGAGG",
            "alt_transcript_length": 28,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 17, 2, "AGG", False, "TGA", 23, 1, [(23, "TGA")], [2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 15, 18, "s", "the annotated stop codon"),
            Mark("alt", 17, 20, "a", "ATG of the scan"),
            Mark("alt", 23, 26, "s", "the stop codon of the scan"),
        ),
        ruler=Ruler((0, 3, 8, 17, 23, 28)),
    ),
]

# SC-06: ATG and 12 C, so that the frame +1 of the CDS reads past the stop codon to a TGA in the 3' UTR
ATG_C12 = Layout(
    Transcript(("gggATGCCCCC", "CCCCCCCTAAggtgacc")),
    {
        **ATG_TAA_BASE,
        "ref_cds_start": per_strand(13, 17),
        "ref_cds_stop": per_strand(51, 55),
        "ref_cds_seq": "ATG" + "C" * 12 + "TAA",
        "ref_cds_len": 18,
        "ref_cds_info": [(1, 8), (2, 10)],
        "ref_first_stop_pos": 15,
        "ref_all_stop_codons": [(15, "TAA")],
        "ref_stop_codon_exons": [2],
        "transcript_end": 58,
        "transcript_seq": "GGGATG" + "C" * 12 + "TAAGGTGACC",
        "transcript_length": 28,
        "cds_end_in_transcript": 21,
        "transcript_exon_info": [(1, 11), (2, 17)],
        "utr3_length": 7,
        "total_exon_count": 2,
    },
)

# SC-07, SC-21: a 13 nt 5' UTR with an ATG at tx 2, ATG GCT ATG GCT GCT GCT TGA, and a TAG in the frame of the CDS in
# the 3' UTR. The 5' UTR is one exon or two.
SCAN_START_TRANSCRIPT = "CCATGCCGCCGCC" + "ATGGCTATGGCTGCTGCTTGA" + "GCTGCTTAGCCTAACCC"
SCAN_START_BASE = {
    **ATG_TAA_BASE,
    "ref_cds_start": per_strand(23, 47),
    "ref_cds_seq": "ATGGCTATGGCTGCTGCTTGA",
    "ref_cds_len": 21,
    "ref_last_codon": "TGA",
    "ref_first_stop_codon": "TGA",
    "ref_first_stop_pos": 18,
    "ref_all_stop_codons": [(18, "TGA")],
    "transcript_seq": SCAN_START_TRANSCRIPT,
    "transcript_length": 51,
    "cds_start_in_transcript": 13,
    "cds_end_in_transcript": 34,
    "utr3_length": 17,
    "utr5_length": 13,
}
SCAN_START = Layout(
    Transcript(("ccatgccgccgccATGGCTATGGCTGCTGCTTGA", "gctgcttagcctaaccc")),
    {
        **SCAN_START_BASE,
        "ref_cds_stop": per_strand(44, 68),
        "ref_cds_info": [(1, 21)],
        "ref_stop_codon_exons": [1],
        "transcript_end": 81,
        "transcript_exon_info": [(1, 34), (2, 17)],
        "total_exon_count": 2,
    },
)
SCAN_START_AFTER_A_5UTR_INTRON = Layout(
    Transcript(("ccatg", "ccgccgccATGGCTATGGCTGCTGCTTGA", "gctgcttagcctaaccc")),
    {
        **SCAN_START_BASE,
        "ref_cds_start": per_strand(43, 47),
        "ref_cds_stop": per_strand(64, 68),
        "ref_cds_info": [(2, 21)],
        "ref_start_codon_exon": 2,
        "ref_stop_codon_exons": [2],
        "transcript_end": 101,
        "transcript_exon_info": [(1, 5), (2, 29), (3, 17)],
        "total_exon_count": 3,
    },
)

# SC-14: one exon; the 3' UTR holds the next ATG and a TGA in its frame
STOP_THEN_ATG = Layout(
    Transcript(("ccATGAAACCCTAAgatgccctgacc",)),
    {
        **ATG_TAA_BASE,
        "ref_cds_start": per_strand(12, 22),
        "ref_cds_stop": per_strand(24, 34),
        "ref_cds_seq": "ATGAAACCCTAA",
        "ref_cds_len": 12,
        "ref_cds_info": [(1, 12)],
        "ref_first_stop_pos": 9,
        "ref_all_stop_codons": [(9, "TAA")],
        "ref_stop_codon_exons": [1],
        "transcript_end": 36,
        "transcript_seq": "CCATGAAACCCTAAGATGCCCTGACC",
        "transcript_length": 26,
        "cds_start_in_transcript": 2,
        "cds_end_in_transcript": 14,
        "transcript_exon_info": [(1, 26)],
        "utr3_length": 12,
        "utr5_length": 2,
        "total_exon_count": 1,
    },
)

# SC-15: a cds_end_NF transcript without stop_codon rows, whose CDS ends in GAC. Its ATG at CDS 4 reads CTA AAA as TAA,
# or its frame reads CCA AAA and finds a TAA only in the 3' region after the CDS.
NO_STOP_CODON_BASE = {
    **ATG_TAA_BASE,
    "ref_cds_len": 24,
    "has_stop_codon": False,
    "ref_cds_info": [(1, 15), (2, 9)],
    "ref_last_codon": "GAC",
    "ref_valid_stop": False,
    "ref_first_stop_codon": None,
    "ref_first_stop_pos": None,
    "ref_num_stop_codons": 0,
    "ref_all_stop_codons": [],
    "ref_stop_codon_exons": [],
    "cds_end_in_transcript": 27,
    "utr3_length": None,
    "total_exon_count": 2,
    "likely_misannotated": True,
}
NO_STOP_CODON_PTC = Layout(
    Transcript(("gggATGGATGCCCTAAAA", "AAAAAAGACcccc"), stop_codon=False, tags=("cds_end_NF",)),
    {
        **NO_STOP_CODON_BASE,
        "ref_cds_start": per_strand(13, 14),
        "ref_cds_stop": per_strand(57, 58),
        "ref_cds_seq": "ATGGATGCCCTA" + "AAA" * 3 + "GAC",
        "transcript_end": 61,
        "transcript_seq": "GGGATGGATGCCCTA" + "AAA" * 3 + "GACCCCC",
        "transcript_length": 31,
        "transcript_exon_info": [(1, 18), (2, 13)],
    },
)
NO_STOP_CODON_STOP_PAST_THE_CDS = Layout(
    Transcript(("gggATGGATGCCCCAAAA", "AAAAAAGACctaacc"), stop_codon=False, tags=("cds_end_NF",)),
    {
        **NO_STOP_CODON_BASE,
        "ref_cds_start": per_strand(13, 16),
        "ref_cds_stop": per_strand(57, 60),
        "ref_cds_seq": "ATGGATGCCCCA" + "AAA" * 3 + "GAC",
        "transcript_end": 63,
        "transcript_seq": "GGGATGGATGCCCCA" + "AAA" * 3 + "GACCTAACC",
        "transcript_length": 33,
        "transcript_exon_info": [(1, 18), (2, 15)],
    },
)

# SC-17: as OUT_OF_FRAME_RESCUE, with the CTA 150 nt after the ATG at CDS 4, in an exon 1 of 173 nt
START_PROXIMAL_150_CDS = "ATGGATGCC" + "AAA" * 48 + "CTA" + "AAA" * 30 + "GACTAA"
START_PROXIMAL_150 = Layout(
    Transcript(
        ("ggg" + START_PROXIMAL_150_CDS[:170], START_PROXIMAL_150_CDS[170:240], START_PROXIMAL_150_CDS[240:] + "ggggg")
    ),
    {
        **ATG_TAA_BASE,
        "ref_cds_stop": per_strand(305, 307),
        "ref_cds_seq": START_PROXIMAL_150_CDS,
        "ref_cds_len": 252,
        "ref_cds_info": [(1, 170), (2, 70), (3, 12)],
        "ref_first_stop_pos": 249,
        "ref_all_stop_codons": [(249, "TAA")],
        "ref_stop_codon_exons": [3],
        "transcript_end": 310,
        "transcript_seq": "GGG" + START_PROXIMAL_150_CDS + "GGGGG",
        "transcript_length": 260,
        "cds_end_in_transcript": 255,
        "transcript_exon_info": [(1, 173), (2, 70), (3, 17)],
        "total_exon_count": 3,
    },
)


# SC-19: a 5' UTR gggtca, and the ATG at CDS 4 that reads CTA AAA as TAA
UTR5_AND_START_CODON = "ATGGATGCCCTA" + "AAA" * 3 + "GACTAA"
UTR5_AND_START = Layout(
    Transcript(("gggtca" + UTR5_AND_START_CODON[:16], UTR5_AND_START_CODON[16:] + "ggggg")),
    {
        **ATG_TAA_BASE,
        "ref_cds_start": per_strand(16, 15),
        "ref_cds_stop": per_strand(63, 62),
        "ref_cds_seq": UTR5_AND_START_CODON,
        "ref_cds_len": 27,
        "ref_cds_info": [(1, 16), (2, 11)],
        "ref_first_stop_pos": 24,
        "ref_all_stop_codons": [(24, "TAA")],
        "ref_stop_codon_exons": [2],
        "transcript_end": 68,
        "transcript_seq": "GGGTCA" + UTR5_AND_START_CODON + "GGGGG",
        "transcript_length": 38,
        "cds_start_in_transcript": 6,
        "cds_end_in_transcript": 33,
        "transcript_exon_info": [(1, 22), (2, 16)],
        "utr5_length": 6,
        "total_exon_count": 2,
    },
)

CASES += [
    # SC-06
    Case(
        "insertion_inside_the_start_codon_that_keeps_an_atg_at_the_cds_start_is_no_start_loss",
        """
        ref    5' [ggg A----TG CCC CC]|[C CCC CCC TAA ggtgacc] 3'
        alt    5' [ggg ATGATTG CCC CC]|[C CCC CCC TAA ggtgacc] 3'
                        ^^^^ ->TGAT
        alt tx     0   3                          22    27   32
                                                  sss the annotated stop codon, shifted
                                                        sss the stop codon of the scan
        A>ATGAT keeps an ATG at the CDS start, so it is no start loss. The 4 nt frameshift reads past the shifted stop
        codon at tx 22 to the TGA at tx 27 in the 3' UTR: a stop loss, stop_codon_distance = 22 - 27 = -5.
        The insertion also fits after AT (GATT) and after ATG (ATTG).
        """,
        ATG_C12,
        Change("gggA[>TGAT]TGCCC"),
        {
            **variant(per_strand("A", "A"), per_strand("ATGAT", "AATCA"), per_strand(13, 53), per_strand(14, 54)),
            "alt_cds_start": per_strand(13, 17),
            "alt_cds_stop": per_strand(51, 55),
            "alt_cds_seq": "ATGATTG" + "C" * 12 + "TAA",
            "alt_cds_len": 22,
            "alt_cds_info": [(1, 12), (2, 10)],
            **columns(ALT_CODON_COLUMNS, 0, 1, "TAA", True, None, None, 0, [], []),
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": "GGGATGATTG" + "C" * 12 + "TAAGGTGACC",
            "alt_transcript_length": 32,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 3, 1, "ACC", False, "TGA", 27, 1, [(27, "TGA")], [2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": -5,
            **NO_RULE,
            "alt_transcript_exon_info": [(1, 15), (2, 17)],
        },
        equivalent=(Change("gggAT[>GATT]GCCC"), Change("gggATG[>ATTG]CCCC")),
        marks=(
            Ruler((0, 3, 22, 27, 32), "tx", "alt"),
            Mark("alt", 22, 25, "s", "the annotated stop codon, shifted"),
            Mark("alt", 27, 30, "s", "the stop codon of the scan"),
        ),
    ),
    # SC-07, SC-21
    Case(
        "start_loss_scan_skips_the_atg_in_the_5utr_and_takes_the_first_atg_from_the_cds_start",
        """
        tx      0 2           13      19              31    34    40         51
        ref 5' [ccatgccgccgcc ATG GCT ATG GCT GCT GCT TGA]|[gctgcttagcctaaccc] 3'
        alt 5' [ccatgccgccgcc ATA GCT ATG GCT GCT GCT TGA]|[gctgcttagcctaaccc] 3'
                                ^ G>A
                                      aaa ATG of the scan
                                                      sss stop codon of the scan
                                                                  sss stop codon of the scan
        exon 1: 34 nt, exon 2: 17 nt
        ATG>ATA is a start loss. The scan starts at the CDS start, tx 13, so it skips the ATG at tx 2 in the 5' UTR.
        It finds the ATG at tx 19, in frame with the stop codon, and reads a stop codon in exon 1 and one in exon 2.
        """,
        SCAN_START,
        Change("gccAT[G>A]GCTATG"),
        {
            **variant(per_strand("G", "C"), per_strand("A", "T"), per_strand(25, 65), per_strand(26, 66)),
            "alt_cds_start": per_strand(23, 47),
            "alt_cds_stop": per_strand(44, 68),
            "alt_cds_seq": "ATAGCTATGGCTGCTGCTTGA",
            "alt_cds_len": 21,
            "alt_cds_info": [(1, 21)],
            **columns(ALT_CODON_COLUMNS, None, None, "TGA", True, "TGA", 18, 1, [(18, "TGA")], [1]),
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "CCATGCCGCCGCC" + "ATAGCTATGGCTGCTGCTTGA" + "GCTGCTTAGCCTAACCC",
            "alt_transcript_length": 51,
            "alt_cds_start_in_transcript": 13,
            **columns(SCAN_COLUMNS, 19, 1, "CCC", False, "TGA", 31, 2, [(31, "TGA"), (40, "TAG")], [1, 2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 19, 22, "a", "ATG of the scan"),
            Mark("alt", 31, 34, "s", "stop codon of the scan"),
            Mark("alt", 40, 43, "s", "stop codon of the scan"),
        ),
        ruler=Ruler((0, 2, 13, 19, 31, 34, 40, 51)),
    ),
    # SC-07, SC-21
    Case(
        "start_loss_scan_after_an_intron_in_the_5utr_takes_the_first_atg_from_the_cds_start",
        """
        tx      0 2     5        13      19              31    34    40         51
        ref 5' [ccatg]|[ccgccgcc ATG GCT ATG GCT GCT GCT TGA]|[gctgcttagcctaaccc] 3'
        alt 5' [ccatg]|[ccgccgcc ATA GCT ATG GCT GCT GCT TGA]|[gctgcttagcctaaccc] 3'
                                   ^ G>A
                                         aaa ATG of the scan
                                                         sss stop codon of the scan
                                                                     sss stop codon of the scan
        exon 1: 5 nt, exon 2: 29 nt, exon 3: 17 nt
        As without the intron, in transcript positions. The ATG of the scan lies in exon 2, its stop codons in exons 2
        and 3.
        """,
        SCAN_START_AFTER_A_5UTR_INTRON,
        Change("gccAT[G>A]GCTATG"),
        {
            **variant(per_strand("G", "C"), per_strand("A", "T"), per_strand(45, 65), per_strand(46, 66)),
            "alt_cds_start": per_strand(43, 47),
            "alt_cds_stop": per_strand(64, 68),
            "alt_cds_seq": "ATAGCTATGGCTGCTGCTTGA",
            "alt_cds_len": 21,
            "alt_cds_info": [(2, 21)],
            **columns(ALT_CODON_COLUMNS, None, None, "TGA", True, "TGA", 18, 1, [(18, "TGA")], [2]),
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "CCATGCCGCCGCC" + "ATAGCTATGGCTGCTGCTTGA" + "GCTGCTTAGCCTAACCC",
            "alt_transcript_length": 51,
            "alt_cds_start_in_transcript": 13,
            **columns(SCAN_COLUMNS, 19, 2, "CCC", False, "TGA", 31, 2, [(31, "TGA"), (40, "TAG")], [2, 3]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 19, 22, "a", "ATG of the scan"),
            Mark("alt", 31, 34, "s", "stop codon of the scan"),
            Mark("alt", 40, 43, "s", "stop codon of the scan"),
        ),
        ruler=Ruler((0, 2, 5, 13, 19, 31, 34, 40, 51)),
    ),
    # SC-14
    Case(
        "start_loss_with_a_deleted_stop_codon_reads_from_the_next_atg_in_the_former_3utr",
        """
        ref    5' [cc ATG AAA CCC TAA gatgccctgacc] 3'
        alt    5' [cc A-- --- --- --- gatgccctgacc] 3'
                       ^^^^^^^^^^^^^^ 11 nt>-
        alt tx     0  2                4     10   15
                                       aaa ATG of the scan
                                             sss the stop codon of the scan
        The deletion of TGAAACCCTAA removes the start codon and the stop codon; the alt CDS is the A at tx 2. The scan
        finds the ATG at tx 4, in the former 3' UTR, downstream of the stop codon. Its ORF does not overlap the CDS:
        neither flag. The deletion also fits 1 nt to the left, as ATGAAACCCTA.
        """,
        STOP_THEN_ATG,
        Change("ccA[TGAAACCCTAA>]gatg"),
        {
            **variant(
                per_strand("ATGAAACCCTAA", "CTTAGGGTTTCA"), per_strand("A", "C"), per_strand(12, 21), per_strand(24, 33)
            ),
            "alt_cds_start": per_strand(12, 22),
            "alt_cds_stop": per_strand(24, 34),
            "alt_cds_seq": "A",
            "alt_cds_len": 1,
            "alt_cds_info": [(1, 1)],
            **columns(ALT_CODON_COLUMNS, *[None] * 9),
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "CCAGATGCCCTGACC",
            "alt_transcript_length": 15,
            "alt_cds_start_in_transcript": 2,
            **columns(SCAN_COLUMNS, 4, 1, "ACC", False, "TGA", 10, 1, [(10, "TGA")], [1]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": [(1, 15)],
        },
        equivalent=(Change("cc[ATGAAACCCTA>]Agatg"),),
        marks=(
            Ruler((0, 2, 4, 10, 15), "tx", "alt"),
            Mark("alt", 4, 7, "a", "ATG of the scan"),
            Mark("alt", 10, 13, "s", "the stop codon of the scan"),
        ),
    ),
    # SC-15, the side of a rescued stop codon inside the alt CDS
    Case(
        "start_loss_without_an_annotated_stop_codon_and_a_rescued_stop_inside_the_alt_cds_is_a_ptc",
        """
        tx      0   3    7       13       18          27  31
        ref 5' [ggg ATG GAT GCC CTA AAA]|[AAA AAA GAC cccc] 3'
        alt 5' [ggg ACG GAT GCC CTA AAA]|[AAA AAA GAC cccc] 3'
                     ^ T>C
                         aaaa ATG of the scan
                                 **** PTC
                         <------> ptc_to_start_codon = 6
                                 <----> ptc_to_intron = 5
        exon 1: 18 nt, exon 2: 13 nt; cds_end_NF, no stop_codon rows
        Without an annotated stop codon, the TAA at tx 13 of the rescued ORF is a PTC because it lies inside the alt
        CDS, which ends at tx 27. stop_codon_distance = null.
        """,
        NO_STOP_CODON_PTC,
        Change("gggA[T>C]GGATG"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(14, 56), per_strand(15, 57)),
            "alt_cds_start": per_strand(13, 14),
            "alt_cds_stop": per_strand(57, 58),
            "alt_cds_seq": "ACGGATGCCCTA" + "AAA" * 3 + "GAC",
            "alt_cds_len": 24,
            "alt_cds_info": [(1, 15), (2, 9)],
            **columns(ALT_CODON_COLUMNS, None, None, "GAC", False, None, None, 0, [], []),
            "alt_is_premature": True,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGACGGATGCCCTA" + "AAA" * 3 + "GACCCCC",
            "alt_transcript_length": 31,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 7, 1, "CCC", False, "TAA", 13, 1, [(13, "TAA")], [1]),
            "unknown_reason": None,
            **columns(PTC_FEATURE_COLUMNS, 0, 1, 6, True, 18, 5),
            "stop_codon_distance": None,
            **columns(RULE_COLUMNS, False, True, False, True, False, True),
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 7, 10, "a", "ATG of the scan"),
            Mark("alt", 13, 16, "*", "PTC"),
            Span("alt", 7, 13, "ptc_to_start_codon = 6"),
            Span("alt", 13, 18, "ptc_to_intron = 5"),
        ),
        ruler=Ruler((0, 3, 7, 13, 18, 27, 31)),
    ),
    # SC-15, the side of a rescued stop codon past the end of the alt CDS
    Case(
        "start_loss_without_an_annotated_stop_codon_and_a_rescued_stop_past_the_alt_cds_is_neither",
        """
        tx      0   3    7                18           28   33
        ref 5' [ggg ATG GAT GCC CCA AAA]|[AAA AAA GAC ctaacc] 3'
        alt 5' [ggg ACG GAT GCC CCA AAA]|[AAA AAA GAC ctaacc] 3'
                     ^ T>C
                         aaaa ATG of the scan
                                                       sss the stop codon of the scan
        exon 1: 18 nt, exon 2: 15 nt; cds_end_NF, no stop_codon rows
        The frame of the ATG at tx 7 reads its first stop codon at tx 28, past the end of the alt CDS at tx 27:
        neither flag, stop_codon_distance = null.
        """,
        NO_STOP_CODON_STOP_PAST_THE_CDS,
        Change("gggA[T>C]GGATG"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(14, 58), per_strand(15, 59)),
            "alt_cds_start": per_strand(13, 16),
            "alt_cds_stop": per_strand(57, 60),
            "alt_cds_seq": "ACGGATGCCCCA" + "AAA" * 3 + "GAC",
            "alt_cds_len": 24,
            "alt_cds_info": [(1, 15), (2, 9)],
            **columns(ALT_CODON_COLUMNS, None, None, "GAC", False, None, None, 0, [], []),
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGACGGATGCCCCA" + "AAA" * 3 + "GACCTAACC",
            "alt_transcript_length": 33,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 7, 1, "ACC", False, "TAA", 28, 1, [(28, "TAA")], [2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(Mark("alt", 7, 10, "a", "ATG of the scan"), Mark("alt", 28, 31, "s", "the stop codon of the scan")),
        ruler=Ruler((0, 3, 7, 18, 28, 33)),
    ),
    # SC-17
    Case(
        "start_loss_with_a_ptc_150_nt_after_the_atg_of_the_scan_is_not_start_proximal",
        """
        tx      0   3    7                       157                     173                243         252      260
        ref 5' [ggg ATG GAT GCC AAA ..138.. AAA CTA AAA AAA AAA AAA AA]|[A AAA ..63.. AAA]|[AAA AAA GAC TAA ggggg] 3'
        alt 5' [ggg ACG GAT GCC AAA ..138.. AAA CTA AAA AAA AAA AAA AA]|[A AAA ..63.. AAA]|[AAA AAA GAC TAA ggggg] 3'
                     ^ T>C
                         aaaa ATG of the scan
                                                 **** PTC
                         <----------------------> ptc_to_start_codon = 150, not < 150
                                                 <-------------------> ptc_to_intron = 16
        exon 1: 173 nt, exon 2: 70 nt, exon 3: 17 nt
        The frame of the ATG at tx 7 reads CTA AAA as TAA at tx 157: a PTC 150 nt after the ATG, so no NMD escape rule
        applies.
        """,
        START_PROXIMAL_150,
        Change("gggA[T>C]GGATG"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(14, 305), per_strand(15, 306)),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(305, 307),
            "alt_cds_seq": "ACGGATGCC" + "AAA" * 48 + "CTA" + "AAA" * 30 + "GACTAA",
            "alt_cds_len": 252,
            "alt_cds_info": [(1, 170), (2, 70), (3, 12)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 249, 1, [(249, "TAA")], [3]),
            "alt_is_premature": True,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGACGGATGCC" + "AAA" * 48 + "CTA" + "AAA" * 30 + "GACTAAGGGGG",
            "alt_transcript_length": 260,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 7, 1, "GGG", False, "TAA", 157, 1, [(157, "TAA")], [1]),
            "unknown_reason": None,
            **columns(PTC_FEATURE_COLUMNS, 0, 2, 150, False, 173, 16),
            "stop_codon_distance": 95,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 7, 10, "a", "ATG of the scan"),
            Mark("alt", 157, 160, "*", "PTC"),
            Span("alt", 7, 157, "ptc_to_start_codon = 150, not < 150"),
            Span("alt", 157, 173, "ptc_to_intron = 16"),
        ),
        ruler=Ruler((0, 3, 7, 157, 173, 243, 252, 260)),
    ),
    # SC-19
    Case(
        "deletion_of_5utr_bases_and_the_start_codon_moves_the_alt_cds_start_and_the_scan_reads_from_there",
        """
        ref    5' [gggtca ATG GAT GCC CTA AAA A]|[AA AAA GAC TAA ggggg] 3'
        alt    5' [gggt-- --- GAT GCC CTA AAA A]|[AA AAA GAC TAA ggggg] 3'
                       ^^^^^^ caATG>-
        alt tx     0          4        11         17         25       33
                               aaaa ATG of the scan
                                       **** PTC
                                                             sss the annotated stop codon
                               <------> ptc_to_start_codon = 6
                                       <------> ptc_to_intron = 6
        exon 1: 22 nt, exon 2: 16 nt in the ref transcript; exon 1: 17 nt in the alt transcript, the PTC exon
        The 5' UTR loses 2 nt, so alt_cds_start_in_transcript = 4. From there, the scan finds the ATG at tx 5. Its ORF
        ends at the TAA at tx 11: a PTC, 14 nt upstream of the stop codon at tx 25.
        """,
        UTR5_AND_START,
        Change("gggt[caATG>]GATG"),
        {
            **variant(per_strand("TCAATG", "CCATTG"), per_strand("T", "C"), per_strand(13, 58), per_strand(19, 64)),
            "alt_cds_start": per_strand(16, 15),
            "alt_cds_stop": per_strand(63, 62),
            "alt_cds_seq": "GATGCCCTA" + "AAA" * 3 + "GACTAA",
            "alt_cds_len": 24,
            "alt_cds_info": [(1, 13), (2, 11)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 21, 1, [(21, "TAA")], [2]),
            "alt_is_premature": True,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGTGATGCCCTA" + "AAA" * 3 + "GACTAAGGGGG",
            "alt_transcript_length": 33,
            "alt_cds_start_in_transcript": 4,
            **columns(SCAN_COLUMNS, 5, 1, "GGG", False, "TAA", 11, 1, [(11, "TAA")], [1]),
            "unknown_reason": None,
            **columns(PTC_FEATURE_COLUMNS, 0, 1, 6, True, 17, 6),
            "stop_codon_distance": 14,
            **columns(RULE_COLUMNS, False, True, False, True, False, True),
            "alt_transcript_exon_info": [(1, 17), (2, 16)],
        },
        equivalent=(Change("gggt[caATGG>G]ATG"),),
        marks=(
            Ruler((0, 4, 11, 17, 25, 33), "tx", "alt"),
            Mark("alt", 5, 8, "a", "ATG of the scan"),
            Mark("alt", 11, 14, "*", "PTC"),
            Mark("alt", 25, 28, "s", "the annotated stop codon"),
            Span("alt", 5, 11, "ptc_to_start_codon = 6"),
            Span("alt", 11, 17, "ptc_to_intron = 6"),
        ),
    ),
]

# SC-23, SC-24: Ensembl GFF3. It has no start_codon rows, so a CDS has a start codon only if it starts with ATG in
# phase 0. This is the CDS of ATG_START with the non-ATG start CTG.
ENSEMBL_CTG_START = Layout(
    Transcript(("gggCTGAA", "ACCCGACTAAggggg"), flavor="ensembl"),
    {
        **CTG_START.ref,
        "has_start_codon": False,
        "ref_start_codon_pos": None,
        "ref_start_codon_exon": None,
        "likely_misannotated": True,
    },
)

# SC-23: the CDS A TGC AAA CCC GAC TAA starts with ATG, but its phase is 1
ENSEMBL_ATG_IN_PHASE_1 = Layout(
    Transcript(("gggATGCA", "AACCCGACTAAggggg"), frame=1, start_codon=False, flavor="ensembl", tags=("cds_start_NF",)),
    {
        **ATG_START.ref,
        "ref_cds_stop": per_strand(49, 51),
        "ref_cds_seq": "ATGCAAACCCGACTAA",
        "ref_cds_len": 16,
        "has_start_codon": False,
        "cds_frame": 1,
        "ref_cds_info": [(1, 5), (2, 11)],
        "ref_start_codon_pos": None,
        "ref_start_codon_exon": None,
        "ref_first_stop_pos": 13,
        "ref_all_stop_codons": [(13, "TAA")],
        "transcript_end": 54,
        "transcript_seq": "GGGATGCAAACCCGACTAAGGGGG",
        "transcript_length": 24,
        "cds_end_in_transcript": 19,
        "transcript_exon_info": [(1, 8), (2, 16)],
        "likely_misannotated": True,
    },
)

# SC-23: the CDS starts with ATG in phase 0, and an intron splits the ATG into AT and G
ENSEMBL_SPLIT_ATG = Layout(
    Transcript(("gggAT", "GAAACCCGACTAAggggg"), flavor="ensembl"),
    {
        **ATG_START.ref,
        "ref_cds_info": [(1, 2), (2, 13)],
        "transcript_exon_info": [(1, 5), (2, 18)],
    },
)

# SC-24
ENSEMBL_ATG_START = Layout(Transcript(("gggATGAA", "ACCCGACTAAggggg"), flavor="ensembl"), ATG_START.ref)

CASES += [
    # SC-23
    Case(
        "ensembl_cds_starting_with_ctg_has_no_start_codon_and_a_missense_is_no_start_loss",
        """
        tx      0   3        8         15  18
        ref 5' [ggg CTG AA]|[A CCC GAC TAA ggggg] 3'
        alt 5' [ggg CTG AA]|[A CGC GAC TAA ggggg] 3'
                                ^ C>G
        Ensembl GFF3, no start_codon rows
        An Ensembl CDS has a start codon only if it starts with ATG. CTG is not found: has_start_codon = False, both
        start codon positions are null, likely_misannotated = True. The missense CCC>CGC is no start loss.
        """,
        ENSEMBL_CTG_START,
        Change("AC[C>G]CGAC"),
        {
            **variant(per_strand("C", "G"), per_strand("G", "C"), per_strand(40, 22), per_strand(41, 23)),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(48, 50),
            "alt_cds_seq": "CTGAAACGCGACTAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 5), (2, 10)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 12, 1, [(12, "TAA")], [2]),
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GGGCTGAAACGCGACTAAGGGGG",
            "alt_transcript_length": 23,
            "alt_cds_start_in_transcript": 3,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        ruler=Ruler((0, 3, 8, 15, 18)),
    ),
    # SC-23
    Case(
        "ensembl_cds_starting_with_atg_in_phase_1_has_no_start_codon",
        """
        tx      0   3         8          16  19
        ref 5' [ggg A TGC A]|[AA CCC GAC TAA ggggg] 3'
        alt 5' [ggg A TGC A]|[AA CGC GAC TAA ggggg] 3'
                                  ^ C>G
                                         sss the annotated stop codon
        Ensembl GFF3, phase 1 of the 5'-most CDS row
        The CDS starts with ATG, but its phase is 1, so the first complete codon is the TGC at CDS 1. An Ensembl CDS
        has a start codon only if it starts with ATG in phase 0: has_start_codon = False. The frame 1 reads to the TAA
        at CDS 13. The missense CCC>CGC is no start loss.
        """,
        ENSEMBL_ATG_IN_PHASE_1,
        Change("AC[C>G]CGAC"),
        {
            **variant(per_strand("C", "G"), per_strand("G", "C"), per_strand(41, 22), per_strand(42, 23)),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(49, 51),
            "alt_cds_seq": "ATGCAAACGCGACTAA",
            "alt_cds_len": 16,
            "alt_cds_info": [(1, 5), (2, 11)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 13, 1, [(13, "TAA")], [2]),
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GGGATGCAAACGCGACTAAGGGGG",
            "alt_transcript_length": 24,
            "alt_cds_start_in_transcript": 3,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(Mark("alt", 16, 19, "s", "the annotated stop codon"),),
        ruler=Ruler((0, 3, 8, 16, 19)),
    ),
    # SC-23
    Case(
        "ensembl_cds_starting_with_an_atg_split_by_an_intron_has_a_start_codon",
        """
        tx      0   3    5             15  18
        ref 5' [ggg AT]|[G AAA CCC GAC TAA ggggg] 3'
        alt 5' [ggg AT]|[G AAA CGC GAC TAA ggggg] 3'
                                ^ C>G
        Ensembl GFF3
        The CDS starts with ATG in phase 0, split into AT and G by an intron. An Ensembl CDS has a start codon if it
        starts with ATG in phase 0, wherever the exon junction lies: has_start_codon = True, in exon 1. The missense
        CCC>CGC keeps the start codon.
        """,
        ENSEMBL_SPLIT_ATG,
        Change("AC[C>G]CGAC"),
        {
            **variant(per_strand("C", "G"), per_strand("G", "C"), per_strand(40, 22), per_strand(41, 23)),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(48, 50),
            "alt_cds_seq": "ATGAAACGCGACTAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 2), (2, 13)],
            **columns(ALT_CODON_COLUMNS, 0, 1, "TAA", True, "TAA", 12, 1, [(12, "TAA")], [2]),
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GGGATGAAACGCGACTAAGGGGG",
            "alt_transcript_length": 23,
            "alt_cds_start_in_transcript": 3,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        ruler=Ruler((0, 3, 5, 15, 18)),
    ),
    # SC-24
    Case(
        "ensembl_atg_to_acg_in_a_leading_atg_in_phase_0_is_a_start_loss",
        """
        tx      0   3        8         15  18
        ref 5' [ggg ATG AA]|[A CCC GAC TAA ggggg] 3'
        alt 5' [ggg ACG AA]|[A CCC GAC TAA ggggg] 3'
                     ^ T>C
        Ensembl GFF3, no start_codon rows
        The Ensembl CDS starts with ATG in phase 0, so ATG>ACG is a start loss, as with the start_codon rows of a GENCODE
        GFF3. From tx 3 on, the alt transcript has no ATG, so the scan reads no stop codon.
        """,
        ENSEMBL_ATG_START,
        Change("gggA[T>C]GAA"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(14, 48), per_strand(15, 49)),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(48, 50),
            "alt_cds_seq": "ACGAAACCCGACTAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 5), (2, 10)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 12, 1, [(12, "TAA")], [2]),
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGACGAAACCCGACTAAGGGGG",
            "alt_transcript_length": 23,
            "alt_cds_start_in_transcript": 3,
            **NO_ATG_FOUND,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        ruler=Ruler((0, 3, 8, 15, 18)),
    ),
]

# The reference values that the cds_start_NF layouts of FR-01 to FR-10 share: no start_codon rows, so no start codon
CDS_START_NF_REF = {
    **IDS,
    "has_start_codon": False,
    "has_stop_codon": True,
    "cds_in_transcript": True,
    "ref_start_codon_pos": None,
    "ref_start_codon_exon": None,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_num_stop_codons": 1,
    "ref_is_premature": False,
    "transcript_start": 10,
    "cds_start_in_transcript": 3,
    "utr5_length": 3,
    "total_exon_count": 2,
    "likely_misannotated": True,
}


def exon_2_cds_row_with_phase_0(rows, strand):
    """Set the phase of the exon 2 CDS row to 0. The 5'-most CDS row alone gives cds_frame."""
    for row in rows:
        if row[2] == "CDS" and row[8].startswith("ID=CDS:tx1:2;"):
            row[7] = "0"
    return rows


def exon_1_cds_row_with_phase_dot(rows, strand):
    """Set the phase of the exon 1 CDS row, the 5'-most one, to ".". GFF3 requires the phase 0, 1 or 2 on a CDS row."""
    for row in rows:
        if row[2] == "CDS" and row[8].startswith("ID=CDS:tx1:1;"):
            row[7] = "."
    return rows


def frame_layout(frame, edit_gff3=None):
    """
    The CDS of 15 nt TGC AAA CCC CAA GGC and the stop codon TAA, with `frame` extra bases AC before the first complete
    codon. The transcript has no start_codon rows and the tag cds_start_NF. All positions after the CDS start shift
    by the frame.
    """
    cds = "AC"[:frame] + "TGCAAACCCCAAGGC"
    return Layout(
        Transcript(
            ("ggg" + cds[:5], cds[5:] + "TAA" + "ggggg"),
            frame=frame,
            start_codon=False,
            tags=("cds_start_NF",),
            edit_gff3=edit_gff3,
        ),
        {
            **CDS_START_NF_REF,
            "ref_cds_start": per_strand(13, 15),
            "ref_cds_stop": per_strand(51 + frame, 53 + frame),
            "ref_cds_seq": cds + "TAA",
            "ref_cds_len": 18 + frame,
            "cds_frame": frame,
            "ref_cds_info": [(1, 5), (2, 13 + frame)],
            "ref_first_stop_pos": 15 + frame,
            "ref_all_stop_codons": [(15 + frame, "TAA")],
            "ref_stop_codon_exons": [2],
            "transcript_end": 56 + frame,
            "transcript_seq": "GGG" + cds + "TAAGGGGG",
            "transcript_length": 26 + frame,
            "cds_end_in_transcript": 21 + frame,
            "transcript_exon_info": [(1, 8), (2, 18 + frame)],
            "utr3_length": 5,
        },
    )


def frame_snv_row(frame):
    """The row of TGC>TGT: a TAA out of frame at CDS 2 + frame, and the annotated TAA still in frame."""
    cds = "AC"[:frame] + "TGTAAACCCCAAGGC" + "TAA"
    return {
        **variant(per_strand("C", "G"), per_strand("T", "A"), per_strand(15 + frame, 50), per_strand(16 + frame, 51)),
        "alt_cds_start": per_strand(13, 15),
        "alt_cds_stop": per_strand(51 + frame, 53 + frame),
        "alt_cds_seq": cds,
        "alt_cds_len": 18 + frame,
        "alt_cds_info": [(1, 5), (2, 13 + frame)],
        **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 15 + frame, 1, [(15 + frame, "TAA")], [2]),
        "alt_is_premature": False,
        "start_loss": False,
        "stop_loss": False,
        "alt_transcript_seq": "GGG" + cds + "GGGGG",
        "alt_transcript_length": 26 + frame,
        "alt_cds_start_in_transcript": 3,
        **NOT_SCANNED,
        "unknown_reason": None,
        **NO_PTC_FEATURES,
        "stop_codon_distance": 0,
        **NO_RULE,
    }


def frame_ptc_row(frame):
    """The row of CAA>TAA: a PTC at CDS 9 + frame, 6 nt before the annotated stop codon, in the last exon."""
    cds = "AC"[:frame] + "TGCAAACCCTAAGGC" + "TAA"
    return {
        **variant(per_strand("C", "G"), per_strand("T", "A"), per_strand(42 + frame, 23), per_strand(43 + frame, 24)),
        "alt_cds_start": per_strand(13, 15),
        "alt_cds_stop": per_strand(51 + frame, 53 + frame),
        "alt_cds_seq": cds,
        "alt_cds_len": 18 + frame,
        "alt_cds_info": [(1, 5), (2, 13 + frame)],
        **columns(
            ALT_CODON_COLUMNS,
            None,
            None,
            "TAA",
            True,
            "TAA",
            9 + frame,
            2,
            [(9 + frame, "TAA"), (15 + frame, "TAA")],
            [2, 2],
        ),
        "alt_is_premature": True,
        "start_loss": False,
        "stop_loss": False,
        "alt_transcript_seq": "GGG" + cds + "GGGGG",
        "alt_transcript_length": 26 + frame,
        "alt_cds_start_in_transcript": 3,
        **NOT_SCANNED,
        "unknown_reason": None,
        **columns(PTC_FEATURE_COLUMNS, 1, 0, None, False, 18 + frame, 14),
        "stop_codon_distance": 6,
        **columns(RULE_COLUMNS, True, False, False, False, False, True),
    }


CASES += [
    # FR-01, the SNV that makes a TAA out of frame
    Case(
        "cds_frame_0_snv_that_makes_a_taa_out_of_frame_is_neither_a_ptc_nor_a_stop_loss",
        """
        tx      0   3        8     12      18  21
        ref 5' [ggg TGC AA]|[A CCC CAA GGC TAA ggggg] 3'
        alt 5' [ggg TGT AA]|[A CCC CAA GGC TAA ggggg] 3'
                      ^ C>T
        cds_frame 0, cds_start_NF
        TGC>TGT puts a TAA at CDS 2, out of frame: the T of TGT and the AA after it. The frame 0 scan ignores it. The
        annotated TAA at CDS 15 is the first in-frame stop codon, so the row is neither a PTC nor a stop loss.
        """,
        frame_layout(0),
        Change("gggTG[C>T]AA"),
        {**frame_snv_row(0), "alt_transcript_exon_info": SAME_EXONS},
        ruler=Ruler((0, 3, 8, 12, 18, 21)),
    ),
    # FR-01, the SNV that makes a PTC in frame
    Case(
        "cds_frame_0_nonsense_snv_in_the_last_exon_is_a_ptc_6_nt_before_the_annotated_stop_codon",
        """
        tx      0   3        8     12      18  21
        ref 5' [ggg TGC AA]|[A CCC CAA GGC TAA ggggg] 3'
        alt 5' [ggg TGC AA]|[A CCC TAA GGC TAA ggggg] 3'
                                   ^ C>T
                                   *** PTC
                                           sss the annotated stop codon
                                   <---------------> ptc_to_intron = 14
        cds_frame 0, cds_start_NF
        CAA>TAA is a PTC at CDS 9. The annotated TAA at CDS 15 follows 6 nt later. The PTC lies in the last exon, so
        the last exon rule applies and NMD escapes. ptc_to_intron = 14 is the distance to the end of the last exon. The
        transcript has no annotated start codon: ptc_to_start_codon = null.
        """,
        frame_layout(0),
        Change("CCC[C>T]AAGG"),
        {**frame_ptc_row(0), "alt_transcript_exon_info": SAME_EXONS},
        marks=(
            Mark("alt", 12, 15, "*", "PTC"),
            Mark("alt", 18, 21, "s", "the annotated stop codon"),
            Span("alt", 12, 26, "ptc_to_intron = 14"),
        ),
        ruler=Ruler((0, 3, 8, 12, 18, 21)),
    ),
    # FR-02, FR-04 (the minus strand: the CDS rows have the phases 1 and 2, and exon 1 is the row with the largest Start)
    Case(
        "cds_frame_1_snv_that_makes_a_taa_out_of_frame_is_neither_a_ptc_nor_a_stop_loss",
        """
        tx      0   3         8      13      19  22
        ref 5' [ggg A TGC A]|[AA CCC CAA GGC TAA ggggg] 3'
        alt 5' [ggg A TGT A]|[AA CCC CAA GGC TAA ggggg] 3'
                        ^ C>T
        cds_frame 1, phases 1 and 2 of the CDS rows
        The A at CDS 0 belongs to no codon. TGC>TGT puts a TAA at CDS 3, out of frame. The frame 1 scan ignores it.
        The annotated TAA at CDS 16 is the first in-frame stop codon.
        """,
        frame_layout(1),
        Change("gggATG[C>T]A"),
        {**frame_snv_row(1), "alt_transcript_exon_info": SAME_EXONS},
        ruler=Ruler((0, 3, 8, 13, 19, 22)),
    ),
    # FR-02, FR-04
    Case(
        "cds_frame_1_nonsense_snv_in_the_last_exon_is_a_ptc_6_nt_before_the_annotated_stop_codon",
        """
        tx      0   3         8      13      19  22
        ref 5' [ggg A TGC A]|[AA CCC CAA GGC TAA ggggg] 3'
        alt 5' [ggg A TGC A]|[AA CCC TAA GGC TAA ggggg] 3'
                                     ^ C>T
                                     *** PTC
                                             sss the annotated stop codon
                                     <---------------> ptc_to_intron = 14
        cds_frame 1, phases 1 and 2 of the CDS rows
        The frame 1 codons start at CDS 1. CAA>TAA is a PTC at CDS 10. The annotated TAA at CDS 16 follows 6 nt
        later. ptc_to_intron = 14 is the distance to the end of the last exon.
        """,
        frame_layout(1),
        Change("CCC[C>T]AAGG"),
        {**frame_ptc_row(1), "alt_transcript_exon_info": SAME_EXONS},
        marks=(
            Mark("alt", 13, 16, "*", "PTC"),
            Mark("alt", 19, 22, "s", "the annotated stop codon"),
            Span("alt", 13, 27, "ptc_to_intron = 14"),
        ),
        ruler=Ruler((0, 3, 8, 13, 19, 22)),
    ),
    # FR-03, FR-04 (phases 2 and 0)
    Case(
        "cds_frame_2_snv_that_makes_a_taa_out_of_frame_is_neither_a_ptc_nor_a_stop_loss",
        """
        tx      0   3        8       14      20  23
        ref 5' [ggg AC TGC]|[AAA CCC CAA GGC TAA ggggg] 3'
        alt 5' [ggg AC TGT]|[AAA CCC CAA GGC TAA ggggg] 3'
                         ^ C>T
        cds_frame 2, phases 2 and 0 of the CDS rows
        The AC at CDS 0 and 1 belong to no codon. TGC>TGT puts a TAA at CDS 4, out of frame. The frame 2 scan ignores
        it. The annotated TAA at CDS 17 is the first in-frame stop codon.
        """,
        frame_layout(2),
        Change("gggACTG[C>T]"),
        {**frame_snv_row(2), "alt_transcript_exon_info": SAME_EXONS},
        ruler=Ruler((0, 3, 8, 14, 20, 23)),
    ),
    # FR-03, FR-04
    Case(
        "cds_frame_2_nonsense_snv_in_the_last_exon_is_a_ptc_6_nt_before_the_annotated_stop_codon",
        """
        tx      0   3        8       14      20  23
        ref 5' [ggg AC TGC]|[AAA CCC CAA GGC TAA ggggg] 3'
        alt 5' [ggg AC TGC]|[AAA CCC TAA GGC TAA ggggg] 3'
                                     ^ C>T
                                     *** PTC
                                             sss the annotated stop codon
                                     <---------------> ptc_to_intron = 14
        cds_frame 2, phases 2 and 0 of the CDS rows
        The frame 2 codons start at CDS 2. CAA>TAA is a PTC at CDS 11. The annotated TAA at CDS 17 follows 6 nt
        later. ptc_to_intron = 14 is the distance to the end of the last exon.
        """,
        frame_layout(2),
        Change("CCC[C>T]AAGG"),
        {**frame_ptc_row(2), "alt_transcript_exon_info": SAME_EXONS},
        marks=(
            Mark("alt", 14, 17, "*", "PTC"),
            Mark("alt", 20, 23, "s", "the annotated stop codon"),
            Span("alt", 14, 28, "ptc_to_intron = 14"),
        ),
        ruler=Ruler((0, 3, 8, 14, 20, 23)),
    ),
    # FR-08
    Case(
        "cds_frame_comes_from_the_5prime_most_cds_row_and_not_from_the_phase_of_the_other_rows",
        """
        tx      0   3         8      13      19  22
        ref 5' [ggg A TGC A]|[AA CCC CAA GGC TAA ggggg] 3'
        alt 5' [ggg A TGC A]|[AA CCC TAA GGC TAA ggggg] 3'
                                     ^ C>T
                                     *** PTC
                                             sss the annotated stop codon
                                     <---------------> ptc_to_intron = 14
        phase 1 on the exon 1 CDS row, 0 on the exon 2 CDS row
        The GFF3 gives the exon 2 CDS row the phase 0 instead of 2. Only the 5'-most row counts, so cds_frame = 1 and
        the row is the PTC row of the case above.
        """,
        frame_layout(1, edit_gff3=exon_2_cds_row_with_phase_0),
        Change("CCC[C>T]AAGG"),
        {**frame_ptc_row(1), "alt_transcript_exon_info": SAME_EXONS},
        marks=(
            Mark("alt", 13, 16, "*", "PTC"),
            Mark("alt", 19, 22, "s", "the annotated stop codon"),
            Span("alt", 13, 27, "ptc_to_intron = 14"),
        ),
        ruler=Ruler((0, 3, 8, 13, 19, 22)),
    ),
    # FR-10
    Case(
        "cds_frame_1_change_of_the_base_before_the_first_complete_codon_changes_no_codon",
        """
        tx      0   3         8      13      19  22
        ref 5' [ggg A TGC A]|[AA CCC CAA GGC TAA ggggg] 3'
        alt 5' [ggg G TGC A]|[AA CCC CAA GGC TAA ggggg] 3'
                    ^ A>G
        cds_frame 1, cds_start_NF
        The A at CDS 0 belongs to no codon, so A>G changes no codon. The transcript has no annotated start codon, so
        it is no start loss. The ref and the alt codon scans start at CDS 1.
        """,
        frame_layout(1),
        Change("ggg[A>G]TGCA"),
        {
            **variant(per_strand("A", "T"), per_strand("G", "C"), per_strand(13, 53), per_strand(14, 54)),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(52, 54),
            "alt_cds_seq": "GTGCAAACCCCAAGGCTAA",
            "alt_cds_len": 19,
            "alt_cds_info": [(1, 5), (2, 14)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 16, 1, [(16, "TAA")], [2]),
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GGGGTGCAAACCCCAAGGCTAAGGGGG",
            "alt_transcript_length": 27,
            "alt_cds_start_in_transcript": 3,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        ruler=Ruler((0, 3, 8, 13, 19, 22)),
    ),
    # FR-11
    Case(
        "cds_row_with_phase_dot_is_an_error_that_names_the_transcript",
        """
        ref 5' [ggg ATG AA]|[A CCC GAC TAA ggggg] 3'
        alt 5' [ggg ATG AA]|[A CCC AAC TAA ggggg] 3'
                                   ^ G>A
        The exon 1 CDS row is the 5'-most one, which gives cds_frame, and it has the phase "." instead of 0.
        annotate() rejects the GFF3 with a ValueError that names the transcript tx1.
        """,
        Layout(Transcript(ATG_START.transcript.exons, edit_gff3=exon_1_cds_row_with_phase_dot), ATG_START.ref),
        Change("ACCC[G>A]ACTAA"),
        Raises(ValueError, match="tx1"),
    ),
]

# FR-05: the CDS A TGC AAA CCC TAA of frame 1, and a 3' UTR with a TGA in the CDS frame at tx 25
FRAME_1_STOP_LOSS = Layout(
    Transcript(("gggATGCA", "AACCCTAAgggcccgggtgacc"), frame=1, start_codon=False, tags=("cds_start_NF",)),
    {
        **CDS_START_NF_REF,
        "ref_cds_start": per_strand(13, 24),
        "ref_cds_stop": per_strand(46, 57),
        "ref_cds_seq": "ATGCAAACCCTAA",
        "ref_cds_len": 13,
        "cds_frame": 1,
        "ref_cds_info": [(1, 5), (2, 8)],
        "ref_first_stop_pos": 10,
        "ref_all_stop_codons": [(10, "TAA")],
        "ref_stop_codon_exons": [2],
        "transcript_end": 60,
        "transcript_seq": "GGGATGCAAACCCTAAGGGCCCGGGTGACC",
        "transcript_length": 30,
        "cds_end_in_transcript": 16,
        "transcript_exon_info": [(1, 8), (2, 22)],
        "utr3_length": 14,
    },
)

# FR-06: the CDS CTG AAA CCC GAC TAA of frame 0 without a start codon, and a TAG in the 3' UTR at tx 21
NO_START_STOP_LOSS = Layout(
    Transcript(("gggCTGAA", "ACCCGACTAAgggtagccgg"), start_codon=False, tags=("cds_start_NF",)),
    {
        **CDS_START_NF_REF,
        "ref_cds_start": per_strand(13, 20),
        "ref_cds_stop": per_strand(48, 55),
        "ref_cds_seq": "CTGAAACCCGACTAA",
        "ref_cds_len": 15,
        "cds_frame": 0,
        "ref_cds_info": [(1, 5), (2, 10)],
        "ref_first_stop_pos": 12,
        "ref_all_stop_codons": [(12, "TAA")],
        "ref_stop_codon_exons": [2],
        "transcript_end": 58,
        "transcript_seq": "GGGCTGAAACCCGACTAAGGGTAGCCGG",
        "transcript_length": 28,
        "cds_end_in_transcript": 18,
        "transcript_exon_info": [(1, 8), (2, 20)],
        "utr3_length": 10,
    },
)

# FR-07: a CDS of 15 nt with the phase 1, so the annotated TAA at CDS 12 is out of frame
STOP_OUT_OF_FRAME = Layout(
    Transcript(("gggGCCAA", "ACCCGACTAAggggg"), frame=1, start_codon=False, tags=("cds_start_NF",)),
    {
        **CDS_START_NF_REF,
        "ref_cds_start": per_strand(13, 15),
        "ref_cds_stop": per_strand(48, 50),
        "ref_cds_seq": "GCCAAACCCGACTAA",
        "ref_cds_len": 15,
        "cds_frame": 1,
        "ref_cds_info": [(1, 5), (2, 10)],
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "transcript_end": 53,
        "transcript_seq": "GGGGCCAAACCCGACTAAGGGGG",
        "transcript_length": 23,
        "cds_end_in_transcript": 18,
        "transcript_exon_info": [(1, 8), (2, 15)],
        "utr3_length": 5,
    },
)

# FR-09: the annotated start codon ACG, the CDS frame 1 (CDS A CGC ATG CCC GAC TAA), and an ATG at CDS 4
START_CODON_AND_FRAME_1 = Layout(
    Transcript(("gggACGCATGC", "CCGACTAAggggg"), frame=1),
    {
        **ATG_TAA_BASE,
        "ref_cds_stop": per_strand(49, 51),
        "ref_cds_seq": "ACGCATGCCCGACTAA",
        "ref_cds_len": 16,
        "cds_frame": 1,
        "ref_cds_info": [(1, 8), (2, 8)],
        "ref_first_stop_pos": 13,
        "ref_all_stop_codons": [(13, "TAA")],
        "ref_stop_codon_exons": [2],
        "transcript_end": 54,
        "transcript_seq": "GGGACGCATGCCCGACTAAGGGGG",
        "transcript_length": 24,
        "cds_end_in_transcript": 19,
        "transcript_exon_info": [(1, 11), (2, 13)],
        "total_exon_count": 2,
    },
)

CASES += [
    # FR-05
    Case(
        "cds_frame_1_stop_loss_reads_through_the_3utr_in_the_cds_frame_to_a_tga",
        """
        tx      0   3         8      13  16       25
        ref 5' [ggg A TGC A]|[AA CCC TAA gggcccgggtgacc] 3'
        alt 5' [ggg A TGC A]|[AA CCC CAA gggcccgggtgacc] 3'
                                     ^ T>C
                                     sss the annotated stop codon
                                                  sss the stop codon of the scan
        cds_frame 1, cds_start_NF
        TAA>CAA loses the annotated stop codon at tx 13. The scan starts at the first complete codon, tx 4, so it reads
        the codons at 13, 16, 19, 22 and reaches the TGA at tx 25 in the 3' UTR. The transcript has no annotated start
        codon, so transcript_start_codon_pos = null. stop_codon_distance = 13 - 25 = -12.
        """,
        FRAME_1_STOP_LOSS,
        Change("CCC[T>C]AAggg"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(43, 26), per_strand(44, 27)),
            "alt_cds_start": per_strand(13, 24),
            "alt_cds_stop": per_strand(46, 57),
            "alt_cds_seq": "ATGCAAACCCCAA",
            "alt_cds_len": 13,
            "alt_cds_info": [(1, 5), (2, 8)],
            **columns(ALT_CODON_COLUMNS, None, None, "CAA", False, None, None, 0, [], []),
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": "GGGATGCAAACCCCAAGGGCCCGGGTGACC",
            "alt_transcript_length": 30,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, None, None, "ACC", False, "TGA", 25, 1, [(25, "TGA")], [2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": -12,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 13, 16, "s", "the annotated stop codon"),
            Mark("alt", 25, 28, "s", "the stop codon of the scan"),
        ),
        ruler=Ruler((0, 3, 8, 13, 16, 25)),
    ),
    # FR-06
    Case(
        "stop_loss_in_a_cds_without_a_start_codon_reads_through_the_3utr_to_a_tag",
        """
        tx      0   3        8         15  18 21
        ref 5' [ggg CTG AA]|[A CCC GAC TAA gggtagccgg] 3'
        alt 5' [ggg CTG AA]|[A CCC GAC CAA gggtagccgg] 3'
                                       ^ T>C
                                       sss the annotated stop codon
                                              sss the stop codon of the scan
        cds_frame 0, cds_start_NF
        TAA>CAA loses the annotated stop codon at tx 15. The CDS has no ATG at its start and no annotated start codon.
        The scan reads in the CDS frame from tx 3 and reaches the TAG at tx 21 in the 3' UTR.
        transcript_start_codon_pos = null. stop_codon_distance = 15 - 21 = -6.
        """,
        NO_START_STOP_LOSS,
        Change("GAC[T>C]AAggg"),
        {
            **variant(per_strand("T", "A"), per_strand("C", "G"), per_strand(45, 22), per_strand(46, 23)),
            "alt_cds_start": per_strand(13, 20),
            "alt_cds_stop": per_strand(48, 55),
            "alt_cds_seq": "CTGAAACCCGACCAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 5), (2, 10)],
            **columns(ALT_CODON_COLUMNS, None, None, "CAA", False, None, None, 0, [], []),
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": "GGGCTGAAACCCGACCAAGGGTAGCCGG",
            "alt_transcript_length": 28,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, None, None, "CGG", False, "TAG", 21, 1, [(21, "TAG")], [2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": -6,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(
            Mark("alt", 15, 18, "s", "the annotated stop codon"),
            Mark("alt", 21, 24, "s", "the stop codon of the scan"),
        ),
        ruler=Ruler((0, 3, 8, 15, 18, 21)),
    ),
    # FR-07
    Case(
        "annotated_stop_codon_out_of_frame_in_the_cds_frame_keeps_the_flags_from_the_cds",
        """
        tx      0   3         8        15   18
        ref 5' [ggg G CCA A]|[AC CCG ACT AA ggggg] 3'
        alt 5' [ggg G CCA A]|[AC GCG ACT AA ggggg] 3'
                                 ^ C>G
                                       ssss the annotated TAA, out of frame
        cds_frame 1, cds_start_NF, the CDS has 15 nt
        The CDS of 15 nt does not fit the phase 1: the codons CCA AAC CCG ACT end with AA. The annotated TAA at CDS 12
        is out of frame, so the ref CDS has no in-frame stop codon. The missense CCG>GCG reads no stop codon either.
        The row keeps the flags from the CDS: alt_is_premature = False, stop_loss = False. It is not scanned, and
        stop_codon_distance = null.
        """,
        STOP_OUT_OF_FRAME,
        Change("AC[C>G]CGAC"),
        {
            **variant(per_strand("C", "G"), per_strand("G", "C"), per_strand(40, 22), per_strand(41, 23)),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(48, 50),
            "alt_cds_seq": "GCCAAACGCGACTAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 5), (2, 10)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, None, None, 0, [], []),
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GGGGCCAAACGCGACTAAGGGGG",
            "alt_transcript_length": 23,
            "alt_cds_start_in_transcript": 3,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(Mark("ref", 15, 18, "s", "the annotated TAA, out of frame"),),
        ruler=Ruler((0, 3, 8, 15, 18)),
    ),
    # FR-09
    Case(
        "start_loss_in_a_cds_with_cds_frame_1_scans_from_the_first_complete_codon",
        """
        tx      0   3 4   7       11     16
        ref 5' [ggg A CGC ATG C]|[CC GAC TAA ggggg] 3'
        alt 5' [ggg A TGC ATG C]|[CC GAC TAA ggggg] 3'
                      ^ C>T
                          aaa ATG of the scan
                                         sss the annotated stop codon
        cds_frame 1, the start_codon rows mark ACG
        ACG>ATG changes the annotated start codon: start_loss = True. The scan starts at tx 4, the first complete codon,
        so it skips the new ATG at tx 3 and takes the ATG at tx 7. Its frame is the frame of the annotated stop
        codon: the TAA at tx 16 is the first in-frame stop codon, so the row is neither a PTC nor a stop loss.
        """,
        START_CODON_AND_FRAME_1,
        Change("gggA[C>T]GCATG"),
        {
            **variant(per_strand("C", "G"), per_strand("T", "A"), per_strand(14, 49), per_strand(15, 50)),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(49, 51),
            "alt_cds_seq": "ATGCATGCCCGACTAA",
            "alt_cds_len": 16,
            "alt_cds_info": [(1, 8), (2, 8)],
            **columns(ALT_CODON_COLUMNS, None, None, "TAA", True, "TAA", 13, 1, [(13, "TAA")], [2]),
            "alt_is_premature": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGGATGCATGCCCGACTAAGGGGG",
            "alt_transcript_length": 24,
            "alt_cds_start_in_transcript": 3,
            **columns(SCAN_COLUMNS, 7, 1, "GGG", False, "TAA", 16, 1, [(16, "TAA")], [2]),
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
        },
        marks=(Mark("alt", 7, 10, "a", "ATG of the scan"), Mark("alt", 16, 19, "s", "the annotated stop codon")),
        ruler=Ruler((0, 3, 4, 7, 11, 16)),
    ),
]
