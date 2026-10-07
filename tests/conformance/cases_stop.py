"""
Conformance cases of the stop codon: the stop codon classification (alt_has_ptc, stop_loss and
annotated_stop_distance), the readthrough after a stop loss, the nonstop, indels in and next to the stop codon, rows that
keep the flags from the CDS, transcripts without stop_codon rows, the scan columns of the alt transcript, and codons
split over an exon junction ("Technical Notes.md", section "Stop codon classification").

Most layouts are variants of one transcript of two exons: a 5'UTR of 13 nt, a CDS of 27 nt with the stop codon TAA at
t40, and a 3'UTR of 32 nt. Its 3'UTR holds an in-frame TAG at t52. A 1 nt deletion upstream shifts the frame onto
the TGA at t26 and the TGA at t47, and a 1 nt insertion onto a frame without a stop codon.

The drawings show the transcript 5' to 3', also for the minus strand. The runner renders their layout block from
the case. The line ref is the transcript and the line alt is the transcript with the change. The two lines are
aligned, and `-` fills a gap. `[...]` is an exon and `|` is an exon junction. Upper case is the coding region, in
codons of the annotated frame; lower case is UTR. `..N..` leaves out N bases. `^` marks the change as ref>alt. A mark
line puts a character under bases, and `<-- label -->` spans a length. The ruler tx gives transcript positions,
counted from 0 at the 5' end. Marks and spans on the alt line, and the prose "alt tx", count positions in
`alt_transcript_seq`. `*` marks the PTC, `s` the annotated stop codon at its position in `alt_transcript_seq`, and
`f` a first in-frame stop codon. A few cases mark more: `u` an in-frame stop codon of the ref CDS upstream of the
annotated one, `x` stop_codon rows on a sense codon, and `e` the last base of the alt CDS. A codon mark stands under
the first base of the codon only where the drawing splits the codon: by a gap, an exon junction or the codon spacing
of another frame.
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
    Ruler,
    Span,
    Transcript,
    per_strand,
)

UTR5 = "ccgccgccaccgc"
UTR3 = "cccctgacctag" + "c" * 20
# exon 1 of most layouts: the 5'UTR and the first 12 CDS bases
EXON_1 = UTR5 + "ATGGCCGCCGCC"
# upper case pieces of the expected sequences
REF_UTR5 = "CCGCCGCCACCGC"
REF_UTR3 = "CCCCTGACCTAG" + "C" * 20


# The transcript of most cases: the stop codon TAA at t40, an in-frame TAG at t52 in the 3'UTR
STOP_REF = {
    **IDS,
    "cds_start": per_strand(23, 42),
    "cds_end": per_strand(73, 92),
    "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTAA",
    "ref_cds_length": 30,
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "ref_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
    "cds_in_transcript": True,
    "start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_first_stop_pos": 27,
    "ref_stop_codon_count": 1,
    "ref_stop_codons": [{"position": 27, "codon": "TAA"}],
    "ref_stop_codon_exons": [2],
    "ref_has_ptc": False,
    "transcript_start": per_strand(10, 10),
    "transcript_end": per_strand(105, 105),
    "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTAA" + REF_UTR3,
    "transcript_length": 75,
    "cds_start_in_transcript": 13,
    "cds_end_in_transcript": 43,
    "transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 50}],
    "utr3_length": 32,
    "utr5_length": 13,
    "total_exon_count": 2,
    "likely_misannotated": False,
}
STOP = Layout(Transcript((EXON_1, "CTGACCGCCGCCGCC" + "TAA" + UTR3)), STOP_REF)

# The 3'UTR has no in-frame stop codon
NONSTOP = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCGCC" + "TAA" + "c" * 12)),
    {
        **STOP_REF,
        "cds_start": per_strand(23, 22),
        "cds_end": per_strand(73, 72),
        "transcript_end": per_strand(85, 85),
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTAA" + "C" * 12,
        "transcript_length": 55,
        "transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 30}],
        "utr3_length": 12,
    },
)
# The stop codon ends the transcript: no 3'UTR
NO_UTR3 = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCGCC" + "TAA")),
    {
        **STOP_REF,
        "cds_start": per_strand(23, 10),
        "cds_end": per_strand(73, 60),
        "transcript_end": per_strand(73, 73),
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTAA",
        "transcript_length": 43,
        "transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 18}],
        "utr3_length": 0,
    },
)
# The transcript ends in TGA, out of frame
TGA_END = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCGCC" + "TAA" + "cccctgacctag" + "c" * 17 + "tga")),
    {**STOP_REF, "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTAA" + "CCCCTGACCTAG" + "C" * 17 + "TGA"},
)
# The last sense codon TGG, then the stop codon TAA in the A run TAAA
TGG_TAA_A = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCTGG" + "TAA" + "actg" + "c" * 20)),
    {
        **STOP_REF,
        "cds_start": per_strand(23, 34),
        "cds_end": per_strand(73, 84),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTGGTAA",
        "transcript_end": per_strand(97, 97),
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTGGTAA" + "ACTG" + "C" * 20,
        "transcript_length": 67,
        "transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 42}],
        "utr3_length": 24,
    },
)
# The last sense codon TCC, as the stop codon TAA starts with T
TCC_TAA = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCTCC" + "TAA" + UTR3)),
    {
        **STOP_REF,
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTCCTAA",
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTCCTAA" + REF_UTR3,
    },
)
# The last sense codon TCT
TCT_TAA = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCTCT" + "TAA" + UTR3)),
    {
        **STOP_REF,
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTCTTAA",
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTCTTAA" + REF_UTR3,
    },
)
# The last sense codon TCC, the stop codon TAA and a second TAA in the 3'UTR
TCC_TAA_TAA = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCTCC" + "TAA" + "taacc" + "c" * 20)),
    {
        **STOP_REF,
        "cds_start": per_strand(23, 35),
        "cds_end": per_strand(73, 85),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTCCTAA",
        "transcript_end": per_strand(98, 98),
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTCCTAA" + "TAACC" + "C" * 20,
        "transcript_length": 68,
        "transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 43}],
        "utr3_length": 25,
    },
)
# The stop codon repeat GTA TAG TAG CAT
GTA_TAG_TAG = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCGTA" + "TAG" + "tagcat" + "c" * 20)),
    {
        **STOP_REF,
        "cds_start": per_strand(23, 36),
        "cds_end": per_strand(73, 86),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGTATAG",
        "ref_last_codon": "TAG",
        "ref_first_stop_codon": "TAG",
        "ref_stop_codons": [{"position": 27, "codon": "TAG"}],
        "transcript_end": per_strand(99, 99),
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGTATAG" + "TAGCAT" + "C" * 20,
        "transcript_length": 69,
        "transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 44}],
        "utr3_length": 26,
    },
)
# The last sense codon GCT, then TAA, and a TAA 1 nt into the 3'UTR
GCT_TAA_CTAA = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCGCT" + "TAA" + "ctaagcc" + "c" * 20)),
    {
        **STOP_REF,
        "cds_start": per_strand(23, 37),
        "cds_end": per_strand(73, 87),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCTTAA",
        "transcript_end": per_strand(100, 100),
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCTTAA" + "CTAAGCC" + "C" * 20,
        "transcript_length": 70,
        "transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 45}],
        "utr3_length": 27,
    },
)
# The last sense codon TAC, then the stop codon TAA, and a 3'UTR that starts with g
TAC_TAA_G = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCTAC" + "TAA" + "g" + "c" * 20)),
    {
        **STOP_REF,
        "cds_start": per_strand(23, 31),
        "cds_end": per_strand(73, 81),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTACTAA",
        "transcript_end": per_strand(94, 94),
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTACTAA" + "G" + "C" * 20,
        "transcript_length": 64,
        "transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 39}],
        "utr3_length": 21,
    },
)
# An in-frame TGA at t28, upstream of the annotated stop codon: a selenocysteine codon or a misannotation
INTERNAL_TGA = Layout(
    Transcript((EXON_1, "CTGTGAGCCGCCGCC" + "TAA" + UTR3)),
    {
        **STOP_REF,
        "ref_cds_seq": "ATGGCCGCCGCCCTGTGAGCCGCCGCCTAA",
        "ref_first_stop_codon": "TGA",
        "ref_first_stop_pos": 15,
        "ref_stop_codon_count": 2,
        "ref_stop_codons": [{"position": 15, "codon": "TGA"}, {"position": 27, "codon": "TAA"}],
        "ref_stop_codon_exons": [2, 2],
        "ref_has_ptc": True,
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGTGAGCCGCCGCCTAA" + REF_UTR3,
    },
)
# The stop_codon rows lie on the sense codon TCA
TCA_ROWS = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCGCC" + "TCA" + UTR3)),
    {
        **STOP_REF,
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTCA",
        "ref_last_codon": "TCA",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_stop_codon_count": 0,
        "ref_stop_codons": [],
        "ref_stop_codon_exons": [],
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTCA" + REF_UTR3,
        "likely_misannotated": True,
    },
)
# The CDS has 31 nt. Its stop_codon rows lie on AAC, its last 3 nt. Its in-frame TAA at t40 starts 1 nt before them,
# so ref_has_ptc is True
TAA_C_ROWS = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCGCC" + "TAAC" + UTR3)),
    {
        **STOP_REF,
        "cds_start": per_strand(23, 42),
        "cds_end": per_strand(74, 93),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTAAC",
        "ref_cds_length": 31,
        "ref_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 19}],
        "ref_last_codon": "AAC",
        "ref_valid_stop": False,
        "ref_has_ptc": True,
        "transcript_end": per_strand(106, 106),
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTAAC" + REF_UTR3,
        "transcript_length": 76,
        "cds_end_in_transcript": 44,
        "transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 51}],
        "likely_misannotated": True,
    },
)
# The CDS starts 1 nt after the ATG, with phase 0 and without start_codon rows: the annotated stop codon at t40 is out
# of frame, and the frame of the CDS reads TGA at t26
OUT_OF_FRAME = Layout(
    Transcript((UTR5 + "a" + "TGGCCGCCGCC", "CTGACCGCCGCCGCC" + "TAA" + UTR3), start_codon=False),
    {
        **STOP_REF,
        "cds_start": per_strand(24, 42),
        "cds_end": per_strand(73, 91),
        "ref_cds_seq": "TGGCCGCCGCCCTGACCGCCGCCGCCTAA",
        "ref_cds_length": 29,
        "has_start_codon": False,
        "ref_cds_exons": [{"exon_number": 1, "length": 11}, {"exon_number": 2, "length": 18}],
        "start_codon_exon": None,
        "ref_first_stop_codon": "TGA",
        "ref_first_stop_pos": 12,
        "ref_stop_codons": [{"position": 12, "codon": "TGA"}],
        "ref_has_ptc": True,
        "cds_start_in_transcript": 14,
        "utr5_length": 14,
        "likely_misannotated": True,
    },
)
# No stop_codon rows: the CDS ends with GCC at t39, and the TAA at t40 lies in the 3'UTR
NO_STOP_ROWS_REF = {
    **STOP_REF,
    "cds_end": per_strand(70, 92),
    "cds_start": per_strand(23, 45),
    "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCC",
    "ref_cds_length": 27,
    "has_stop_codon": False,
    "ref_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 15}],
    "ref_last_codon": "GCC",
    "ref_valid_stop": False,
    "ref_first_stop_codon": None,
    "ref_first_stop_pos": None,
    "ref_stop_codon_count": 0,
    "ref_stop_codons": [],
    "ref_stop_codon_exons": [],
    "cds_end_in_transcript": 40,
    "utr3_length": None,
    "likely_misannotated": True,
}
NO_STOP_ROWS = Layout(Transcript((EXON_1, "CTGACCGCCGCCGCC" + "taa" + UTR3), stop_codon=False), NO_STOP_ROWS_REF)
# No stop_codon rows: the CDS ends with TGG at t37, and the TAA at t40 right after it lies in the 3'UTR
NO_STOP_ROWS_TGG = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCTGG" + "taa" + "actg" + "c" * 20), stop_codon=False),
    {
        **NO_STOP_ROWS_REF,
        "cds_start": per_strand(23, 37),
        "cds_end": per_strand(70, 84),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTGG",
        "ref_last_codon": "TGG",
        "transcript_end": per_strand(97, 97),
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTGG" + "TAAACTG" + "C" * 20,
        "transcript_length": 67,
        "transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 42}],
    },
)
# The CDS lies inside exon 2 of 3
INTERNAL_EXON = Layout(
    Transcript(("gcgc", "cacc" + "ATGGCCCAGGCCTAA" + "gcg", "ccgcgc")),
    {
        **IDS,
        "cds_start": per_strand(38, 39),
        "cds_end": per_strand(53, 54),
        "ref_cds_seq": "ATGGCCCAGGCCTAA",
        "ref_cds_length": 15,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [{"exon_number": 2, "length": 15}],
        "cds_in_transcript": True,
        "start_codon_exon": 2,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 12,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [{"position": 12, "codon": "TAA"}],
        "ref_stop_codon_exons": [2],
        "ref_has_ptc": False,
        "transcript_start": per_strand(10, 10),
        "transcript_end": per_strand(82, 82),
        "transcript_seq": "GCGC" + "CACCATGGCCCAGGCCTAAGCG" + "CCGCGC",
        "transcript_length": 32,
        "cds_start_in_transcript": 8,
        "cds_end_in_transcript": 23,
        "transcript_exons": [
            {"exon_number": 1, "length": 4},
            {"exon_number": 2, "length": 22},
            {"exon_number": 3, "length": 6},
        ],
        "utr3_length": 9,
        "utr5_length": 8,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)
# The start codon is split over the exon junction: AT|G
SPLIT_START_REF = {
    **STOP_REF,
    "ref_cds_exons": [{"exon_number": 1, "length": 2}, {"exon_number": 2, "length": 28}],
    "transcript_exons": [{"exon_number": 1, "length": 15}, {"exon_number": 2, "length": 60}],
}
SPLIT_START_EXONS = (UTR5 + "AT", "GGCCGCCGCCCTGACCGCCGCCGCC" + "TAA" + UTR3)
SPLIT_START = Layout(Transcript(SPLIT_START_EXONS), SPLIT_START_REF)
SPLIT_START_ENSEMBL = Layout(Transcript(SPLIT_START_EXONS, flavor="ensembl"), SPLIT_START_REF)
# The stop codon is split over the exon junction: TA|A. Its last base is the only CDS base of exon 3.
SPLIT_STOP = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCGCC" + "TA", "A" + UTR3)),
    {
        **STOP_REF,
        "cds_end": per_strand(93, 112),
        "ref_cds_exons": [
            {"exon_number": 1, "length": 12},
            {"exon_number": 2, "length": 17},
            {"exon_number": 3, "length": 1},
        ],
        "transcript_end": per_strand(125, 125),
        "transcript_exons": [
            {"exon_number": 1, "length": 25},
            {"exon_number": 2, "length": 17},
            {"exon_number": 3, "length": 33},
        ],
        "total_exon_count": 3,
    },
)
# Codon 4 (CDS 12 to 14) is split over the exon junction: C|AG
SPLIT_CODON_1 = Layout(
    Transcript((UTR5 + "ATGGCCGCCGCCC", "AGGCCGCC" + "TAA" + "cccc")),
    {
        **STOP_REF,
        "cds_start": per_strand(23, 14),
        "cds_end": per_strand(67, 58),
        "ref_cds_seq": "ATGGCCGCCGCCCAGGCCGCCTAA",
        "ref_cds_length": 24,
        "ref_cds_exons": [{"exon_number": 1, "length": 13}, {"exon_number": 2, "length": 11}],
        "ref_first_stop_pos": 21,
        "ref_stop_codons": [{"position": 21, "codon": "TAA"}],
        "transcript_end": per_strand(71, 71),
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCAGGCCGCCTAA" + "CCCC",
        "transcript_length": 41,
        "cds_end_in_transcript": 37,
        "transcript_exons": [{"exon_number": 1, "length": 26}, {"exon_number": 2, "length": 15}],
        "utr3_length": 4,
    },
)
# Codon 4 (CDS 12 to 14) is split over the exon junction: TA|C
SPLIT_CODON_2 = Layout(
    Transcript((UTR5 + "ATGGCCGCCGCCTA", "CGCCGCC" + "TAA" + "cccc")),
    {
        **SPLIT_CODON_1.ref,
        "ref_cds_seq": "ATGGCCGCCGCCTACGCCGCCTAA",
        "ref_cds_exons": [{"exon_number": 1, "length": 14}, {"exon_number": 2, "length": 10}],
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCTACGCCGCCTAA" + "CCCC",
        "transcript_exons": [{"exon_number": 1, "length": 27}, {"exon_number": 2, "length": 14}],
    },
)
# The in-frame TAG at t46 in the 3'UTR is split over the exon junction: TA|G
SPLIT_READTHROUGH = Layout(
    Transcript((EXON_1, "CTGACCGCCGCCGCC" + "TAA" + "cccta", "g" + "c" * 20)),
    {
        **STOP_REF,
        "cds_start": per_strand(23, 56),
        "cds_end": per_strand(73, 106),
        "transcript_end": per_strand(119, 119),
        "transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTAA" + "CCCTA" + "G" + "C" * 20,
        "transcript_length": 69,
        "transcript_exons": [
            {"exon_number": 1, "length": 25},
            {"exon_number": 2, "length": 23},
            {"exon_number": 3, "length": 21},
        ],
        "total_exon_count": 3,
        "utr3_length": 26,
    },
)

# 5' [gacc ATG GAT G]|[TA AGC TAA gc] 3'
#     0    4         11           21   tx
#          0          7       12       CDS
# The ATG at CDS 4 is out of frame. It is the ATG that the scan finds after a start loss.
TWO_EXONS_SEQUENCES = ("gaccATGGATG", "TAAGCTAAgc")
TWO_EXONS_CDS = {
    **IDS,
    "cds_start": per_strand(14, 12),
    "cds_end": per_strand(49, 47),
    "ref_cds_seq": "ATGGATGTAAGCTAA",
    "ref_cds_length": 15,
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "ref_cds_exons": [{"exon_number": 1, "length": 7}, {"exon_number": 2, "length": 8}],
    "start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_first_stop_pos": 12,
    "ref_stop_codon_count": 1,
    "ref_stop_codons": [{"position": 12, "codon": "TAA"}],
    "ref_stop_codon_exons": [2],
    "ref_has_ptc": False,
}
TWO_EXONS = Layout(
    Transcript(TWO_EXONS_SEQUENCES),
    {
        **TWO_EXONS_CDS,
        "cds_in_transcript": True,
        "transcript_start": 10,
        "transcript_end": 51,
        "transcript_seq": "GACCATGGATGTAAGCTAAGC",
        "transcript_length": 21,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 19,
        "transcript_exons": [{"exon_number": 1, "length": 11}, {"exon_number": 2, "length": 10}],
        "utr3_length": 2,
        "utr5_length": 4,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)

CASES = [
    # ST-01
    Case(
        "nonsense_snv_in_an_internal_exon_is_a_ptc",
        (
            """
            tx      0      4    8       14      20        26
            ref 5' [gcgc]|[cacc ATG GCC CAG GCC TAA gcg]|[ccgcgc] 3'
            alt 5' [gcgc]|[cacc ATG GCC TAG GCC TAA gcg]|[ccgcgc] 3'
                                        ^ C>T
                                        *** PTC
                                                sss annotated stop codon
                                        <-----> annotated_stop_distance = 6
                                        <-------------> ptc_to_exon_end = 12
            C>T turns CAG into the PTC TAG in exon 2, the middle exon of three.
            """
        ),
        INTERNAL_EXON,
        Change("GCC[C>T]AGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "variant_start": per_strand(44, 47),
            "variant_end": per_strand(45, 48),
            "alt_cds_seq": "ATGGCCTAGGCCTAA",
            "alt_cds_length": 15,
            "alt_cds_exons": [{"exon_number": 2, "length": 15}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 6,
            "alt_stop_codon_count": 2,
            "alt_stop_codons": [{"position": 6, "codon": "TAG"}, {"position": 12, "codon": "TAA"}],
            "alt_stop_codon_exons": [2, 2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GCGC" + "CACCATGGCCTAGGCCTAAGCG" + "CCGCGC",
            "alt_transcript_length": 32,
            "alt_cds_start_in_transcript": 8,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 6,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 22,
            "annotated_stop_distance": 6,
            "ptc_to_exon_end": 12,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": True,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Mark("alt", 14, 17, "*", "PTC"),
            Mark("alt", 20, 23, "s", "annotated stop codon"),
            Span("alt", 14, 20, "annotated_stop_distance = 6"),
            Span("alt", 14, 26, "ptc_to_exon_end = 12"),
        ),
        ruler=Ruler((0, 4, 8, 14, 20, 26)),
    ),
    # ST-02
    Case(
        "stop_codon_snv_taa_to_tag_is_neither_ptc_nor_stop_loss",
        (
            """
            tx      0             13                25                  40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA cccctg..23..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAG cccctg..23..ccc] 3'
                                                                          ^ A>G
            A>G turns the stop codon TAA into the stop codon TAG.
            """
        ),
        STOP,
        Change("GCCTA[A>G]CCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(72, 42),
            "variant_end": per_strand(73, 43),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTAG",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "TAG",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 27,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 27, "codon": "TAG"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTAG" + REF_UTR3,
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 13, 25, 40)),
    ),
    # ST-02
    Case(
        "stop_codon_snv_taa_to_tga_is_neither_ptc_nor_stop_loss",
        (
            """
            tx      0             13                25                  40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA cccct..24..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TGA cccct..24..ccc] 3'
                                                                         ^ A>G
            A>G turns the stop codon TAA into the stop codon TGA.
            """
        ),
        STOP,
        Change("GCCT[A>G]ACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(71, 43),
            "variant_end": per_strand(72, 44),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTGA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "TGA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 27,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 27, "codon": "TGA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTGA" + REF_UTR3,
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 13, 25, 40)),
    ),
    # ST-03
    Case(
        "stop_codon_snv_reads_through_to_an_in_frame_stop_codon_in_the_3utr",
        (
            """
            tx      0             13                25                  40           52
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA cccctgacctagccc..14..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC CAA cccctgacctagccc..14..ccc] 3'
                                                                        ^ T>C
                                                                        sss annotated stop codon
                                                                                     fff first in-frame stop codon
                                                                        <-----------> annotated_stop_distance = -12
            T>C turns the stop codon TAA into CAA. The scan reads on to the in-frame TAG at alt tx 52 in the 3' UTR.
            """
        ),
        STOP,
        Change("GCC[T>C]AACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(70, 44),
            "variant_end": per_strand(71, 45),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "CAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA" + REF_UTR3,
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAG",
            "alt_scan_first_stop_pos": 52,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 52, "codon": "TAG"}],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -12,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        marks=(
            Mark("alt", 40, 43, "s", "annotated stop codon"),
            Mark("alt", 52, 55, "f", "first in-frame stop codon"),
            Span("alt", 40, 52, "annotated_stop_distance = -12"),
        ),
        ruler=Ruler((0, 13, 25, 40, 52)),
    ),
    # ST-04
    Case(
        "stop_codon_snv_without_an_in_frame_stop_codon_up_to_the_transcript_end_is_a_nonstop",
        (
            """
            tx      0             13                25                  40              55
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA cccccccccccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC CAA cccccccccccc] 3'
                                                                        ^ T>C
                                                                        sss annotated stop codon
            T>C turns the stop codon TAA into CAA. No in-frame stop codon follows up to the transcript end.
            """
        ),
        NONSTOP,
        Change("GCC[T>C]AACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(70, 24),
            "variant_end": per_strand(71, 25),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "CAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA" + "C" * 12,
            "alt_transcript_length": 55,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": None,
            "alt_scan_first_stop_pos": None,
            "alt_scan_stop_codon_count": 0,
            "alt_scan_stop_codons": [],
            "alt_scan_stop_codon_exons": [],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        marks=(Mark("alt", 40, 43, "s", "annotated stop codon"),),
        ruler=Ruler((0, 13, 25, 40, 55)),
    ),
    # ST-05
    Case(
        "frameshift_whose_first_stop_codon_lies_inside_the_cds_is_a_ptc",
        (
            """
            tx      0             13    18          25                  40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA ccc..25..cccc] 3'
            alt 5' [ccgccgccaccgc ATG GC- GCC GCC]|[CTG ACC GCC GCC GCC TAA ccc..25..cccc] 3'
                                        ^ C>-
                                                     * PTC TGA, in the alt frame
                                                                        sss annotated stop codon
                                                     <----------------> annotated_stop_distance = 14
                                                     <----------------------------------> ptc_to_exon_end = 49, to the transcript end
            C>- at tx 18 shifts the frame: the alt CDS reads ATG GCG CCG CCC TGA. Its TGA at alt tx 25 is the PTC.
            The PTC lies in the last exon, so ptc_to_exon_end counts to the transcript end.
            """
        ),
        STOP,
        Change("ATGGC[C>]GCCGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CC", "CG"),
            "alt": per_strand("C", "C"),
            "variant_start": per_strand(27, 85),
            "variant_end": per_strand(29, 87),
            "alt_cds_seq": "ATGGCGCCGCCCTGACCGCCGCCGCCTAA",
            "alt_cds_length": 29,
            "alt_cds_exons": [{"exon_number": 1, "length": 11}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 12,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 12, "codon": "TGA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCGCCGCCCTGACCGCCGCCGCCTAA" + REF_UTR3,
            "alt_transcript_length": 74,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 12,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 50,
            "annotated_stop_distance": 14,
            "ptc_to_exon_end": 49,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": [{"exon_number": 1, "length": 24}, {"exon_number": 2, "length": 50}],
            "nmd_model_status": "ok",
        },
        equivalent=(Change("ATGG[C>]CGCCGCC"), Change("ATG[GCC>GC]GCCGCC")),
        marks=(
            Mark("alt", 25, 26, "*", "PTC TGA, in the alt frame"),
            Mark("alt", 39, 42, "s", "annotated stop codon"),
            Span("alt", 25, 39, "annotated_stop_distance = 14"),
            Span("alt", 25, 74, "ptc_to_exon_end = 49, to the transcript end"),
        ),
        ruler=Ruler((0, 13, 18, 25, 40)),
    ),
    # ST-06
    Case(
        "frameshift_whose_first_stop_codon_lies_in_the_3utr_is_a_stop_loss",
        (
            """
            tx      0             13                25        33        40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA cccctgacct..19..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GC- GCC GCC TAA cccctgacct..19..ccc] 3'
                                                              ^ C>-
                                                                        sss annotated stop codon
                                                                                fff first in-frame stop codon
                                                                        <------> annotated_stop_distance = -7
            C>- at tx 33 shifts the frame past the stop codon. The alt CDS has no stop codon, and the scan reads on to
            the TGA at alt tx 46 in the 3' UTR.
            """
        ),
        STOP,
        Change("ACCGC[C>]GCCGCCTAA"),
        {
            "variant_id": "var1",
            "ref": per_strand("CC", "CG"),
            "alt": per_strand("C", "C"),
            "variant_start": per_strand(62, 50),
            "variant_end": per_strand(64, 52),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCGCCGCCTAA",
            "alt_cds_length": 29,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 17}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCGCCGCCTAA" + REF_UTR3,
            "alt_transcript_length": 74,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TGA",
            "alt_scan_first_stop_pos": 46,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 46, "codon": "TGA"}],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -7,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 49}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("ACCG[C>]CGCCGCCTAA"),),
        marks=(
            Mark("alt", 39, 42, "s", "annotated stop codon"),
            Mark("alt", 46, 49, "f", "first in-frame stop codon"),
            Span("alt", 39, 46, "annotated_stop_distance = -7"),
        ),
        ruler=Ruler((0, 13, 25, 33, 40)),
    ),
    # ST-07
    Case(
        "frameshift_without_an_in_frame_stop_codon_up_to_the_transcript_end_is_a_nonstop",
        (
            """
            tx      0             13                25         33        40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GC-C GCC GCC TAA c..28..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCGC GCC GCC TAA c..28..ccc] 3'
                                                              ^ ->G
            G inserted after tx 32 shifts the frame. No in-frame stop codon follows up to the transcript end.
            """
        ),
        STOP,
        Change("ACCGC[>G]CGCCGCCTAA"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("CG", "GC"),
            "variant_start": per_strand(62, 51),
            "variant_end": per_strand(63, 52),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCGCGCCGCCTAA",
            "alt_cds_length": 31,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 19}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCGCGCCGCCTAA" + REF_UTR3,
            "alt_transcript_length": 76,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": None,
            "alt_scan_first_stop_pos": None,
            "alt_scan_stop_codon_count": 0,
            "alt_scan_stop_codons": [],
            "alt_scan_stop_codon_exons": [],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 51}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("ACCGC[C>GC]GCCGCCTAA"),),
        ruler=Ruler((0, 13, 25, 33, 40)),
    ),
    # ST-08
    Case(
        "insertion_inside_the_stop_codon_taa_to_tgaa_keeps_the_stop_codon",
        (
            """
            tx      0             13                25                  40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC T-AA cccct..24..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TGAA cccct..24..ccc] 3'
                                                                         ^ ->G
                                                                        sss annotated stop codon TGA, annotated_stop_distance = 0
            G inserted after tx 40 turns TAA into TGAA. The alt CDS stops at the TGA, at the place of the annotated stop
            codon.
            """
        ),
        STOP,
        Change("GCCT[>G]AACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "T"),
            "alt": per_strand("TG", "TC"),
            "variant_start": per_strand(70, 43),
            "variant_end": per_strand(71, 44),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTGAA",
            "alt_cds_length": 31,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 19}],
            "alt_last_codon": "GAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 27,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 27, "codon": "TGA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTGAA" + REF_UTR3,
            "alt_transcript_length": 76,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 51}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCT[A>GA]ACCCC"),),
        marks=(Mark("alt", 40, 43, "s", "annotated stop codon TGA, annotated_stop_distance = 0"),),
        ruler=Ruler((0, 13, 25, 40)),
    ),
    # ST-09
    Case(
        "deletion_of_one_a_of_the_stop_codon_in_an_a_run_keeps_the_stop_codon",
        (
            """
            tx      0             13                25                  40  43
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TGG TAA actgccc..14..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TGG T-A actgccc..14..ccc] 3'
                                                                         ^ A>-
                                                                        s annotated stop codon TAA, annotated_stop_distance = 0
            A>- deletes one A of the A run TAAA, at tx 41, 42 or 43. Each placement gives the same alt transcript, and
            the stop codon stays TAA.
            """
        ),
        TGG_TAA_A,
        Change("TGGT[A>]AACTG"),
        {
            "variant_id": "var1",
            "ref": per_strand("TA", "TT"),
            "alt": per_strand("T", "T"),
            "variant_start": per_strand(70, 34),
            "variant_end": per_strand(72, 36),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTGGTAA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 27,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 27, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTGGTAA" + "CTG" + "C" * 20,
            "alt_transcript_length": 66,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 41}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("TGGTA[A>]ACTG"), Change("TGGTAA[A>]CTG")),
        marks=(Mark("alt", 40, 41, "s", "annotated stop codon TAA, annotated_stop_distance = 0"),),
        ruler=Ruler((0, 13, 25, 40, 43)),
    ),
    # ST-10
    Case(
        "in_frame_insertion_right_before_the_stop_codon_keeps_the_stop_codon",
        (
            """
            tx      0             13                25   29                40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC ---TAA ccc..26..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC GCCTAA ccc..26..ccc] 3'
                                                                        ^^^ ->GCC
                                                                           sss annotated stop codon, annotated_stop_distance = 0
            GCC inserted after tx 39 keeps the frame. It extends the repeat CCGCCGCCGCC at tx 29 to 39, so CCG inserted
            after tx 28 is the same change.
            """
        ),
        STOP,
        Change("GCCGCCGCC[>GCC]TAACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "A"),
            "alt": per_strand("CGCC", "AGGC"),
            "variant_start": per_strand(69, 44),
            "variant_end": per_strand(70, 45),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCGCCTAA",
            "alt_cds_length": 33,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 21}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 30,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 30, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCGCCTAA" + REF_UTR3,
            "alt_transcript_length": 78,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 53}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GA[>CCG]CCGCCGCCGCCTAA"), Change("GCCGCCGC[C>CGCC]TAACCCC")),
        marks=(Mark("alt", 43, 46, "s", "annotated stop codon, annotated_stop_distance = 0"),),
        ruler=Ruler((0, 13, 25, 29, 40)),
    ),
    # ST-11
    Case(
        "in_frame_deletion_right_before_the_stop_codon_keeps_the_stop_codon",
        (
            """
            tx      0             13                25   29         37  40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA ccc..26..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC --- TAA ccc..26..ccc] 3'
                                                                    ^^^ GCC>-
                                                                        sss annotated stop codon, annotated_stop_distance = 0
            GCC>- at tx 37 to 39 keeps the frame. In the repeat CCGCCGCCGCC at tx 29 to 39, CCG>- at tx 29 to 31 is the
            same change.
            """
        ),
        STOP,
        Change("GCCGCC[GCC>]TAACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CGCC", "AGGC"),
            "alt": per_strand("C", "A"),
            "variant_start": per_strand(66, 44),
            "variant_end": per_strand(70, 48),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTAA",
            "alt_cds_length": 27,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 15}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 24,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 24, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTAA" + REF_UTR3,
            "alt_transcript_length": 72,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 47}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("CTGA[CCG>]CCGCCGCCTAA"),),
        marks=(Mark("alt", 37, 40, "s", "annotated stop codon, annotated_stop_distance = 0"),),
        ruler=Ruler((0, 13, 25, 29, 37, 40)),
    ),
    # ST-12
    Case(
        "tcc_inserted_right_before_the_stop_codon_taa_keeps_the_stop_codon",
        (
            """
            tx      0             13                25              37     40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC ---TAA cccct..24..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TCCTAA cccct..24..ccc] 3'
                                                                        ^^^ ->TCC
                                                                           sss annotated stop codon, annotated_stop_distance = 0
            TCC inserted after tx 39 keeps the frame. CCT inserted after tx 37 or after tx 40 is the same change.
            """
        ),
        STOP,
        Change("GCCGCCGCC[>TCC]TAACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "A"),
            "alt": per_strand("CTCC", "AGGA"),
            "variant_start": per_strand(69, 44),
            "variant_end": per_strand(70, 45),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTCCTAA",
            "alt_cds_length": 33,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 21}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 30,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 30, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTCCTAA" + REF_UTR3,
            "alt_transcript_length": 78,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 53}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("CCGCCG[>CCT]CCTAACC"), Change("GCCGCCGCCT[>CCT]AACCCC")),
        marks=(Mark("alt", 43, 46, "s", "annotated stop codon, annotated_stop_distance = 0"),),
        ruler=Ruler((0, 13, 25, 37, 40)),
    ),
    # ST-12
    Case(
        "tggccc_inserted_right_before_the_stop_codon_taa_keeps_the_stop_codon",
        (
            """
            tx      0             13                25              37        40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC ------TAA cccct..24..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TGGCCCTAA cccct..24..ccc] 3'
                                                                        ^^^^^^ ->TGGCCC
                                                                              sss annotated stop codon, annotated_stop_distance = 0
            TGGCCC inserted after tx 39 keeps the frame. GGCCCT inserted after tx 40 is the same change.
            """
        ),
        STOP,
        Change("GCCGCCGCC[>TGGCCC]TAACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "A"),
            "alt": per_strand("CTGGCCC", "AGGGCCA"),
            "variant_start": per_strand(69, 44),
            "variant_end": per_strand(70, 45),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTGGCCCTAA",
            "alt_cds_length": 36,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 24}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 33,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 33, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTGGCCCTAA" + REF_UTR3,
            "alt_transcript_length": 81,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 56}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCGCCGCCT[>GGCCCT]AACCCC"),),
        marks=(Mark("alt", 46, 49, "s", "annotated stop codon, annotated_stop_distance = 0"),),
        ruler=Ruler((0, 13, 25, 37, 40)),
    ),
    # ST-12
    Case(
        "last_sense_codon_tcc_deleted_before_the_stop_codon_taa_keeps_the_stop_codon",
        (
            """
            tx      0             13                25              37  40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TCC TAA cccc..25..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC --- TAA cccc..25..ccc] 3'
                                                                    ^^^ TCC>-
                                                                        sss annotated stop codon, annotated_stop_distance = 0
            TCC>- at tx 37 to 39 keeps the frame. CCT>- at tx 38 to 40 is the same change.
            """
        ),
        TCC_TAA,
        Change("GCCGCC[TCC>]TAACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CTCC", "AGGA"),
            "alt": per_strand("C", "A"),
            "variant_start": per_strand(66, 44),
            "variant_end": per_strand(70, 48),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTAA",
            "alt_cds_length": 27,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 15}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 24,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 24, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTAA" + REF_UTR3,
            "alt_transcript_length": 72,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 47}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCGCCT[CCT>]AACCCC"),),
        marks=(Mark("alt", 37, 40, "s", "annotated stop codon, annotated_stop_distance = 0"),),
        ruler=Ruler((0, 13, 25, 37, 40)),
    ),
    # ST-13
    Case(
        "stop_codon_gained_right_before_the_stop_codon_is_a_ptc",
        (
            """
            tx      0             13                25                        40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC ------TAA cccct..23..cccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TCCTAGTAA cccct..23..cccc] 3'
                                                                        ^^^^^^ ->TCCTAG
                                                                           *** PTC
                                                                              sss annotated stop codon
                                                                           <-> annotated_stop_distance = 3
                                                                           <--------------------> ptc_to_exon_end = 38, to the transcript end
            TCCTAG inserted after tx 39 puts the PTC TAG right before the stop codon. CCTAGT inserted after tx 40 is the
            same change.
            The PTC lies in the last exon, so ptc_to_exon_end counts to the transcript end.
            """
        ),
        STOP,
        Change("GCCGCCGCC[>TCCTAG]TAACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "A"),
            "alt": per_strand("CTCCTAG", "ACTAGGA"),
            "variant_start": per_strand(69, 44),
            "variant_end": per_strand(70, 45),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTCCTAGTAA",
            "alt_cds_length": 36,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 24}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 30,
            "alt_stop_codon_count": 2,
            "alt_stop_codons": [{"position": 30, "codon": "TAG"}, {"position": 33, "codon": "TAA"}],
            "alt_stop_codon_exons": [2, 2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTCCTAGTAA" + REF_UTR3,
            "alt_transcript_length": 81,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 30,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 56,
            "annotated_stop_distance": 3,
            "ptc_to_exon_end": 38,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 56}],
            "nmd_model_status": "ok",
        },
        equivalent=(Change("GCCGCCGCCT[>CCTAGT]AACCCC"),),
        marks=(
            Mark("alt", 43, 46, "*", "PTC"),
            Mark("alt", 46, 49, "s", "annotated stop codon"),
            Span("alt", 43, 46, "annotated_stop_distance = 3"),
            Span("alt", 43, 81, "ptc_to_exon_end = 38, to the transcript end"),
        ),
        ruler=Ruler((0, 13, 25, 40)),
    ),
    # Pins "It is placed 3'-most, as HGVS does, and also 5'-most if that placement ends at or before the stop codon"
    # ("Stop codon classification") for a 3'-most placement that ends right at the stop codon. Both placements of the
    # delins T>GTGAC end at the stop codon, so the annotated stop codon shifts by +4 to alt tx 44. The TGA at its old
    # position, alt tx 40, is a PTC. In the closest case, stop_codon_gained_right_before_the_stop_codon_is_a_ptc, the
    # 3'-most placement reaches into the stop codon, and the PTC does not lie at the old position.
    Case(
        "delins_before_the_stop_codon_that_puts_a_tga_at_its_old_position_is_a_ptc",
        """
        tx      0             13                25                      40
        ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TCT---- TAA cccc..24..cccc] 3'
        alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TCGTGAC TAA cccc..24..cccc] 3'
                                                                  ^^^^^ T>GTGAC
                                                                   *** PTC TGA at the old position of the stop codon
                                                                        sss annotated stop codon, shifted by +4
                                                                   <--> annotated_stop_distance = 4
                                                                   <---------------------> ptc_to_exon_end = 39, to the transcript end
        T>GTGAC at tx 39, the last base of the last sense codon TCT, gives TCG TGA C TAA. Both placements of the
        delins end at the stop codon, so the annotated stop codon shifts by +4 to alt tx 44. The TGA at alt tx 40,
        the old position of the stop codon, is the PTC.
        The PTC lies in the last exon, so ptc_to_exon_end counts to the transcript end.
        """,
        TCT_TAA,
        Change("GCCTC[T>GTGAC]TAACC"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("GTGAC", "GTCAC"),
            "variant_start": per_strand(69, 45),
            "variant_end": per_strand(70, 46),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTCGTGACTAA",
            "alt_cds_length": 34,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 22}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 27,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 27, "codon": "TGA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTCGTGACTAA" + REF_UTR3,
            "alt_transcript_length": 79,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 27,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 54,
            "annotated_stop_distance": 4,
            "ptc_to_exon_end": 39,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 54}],
            "nmd_model_status": "ok",
        },
        equivalent=(Change("GCCT[CT>CGTGAC]TAACC"), Change("GCCTC[TT>GTGACT]AACC")),
        marks=(
            Mark("alt", 40, 43, "*", "PTC TGA at the old position of the stop codon"),
            Mark("alt", 44, 47, "s", "annotated stop codon, shifted by +4"),
            Span("alt", 40, 44, "annotated_stop_distance = 4"),
            Span("alt", 40, 79, "ptc_to_exon_end = 39, to the transcript end"),
        ),
        ruler=Ruler((0, 13, 25, 40)),
    ),
    # ST-14
    Case(
        "delins_ccta_to_g_in_tcc_taa_leaves_a_tga_at_the_stop_codon",
        (
            """
            tx      0             13                25              37  40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TCC TAA cccct..24..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TG- --A cccct..24..ccc] 3'
                                                                     ^^^^^ CCTA>G
                                                                    s annotated stop codon TGA, annotated_stop_distance = 0
            CCTA>G at tx 38 to 41 leaves T, G and the last A of TAA. They read TGA, a stop codon at the place of the
            annotated one.
            """
        ),
        TCC_TAA,
        Change("GCCGCCT[CCTA>G]ACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CCTA", "TAGG"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(68, 43),
            "variant_end": per_strand(72, 47),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTGA",
            "alt_cds_length": 27,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 15}],
            "alt_last_codon": "TGA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 24,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 24, "codon": "TGA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTGA" + REF_UTR3,
            "alt_transcript_length": 72,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 47}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCGCC[TCCTA>TG]ACCCC"),),
        marks=(Mark("alt", 37, 38, "s", "annotated stop codon TGA, annotated_stop_distance = 0"),),
        ruler=Ruler((0, 13, 25, 37, 40)),
    ),
    # ST-14
    Case(
        "delins_ctta_to_g_in_tct_taa_leaves_a_tga_at_the_stop_codon",
        (
            """
            tx      0             13                25              37  40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TCT TAA cccct..24..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TG- --A cccct..24..ccc] 3'
                                                                     ^^^^^ CTTA>G
                                                                    s annotated stop codon TGA, annotated_stop_distance = 0
            CTTA>G at tx 38 to 41 leaves T, G and the last A of TAA. They read TGA, a stop codon at the place of the
            annotated one.
            """
        ),
        TCT_TAA,
        Change("GCCGCCT[CTTA>G]ACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CTTA", "TAAG"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(68, 43),
            "variant_end": per_strand(72, 47),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTGA",
            "alt_cds_length": 27,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 15}],
            "alt_last_codon": "TGA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 24,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 24, "codon": "TGA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTGA" + REF_UTR3,
            "alt_transcript_length": 72,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 47}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCGCC[TCTTA>TG]ACCCC"),),
        marks=(Mark("alt", 37, 38, "s", "annotated stop codon TGA, annotated_stop_distance = 0"),),
        ruler=Ruler((0, 13, 25, 37, 40)),
    ),
    # ST-14
    Case(
        "delins_ccta_to_g_before_a_second_taa_leaves_a_tga_at_the_stop_codon",
        (
            """
            tx      0             13                25              37  40  43
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TCC TAA taacc..17..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TG- --A taacc..17..ccc] 3'
                                                                     ^^^^^ CCTA>G
                                                                    s annotated stop codon TGA, annotated_stop_distance = 0
            CCTA>G at tx 38 to 41 leaves T, G and the last A of TAA. They read TGA, a stop codon at the place of the
            annotated one, before the second TAA at tx 43.
            """
        ),
        TCC_TAA_TAA,
        Change("GCCGCCT[CCTA>G]ATAACC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CCTA", "TAGG"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(68, 36),
            "variant_end": per_strand(72, 40),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTGA",
            "alt_cds_length": 27,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 15}],
            "alt_last_codon": "TGA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 24,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 24, "codon": "TGA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTGA" + "TAACC" + "C" * 20,
            "alt_transcript_length": 65,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 40}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCGCC[TCCTA>TG]ATAACC"),),
        marks=(Mark("alt", 37, 38, "s", "annotated stop codon TGA, annotated_stop_distance = 0"),),
        ruler=Ruler((0, 13, 25, 37, 40, 43)),
    ),
    # ST-15
    Case(
        "last_sense_codon_and_stop_codon_deleted_before_a_taa_in_the_3utr_is_a_stop_loss",
        (
            """
            tx      0             13                25              37  40  43
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TCC TAA taacccc..15..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC --- --- taacccc..15..ccc] 3'
                                                                    ^^^^^^^ TCCTAA>-
                                                                sss annotated stop codon
                                                                            fff first in-frame stop codon
                                                                <-> annotated_stop_distance = -3
            TCCTAA>- at tx 37 to 42 deletes the last sense codon and the stop codon. The annotated stop codon maps to
            alt tx 34, anchored on the sequence after the deletion. The scan reads on to the TAA at alt tx 37.
            """
        ),
        TCC_TAA_TAA,
        Change("GCCGCC[TCCTAA>]TAACC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CTCCTAA", "ATTAGGA"),
            "alt": per_strand("C", "A"),
            "variant_start": per_strand(66, 34),
            "variant_end": per_strand(73, 41),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCT",
            "alt_cds_length": 25,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 13}],
            "alt_last_codon": "CCT",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCT" + "AACC" + "C" * 20,
            "alt_transcript_length": 62,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAA",
            "alt_scan_first_stop_pos": 37,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 37, "codon": "TAA"}],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -3,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 37}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCGCCT[CCTAAT>]AACC"),),
        marks=(
            Mark("alt", 34, 37, "s", "annotated stop codon"),
            Mark("alt", 37, 40, "f", "first in-frame stop codon"),
            Span("alt", 34, 37, "annotated_stop_distance = -3"),
        ),
        ruler=Ruler((0, 13, 25, 37, 40, 43)),
    ),
    # ST-15
    Case(
        "last_sense_codon_and_stop_codon_deleted_in_a_stop_codon_repeat_is_a_stop_loss",
        (
            """
            tx      0             13                25              37  40  43
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GTA TAG tagcat..17..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC --- --- tagcat..17..ccc] 3'
                                                                    ^^^^^^^ GTATAG>-
                                                                sss annotated stop codon
                                                                            fff first in-frame stop codon
                                                                <-> annotated_stop_distance = -3
            GTATAG>- at tx 37 to 42 deletes the last sense codon and the stop codon. The annotated stop codon maps to
            alt tx 34, anchored on the sequence after the deletion. The scan reads on to the TAG at alt tx 37.
            """
        ),
        GTA_TAG_TAG,
        Change("GCCGCC[GTATAG>]TAGCAT"),
        {
            "variant_id": "var1",
            "ref": per_strand("CGTATAG", "ACTATAC"),
            "alt": per_strand("C", "A"),
            "variant_start": per_strand(66, 35),
            "variant_end": per_strand(73, 42),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCC",
            "alt_cds_length": 24,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 12}],
            "alt_last_codon": "GCC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCC" + "TAGCAT" + "C" * 20,
            "alt_transcript_length": 63,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAG",
            "alt_scan_first_stop_pos": 37,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 37, "codon": "TAG"}],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -3,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 38}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCGC[CGTATAG>C]TAGCAT"),),
        marks=(
            Mark("alt", 34, 37, "s", "annotated stop codon"),
            Mark("alt", 37, 40, "f", "first in-frame stop codon"),
            Span("alt", 34, 37, "annotated_stop_distance = -3"),
        ),
        ruler=Ruler((0, 13, 25, 37, 40, 43)),
    ),
    # ST-15
    Case(
        "deletion_of_the_stop_codon_shifted_left_in_a_stop_codon_repeat_is_a_stop_loss",
        (
            """
            tx      0             13                25            36    40  43
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GTA TAG tagcat..17..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GC- --- --G tagcat..17..ccc] 3'
                                                                  ^^^^^^^^ CGTATA>-
                                                                s annotated stop codon
                                                                            fff first in-frame stop codon
                                                                <---------> annotated_stop_distance = -3
            CGTATA>- at tx 36 to 41 removes the first two bases of the stop codon TAG. The annotated stop codon maps to
            alt tx 34, anchored on the sequence after the deletion. The scan reads on to the TAG at alt tx 37.
            """
        ),
        GTA_TAG_TAG,
        Change("GCCGC[CGTATA>]GTAGCAT"),
        {
            "variant_id": "var1",
            "ref": per_strand("CCGTATA", "CTATACG"),
            "alt": per_strand("C", "C"),
            "variant_start": per_strand(65, 36),
            "variant_end": per_strand(72, 43),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCG",
            "alt_cds_length": 24,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 12}],
            "alt_last_codon": "GCG",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCG" + "TAGCAT" + "C" * 20,
            "alt_transcript_length": 63,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAG",
            "alt_scan_first_stop_pos": 37,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 37, "codon": "TAG"}],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -3,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 38}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCG[CCGTATA>C]GTAGCAT"),),
        marks=(
            Mark("alt", 34, 35, "s", "annotated stop codon"),
            Mark("alt", 37, 40, "f", "first in-frame stop codon"),
            Span("alt", 34, 37, "annotated_stop_distance = -3"),
        ),
        ruler=Ruler((0, 13, 25, 36, 40, 43)),
    ),
    # Pins the threshold 0 of the stop codon classification at annotated_stop_distance 1: "Upstream of the annotated stop
    # codon: a PTC". The deletion of CT removes the start of the stop codon, which "maps it to the position anchored on
    # the unchanged sequence downstream": alt tx 40 - 2 = 38. The first in-frame stop codon TAA at alt tx 37 lies 1 nt
    # upstream of it. The closest PTC case has annotated_stop_distance 3:
    # stop_codon_gained_right_before_the_stop_codon_is_a_ptc.
    Case(
        "deletion_of_ct_in_tac_taa_leaves_a_taa_1_nt_before_the_annotated_stop_codon_and_is_a_ptc",
        """
        tx      0             13                25              37  40
        ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TAC TAA gcccc..12..cccc] 3'
        alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TA- -AA gcccc..12..cccc] 3'
                                                                  ^^^ CT>-
                                                                * PTC TAA
                                                                 s annotated stop codon
                                                                <> annotated_stop_distance = 1
                                                                <---------------------> ptc_to_exon_end = 25, to the transcript end
        The deletion of CT at tx 39 and 40 joins TA of the last sense codon TAC with AA of the stop codon TAA. It
        removes the start of the stop codon, so the annotated stop codon maps to the position anchored on the
        unchanged sequence downstream: alt tx 40 - 2 = 38. The PTC TAA at alt tx 37 lies 1 nt upstream of it.
        The PTC lies in the last exon, so ptc_to_exon_end counts to the transcript end.
        """,
        TAC_TAA_G,
        Change("GCCTA[CT>]AAGCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("ACT", "TAG"),
            "alt": per_strand("A", "T"),
            "variant_start": per_strand(68, 32),
            "variant_end": per_strand(71, 35),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTAAA",
            "alt_cds_length": 28,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 16}],
            "alt_last_codon": "AAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 24,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 24, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTAAA" + "G" + "C" * 20,
            "alt_transcript_length": 62,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 24,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 37,
            "annotated_stop_distance": 1,
            "ptc_to_exon_end": 25,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 37}],
            "nmd_model_status": "ok",
        },
        equivalent=(Change("GCCTA[CTA>A]AGCCC"),),
        marks=(
            Mark("alt", 37, 38, "*", "PTC TAA"),
            Mark("alt", 38, 39, "s", "annotated stop codon"),
            Span("alt", 37, 38, "annotated_stop_distance = 1"),
            Span("alt", 37, 62, "ptc_to_exon_end = 25, to the transcript end"),
        ),
        ruler=Ruler((0, 13, 25, 37, 40)),
    ),
    # Pins the threshold 0 of the stop codon classification at annotated_stop_distance -1: "Downstream of the annotated
    # stop codon ...: a stop loss". Both placements of the deletion, ACTA and CTAA, remove the start of the stop codon,
    # so it maps to the position anchored on the unchanged sequence downstream: alt tx 40 - 4 = 36. The first in-frame
    # stop codon TAG at alt tx 37 lies 1 nt downstream of it. The closest stop loss case,
    # last_sense_codon_and_stop_codon_deleted_before_a_taa_in_the_3utr_is_a_stop_loss, has annotated_stop_distance -3.
    Case(
        "deletion_of_ctaa_in_tac_taa_leaves_a_tag_1_nt_after_the_annotated_stop_codon_and_is_a_stop_loss",
        """
        tx      0             13                25              37  40
        ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TAC TAA gccccc..12..ccc] 3'
        alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC T-- --A gccccc..12..ccc] 3'
                                                                 ^^^^^ ACTA>-
                                                              s annotated stop codon
                                                                f first in-frame stop codon TAG
                                                              <> annotated_stop_distance = -1
        The deletion of ACTA at tx 38 to 41, or of CTAA at tx 39 to 42, joins TA of the last sense codon TAC with
        the g of the 3' UTR. Both placements remove the start of the stop codon, so the annotated stop codon maps to
        the position anchored on the unchanged sequence downstream: alt tx 40 - 4 = 36. The first in-frame stop
        codon TAG at alt tx 37 lies 1 nt downstream of it, a stop loss.
        """,
        TAC_TAA_G,
        Change("GCCT[ACTA>]AGCCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("TACTA", "TTAGT"),
            "alt": per_strand("T", "T"),
            "variant_start": per_strand(67, 31),
            "variant_end": per_strand(72, 36),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTA",
            "alt_cds_length": 26,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 14}],
            "alt_last_codon": "CTA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTA" + "G" + "C" * 20,
            "alt_transcript_length": 60,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAG",
            "alt_scan_first_stop_pos": 37,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 37, "codon": "TAG"}],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -1,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 35}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCTA[CTAA>]GCCCC"), Change("GCCT[ACTAA>A]GCCCC")),
        marks=(
            Mark("alt", 36, 37, "s", "annotated stop codon"),
            Mark("alt", 37, 38, "f", "first in-frame stop codon TAG"),
            Span("alt", 36, 37, "annotated_stop_distance = -1"),
        ),
        ruler=Ruler((0, 13, 25, 37, 40)),
    ),
    # ST-16
    Case(
        "delins_atag_to_t_in_a_stop_codon_repeat_is_a_stop_loss",
        (
            """
            tx      0             13                25              37  40  43
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GTA TAG tagcat..17..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GTT --- tagcat..17..ccc] 3'
                                                                      ^^^^^ ATAG>T
                                                                    sss annotated stop codon
                                                                            fff first in-frame stop codon
                                                                    <-> annotated_stop_distance = -3
            ATAG>T at tx 39 to 42 gives GTT TAG CAT. The annotated stop codon maps to alt tx 37, so the TAG at alt tx 40
            is a stop loss, although the protein is unchanged.
            """
        ),
        GTA_TAG_TAG,
        Change("GCCGT[ATAG>T]TAGCAT"),
        {
            "variant_id": "var1",
            "ref": per_strand("ATAG", "CTAT"),
            "alt": per_strand("T", "A"),
            "variant_start": per_strand(69, 36),
            "variant_end": per_strand(73, 40),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGTT",
            "alt_cds_length": 27,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 15}],
            "alt_last_codon": "GTT",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGTT" + "TAGCAT" + "C" * 20,
            "alt_transcript_length": 66,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAG",
            "alt_scan_first_stop_pos": 40,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 40, "codon": "TAG"}],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -3,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 41}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCG[TATAG>TT]TAGCAT"),),
        marks=(
            Mark("alt", 37, 40, "s", "annotated stop codon"),
            Mark("alt", 40, 43, "f", "first in-frame stop codon"),
            Span("alt", 37, 40, "annotated_stop_distance = -3"),
        ),
        ruler=Ruler((0, 13, 25, 37, 40, 43)),
    ),
    # ST-17
    Case(
        "seven_nt_deletion_of_the_stop_codon_reads_on_to_a_taa_in_the_3utr_as_a_stop_loss",
        (
            """
            tx      0             13                25            36    40  43
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCT TAA ctaagcc..17..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GC- --- --- ctaagcc..17..ccc] 3'
                                                                  ^^^^^^^^^ CGCTTAA>-
                                                              s annotated stop codon
                                                                             fff first in-frame stop codon
                                                              <-------------> annotated_stop_distance = -4
            CGCTTAA>- at tx 36 to 42 deletes the stop codon. The annotated stop codon maps to alt tx 33, anchored on the
            sequence after the deletion. The scan reads on to the TAA at alt tx 37.
            """
        ),
        GCT_TAA_CTAA,
        Change("GCCGC[CGCTTAA>]CTAAGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CCGCTTAA", "GTTAAGCG"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(65, 36),
            "variant_end": per_strand(73, 44),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCC",
            "alt_cds_length": 24,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 12}],
            "alt_last_codon": "GCC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCC" + "TAAGCC" + "C" * 20,
            "alt_transcript_length": 63,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAA",
            "alt_scan_first_stop_pos": 37,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 37, "codon": "TAA"}],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -4,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 38}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCGCC[GCTTAAC>]TAAGCC"),),
        marks=(
            Mark("alt", 33, 34, "s", "annotated stop codon"),
            Mark("alt", 37, 40, "f", "first in-frame stop codon"),
            Span("alt", 33, 37, "annotated_stop_distance = -4"),
        ),
        ruler=Ruler((0, 13, 25, 36, 40, 43)),
    ),
    # ST-18
    Case(
        "stop_codon_snv_behind_an_in_frame_tga_keeps_the_flags_from_the_cds",
        (
            """
            tx      0             13                25  28              40           52
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG TGA GCC GCC GCC TAA cccctgacctagc..15..cccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG TGA GCC GCC GCC CAA cccctgacctagc..15..cccc] 3'
                                                                        ^ T>C
                                                        uuu in-frame TGA of the ref CDS
                                                        *** PTC
                                                                        sss annotated stop codon
                                                        <-------------> annotated_stop_distance = 12
                                                        <-----------------------------------------> ptc_to_exon_end = 47, to the transcript end
            T>C turns the stop codon TAA into CAA. The ref CDS reads the in-frame TGA at tx 28 upstream of it, so the
            row keeps the flags from the CDS.
            """
        ),
        INTERNAL_TGA,
        Change("GCC[T>C]AACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(70, 44),
            "variant_end": per_strand(71, 45),
            "alt_cds_seq": "ATGGCCGCCGCCCTGTGAGCCGCCGCCCAA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "CAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 15,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 15, "codon": "TGA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGTGAGCCGCCGCCCAA" + REF_UTR3,
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TGA",
            "alt_scan_first_stop_pos": 28,
            "alt_scan_stop_codon_count": 2,
            "alt_scan_stop_codons": [{"position": 28, "codon": "TGA"}, {"position": 52, "codon": "TAG"}],
            "alt_scan_stop_codon_exons": [2, 2],
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 15,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 50,
            "annotated_stop_distance": 12,
            "ptc_to_exon_end": 47,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "ref_ptc",
        },
        marks=(
            Mark("ref", 28, 31, "u", "in-frame TGA of the ref CDS"),
            Mark("alt", 28, 31, "*", "PTC"),
            Mark("alt", 40, 43, "s", "annotated stop codon"),
            Span("alt", 28, 40, "annotated_stop_distance = 12"),
            Span("alt", 28, 75, "ptc_to_exon_end = 47, to the transcript end"),
        ),
        ruler=Ruler((0, 13, 25, 28, 40, 52)),
    ),
    # ST-18
    Case(
        "synonymous_snv_behind_an_in_frame_tga_keeps_the_flags_from_the_cds",
        (
            """
            tx      0             13                25  28              40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG TGA GCC GCC GCC TAA ccc..25..cccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG TGA GCC GCC GCT TAA ccc..25..cccc] 3'
                                                                      ^ C>T
                                                        uuu in-frame TGA of the ref CDS
                                                        *** PTC
                                                                        sss annotated stop codon
                                                        <-------------> annotated_stop_distance = 12
                                                        <-------------------------------> ptc_to_exon_end = 47, to the transcript end
            C>T at tx 39 is synonymous: GCC>GCT. The ref CDS reads the in-frame TGA at tx 28 upstream of the annotated
            stop codon, so the row keeps the flags from the CDS.
            """
        ),
        INTERNAL_TGA,
        Change("GCCGC[C>T]TAACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "variant_start": per_strand(69, 45),
            "variant_end": per_strand(70, 46),
            "alt_cds_seq": "ATGGCCGCCGCCCTGTGAGCCGCCGCTTAA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 15,
            "alt_stop_codon_count": 2,
            "alt_stop_codons": [{"position": 15, "codon": "TGA"}, {"position": 27, "codon": "TAA"}],
            "alt_stop_codon_exons": [2, 2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGTGAGCCGCCGCTTAA" + REF_UTR3,
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 15,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 50,
            "annotated_stop_distance": 12,
            "ptc_to_exon_end": 47,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "ref_ptc",
        },
        marks=(
            Mark("ref", 28, 31, "u", "in-frame TGA of the ref CDS"),
            Mark("alt", 28, 31, "*", "PTC"),
            Mark("alt", 40, 43, "s", "annotated stop codon"),
            Span("alt", 28, 40, "annotated_stop_distance = 12"),
            Span("alt", 28, 75, "ptc_to_exon_end = 47, to the transcript end"),
        ),
        ruler=Ruler((0, 13, 25, 28, 40)),
    ),
    # ST-19
    Case(
        "missense_snv_with_stop_codon_rows_on_a_sense_codon_keeps_the_flags_from_the_cds",
        (
            """
            tx      0             13       20       25                  40           52
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TCA cccctgacctagccc..14..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GAC GCC]|[CTG ACC GCC GCC GCC TCA cccctgacctagccc..14..ccc] 3'
                                           ^ C>A
                                                                        xxx stop_codon rows on the sense codon TCA
                                                                                     fff first in-frame stop codon of transcript_seq
            C>A at tx 20 is a missense SNV: GCC>GAC. The stop_codon rows lie on the sense codon TCA, so the row keeps
            the flags from the CDS.
            """
        ),
        TCA_ROWS,
        Change("ATGGCCG[C>A]CGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "variant_start": per_strand(30, 84),
            "variant_end": per_strand(31, 85),
            "alt_cds_seq": "ATGGCCGACGCCCTGACCGCCGCCGCCTCA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "TCA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGACGCCCTGACCGCCGCCGCCTCA" + REF_UTR3,
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        marks=(
            Mark("ref", 40, 43, "x", "stop_codon rows on the sense codon TCA"),
            Mark("ref", 52, 55, "f", "first in-frame stop codon of transcript_seq"),
        ),
        ruler=Ruler((0, 13, 20, 25, 40, 52)),
    ),
    # Pins the threshold of ref_has_ptc, "Whether the first in-frame stop codon starts before the last 3 nt of the
    # CDS", with the first in-frame stop codon 1 nt before them. The stop_codon rows lie on AAC, which "is no stop
    # codon", so the row keeps the flags from the CDS. alt_has_ptc is then "as ref_has_ptc, for the alt CDS":
    # True, a PTC row. The CDS of the closest case has no in-frame stop codon:
    # missense_snv_with_stop_codon_rows_on_a_sense_codon_keeps_the_flags_from_the_cds.
    Case(
        "missense_in_a_cds_whose_in_frame_taa_starts_1_nt_before_its_last_3_nt_keeps_the_ptc_from_the_cds",
        """
        tx      0             13       20       25                  40
        ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA C cc..26..cccc] 3'
        alt 5' [ccgccgccaccgc ATG GCC GAC GCC]|[CTG ACC GCC GCC GCC TAA C cc..26..cccc] 3'
                                       ^ C>A
                                                                    uuu in-frame TAA of the ref CDS
                                                                     x stop_codon rows on AAC, the last 3 nt of the CDS
                                                                    *** PTC
                                                                    <> annotated_stop_distance = 1
                                                                    <----------------> ptc_to_exon_end = 36, to the transcript end
        C>A at tx 20 is a missense SNV: GCC>GAC. The CDS has 31 nt. Its in-frame TAA at tx 40 starts 1 nt before its
        last 3 nt, so ref_has_ptc is True. The stop_codon rows lie on AAC, which is no stop codon, so the row
        keeps the flags from the CDS: alt_has_ptc is True, as ref_has_ptc for the alt CDS.
        The PTC lies in the last exon, so ptc_to_exon_end counts to the transcript end.
        """,
        TAA_C_ROWS,
        Change("ATGGCCG[C>A]CGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "variant_start": per_strand(30, 85),
            "variant_end": per_strand(31, 86),
            "alt_cds_seq": "ATGGCCGACGCCCTGACCGCCGCCGCCTAAC",
            "alt_cds_length": 31,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 19}],
            "alt_last_codon": "AAC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 27,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 27, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGACGCCCTGACCGCCGCCGCCTAAC" + REF_UTR3,
            "alt_transcript_length": 76,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 27,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 51,
            "annotated_stop_distance": 1,
            "ptc_to_exon_end": 36,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "ref_ptc",
        },
        marks=(
            Mark("ref", 40, 43, "u", "in-frame TAA of the ref CDS"),
            Mark("ref", 41, 42, "x", "stop_codon rows on AAC, the last 3 nt of the CDS"),
            Mark("alt", 40, 43, "*", "PTC"),
            Span("alt", 40, 41, "annotated_stop_distance = 1"),
            Span("alt", 40, 76, "ptc_to_exon_end = 36, to the transcript end"),
        ),
        ruler=Ruler((0, 13, 20, 25, 40)),
    ),
    # ST-20
    Case(
        "stop_codon_snv_with_the_annotated_stop_codon_out_of_frame_keeps_the_flags_from_the_cds",
        (
            """
            tx      0              14                 26                40       47
            ref 5' [ccgccgccaccgca TGG CCG CCG CC]|[C TGA CCG CCG CCG CCT AA cccctgac..20..cccc] 3'
            alt 5' [ccgccgccaccgca TGG CCG CCG CC]|[C TGA CCG CCG CCG CCC AA cccctgac..20..cccc] 3'
                                                                        ^ T>C
                                                      uuu in-frame TGA of the ref CDS
                                                      *** PTC
                                                                        s annotated stop codon, out of frame
                                                      <----------------> annotated_stop_distance = 14
                                                      <---------------------------------------> ptc_to_exon_end = 49, to the transcript end
            T>C at tx 40 hits the annotated stop codon TAA, which is out of frame. The CDS starts 1 nt after the ATG,
            and its frame reads the TGA at tx 26 first. So the row keeps the flags from the CDS.
            """
        ),
        OUT_OF_FRAME,
        Change("GCC[T>C]AACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(70, 44),
            "variant_end": per_strand(71, 45),
            "alt_cds_seq": "TGGCCGCCGCCCTGACCGCCGCCGCCCAA",
            "alt_cds_length": 29,
            "alt_cds_exons": [{"exon_number": 1, "length": 11}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "CAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 12,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 12, "codon": "TGA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "A" + "TGGCCGCCGCCCTGACCGCCGCCGCCCAA" + REF_UTR3,
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 14,
            "alt_scan_start_codon_pos": None,
            "alt_scan_start_codon_exon": None,
            "alt_scan_first_stop_codon": "TGA",
            "alt_scan_first_stop_pos": 26,
            "alt_scan_stop_codon_count": 2,
            "alt_scan_stop_codons": [{"position": 26, "codon": "TGA"}, {"position": 47, "codon": "TGA"}],
            "alt_scan_stop_codon_exons": [2, 2],
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": None,
            "ptc_less_than_150nt_to_start": None,
            "ptc_exon_length": 50,
            "annotated_stop_distance": 14,
            "ptc_to_exon_end": 49,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": None,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "ref_ptc",
        },
        marks=(
            Mark("ref", 26, 29, "u", "in-frame TGA of the ref CDS"),
            Mark("alt", 26, 29, "*", "PTC"),
            Mark("alt", 40, 41, "s", "annotated stop codon, out of frame"),
            Span("alt", 26, 40, "annotated_stop_distance = 14"),
            Span("alt", 26, 75, "ptc_to_exon_end = 49, to the transcript end"),
        ),
        ruler=Ruler((0, 14, 26, 40, 47)),
    ),
    # ST-22
    Case(
        "frameshift_whose_first_stop_codon_lies_inside_the_cds_without_stop_codon_rows_is_a_ptc",
        (
            """
            tx      0             13    18          25                  40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC taac..27..cccc] 3'
            alt 5' [ccgccgccaccgc ATG GC- GCC GCC]|[CTG ACC GCC GCC GCC taac..27..cccc] 3'
                                        ^ C>-
                                                     * PTC TGA, in the alt frame
                                                                      e last base of the alt CDS
                                                     <-------------------------------> ptc_to_exon_end = 49, to the transcript end
            C>- at tx 18 shifts the frame: the alt CDS reads ATG GCG CCG CCC TGA. Without stop_codon rows, the CDS ends
            with GCC at tx 39, and a first stop codon inside the alt CDS is a PTC.
            The PTC lies in the last exon, so ptc_to_exon_end counts to the transcript end.
            """
        ),
        NO_STOP_ROWS,
        Change("ATGGC[C>]GCCGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CC", "CG"),
            "alt": per_strand("C", "C"),
            "variant_start": per_strand(27, 85),
            "variant_end": per_strand(29, 87),
            "alt_cds_seq": "ATGGCGCCGCCCTGACCGCCGCCGCC",
            "alt_cds_length": 26,
            "alt_cds_exons": [{"exon_number": 1, "length": 11}, {"exon_number": 2, "length": 15}],
            "alt_last_codon": "GCC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 12,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 12, "codon": "TGA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCGCCGCCCTGACCGCCGCCGCC" + "TAA" + REF_UTR3,
            "alt_transcript_length": 74,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 12,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 50,
            "annotated_stop_distance": None,
            "ptc_to_exon_end": 49,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": [{"exon_number": 1, "length": 24}, {"exon_number": 2, "length": 50}],
            "nmd_model_status": "no_annotated_stop",
        },
        equivalent=(Change("ATGG[C>]CGCCGCC"),),
        marks=(
            Mark("alt", 25, 26, "*", "PTC TGA, in the alt frame"),
            Mark("alt", 38, 39, "e", "last base of the alt CDS"),
            Span("alt", 25, 74, "ptc_to_exon_end = 49, to the transcript end"),
        ),
        ruler=Ruler((0, 13, 18, 25, 40)),
    ),
    # ST-23
    Case(
        "frameshift_whose_first_stop_codon_lies_past_the_cds_end_without_stop_codon_rows_is_neither",
        (
            """
            tx      0             13                25        33        40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC taacccctgacct..19..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GC- GCC GCC taacccctgacct..19..ccc] 3'
                                                              ^ C>-
                                                                      e last base of the alt CDS
                                                                               fff first in-frame stop codon
            C>- at tx 33 shifts the frame. Without stop_codon rows, the CDS ends with GCC at tx 39. The first in-frame
            stop codon lies past the end of the alt CDS, so it is neither a PTC nor a stop loss.
            """
        ),
        NO_STOP_ROWS,
        Change("ACCGC[C>]GCCGCCTAA"),
        {
            "variant_id": "var1",
            "ref": per_strand("CC", "CG"),
            "alt": per_strand("C", "C"),
            "variant_start": per_strand(62, 50),
            "variant_end": per_strand(64, 52),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCGCCGCC",
            "alt_cds_length": 26,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 14}],
            "alt_last_codon": "GCC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCGCCGCC" + "TAA" + REF_UTR3,
            "alt_transcript_length": 74,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 49}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("ACCG[C>]CGCCGCCTAA"),),
        marks=(
            Mark("alt", 38, 39, "e", "last base of the alt CDS"),
            Mark("alt", 46, 49, "f", "first in-frame stop codon"),
        ),
        ruler=Ruler((0, 13, 25, 33, 40)),
    ),
    # ST-24
    Case(
        "frameshift_without_an_in_frame_stop_codon_without_stop_codon_rows_is_neither",
        (
            """
            tx      0             13                25         33        40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GC-C GCC GCC taac..28..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCGC GCC GCC taac..28..ccc] 3'
                                                              ^ ->G
            G inserted after tx 32 shifts the frame. No in-frame stop codon follows up to the transcript end. Without
            stop_codon rows, that is neither a PTC nor a stop loss.
            """
        ),
        NO_STOP_ROWS,
        Change("ACCGC[>G]CGCCGCCTAA"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("CG", "GC"),
            "variant_start": per_strand(62, 51),
            "variant_end": per_strand(63, 52),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCGCGCCGCC",
            "alt_cds_length": 28,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 16}],
            "alt_last_codon": "GCC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCGCGCCGCC" + "TAA" + REF_UTR3,
            "alt_transcript_length": 76,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 25}, {"exon_number": 2, "length": 51}],
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 13, 25, 33, 40)),
    ),
    # ST-25
    Case(
        "stop_gained_in_the_last_codon_without_stop_codon_rows_ends_at_the_cds_end_and_is_a_ptc",
        (
            """
            tx      0             13                25              37  40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TGG taaac..18..cccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TAG taaac..18..cccc] 3'
                                                                     ^ G>A
                                                                    *** PTC
                                                                      e last base of the alt CDS
                                                                    <-----------------> ptc_to_exon_end = 30, to the transcript end
            G>A at tx 38 turns the last sense codon TGG into TAG. Without stop_codon rows, the CDS ends with this codon,
            so the TAG ends at the CDS end and is a PTC.
            The PTC lies in the last exon, so ptc_to_exon_end counts to the transcript end.
            """
        ),
        NO_STOP_ROWS_TGG,
        Change("GCCT[G>A]GTAAACTG"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "variant_start": per_strand(68, 38),
            "variant_end": per_strand(69, 39),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTAG",
            "alt_cds_length": 27,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 15}],
            "alt_last_codon": "TAG",
            "alt_valid_stop": False,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 24,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 24, "codon": "TAG"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTAG" + "TAAACTG" + "C" * 20,
            "alt_transcript_length": 67,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 24,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 42,
            "annotated_stop_distance": None,
            "ptc_to_exon_end": 30,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_annotated_stop",
        },
        marks=(
            Mark("alt", 37, 40, "*", "PTC"),
            Mark("alt", 39, 40, "e", "last base of the alt CDS"),
            Span("alt", 37, 67, "ptc_to_exon_end = 30, to the transcript end"),
        ),
        ruler=Ruler((0, 13, 25, 37, 40)),
    ),
    # ST-25
    Case(
        "missense_in_the_last_codon_without_stop_codon_rows_stops_right_after_the_cds_end_and_is_neither",
        (
            """
            tx      0             13                25              37  40
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TGG taaact..18..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC TGC taaact..18..ccc] 3'
                                                                      ^ G>C
                                                                      e last base of the alt CDS
                                                                        fff first in-frame stop codon
            G>C at tx 39 turns the last sense codon TGG into TGC. Without stop_codon rows, the CDS ends with this codon.
            The first in-frame stop codon is the TAA right after the CDS end, so it is neither a PTC nor a stop loss.
            """
        ),
        NO_STOP_ROWS_TGG,
        Change("GCCTG[G>C]TAAACTG"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(69, 37),
            "variant_end": per_strand(70, 38),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTGC",
            "alt_cds_length": 27,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 15}],
            "alt_last_codon": "TGC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCTGC" + "TAAACTG" + "C" * 20,
            "alt_transcript_length": 67,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        marks=(
            Mark("alt", 39, 40, "e", "last base of the alt CDS"),
            Mark("alt", 40, 43, "f", "first in-frame stop codon"),
        ),
        ruler=Ruler((0, 13, 25, 37, 40)),
    ),
    # ST-26
    Case(
        "stop_loss_whose_alt_transcript_ends_in_a_stop_codon_out_of_frame",
        (
            """
            tx      0             13                25                  40           52
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA cccctgacctagccc..11..ccctga] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC CAA cccctgacctagccc..11..ccctga] 3'
                                                                        ^ T>C
                                                                        sss annotated stop codon
                                                                                     fff first in-frame stop codon
                                                                        <-----------> annotated_stop_distance = -12
                                                                                                    lll last 3 nt TGA, out of frame
            T>C turns the stop codon TAA into CAA. The scan reads on to the in-frame TAG at alt tx 52. The transcript
            ends in TGA, which is out of frame.
            """
        ),
        TGA_END,
        Change("GCC[T>C]AACCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(70, 44),
            "variant_end": per_strand(71, 45),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "CAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA" + "CCCCTGACCTAG" + "C" * 17 + "TGA",
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAG",
            "alt_scan_first_stop_pos": 52,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 52, "codon": "TAG"}],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -12,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        marks=(
            Mark("alt", 40, 43, "s", "annotated stop codon"),
            Mark("alt", 52, 55, "f", "first in-frame stop codon"),
            Span("alt", 40, 52, "annotated_stop_distance = -12"),
            Mark("alt", 72, 75, "l", "last 3 nt TGA, out of frame"),
        ),
        ruler=Ruler((0, 13, 25, 40, 52)),
    ),
    # ST-27
    Case(
        "stop_codon_snv_at_the_transcript_end_is_a_nonstop",
        (
            """
            tx      0             13                25                  40 43
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC CAA] 3'
                                                                        ^ T>C
                                                                        sss annotated stop codon
            T>C turns the stop codon TAA into CAA. TAA is the last codon of the transcript, so no in-frame stop codon
            follows up to the transcript end.
            """
        ),
        NO_UTR3,
        Change("GCCGCC[T>C]AA"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(70, 12),
            "variant_end": per_strand(71, 13),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "CAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA",
            "alt_transcript_length": 43,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": None,
            "alt_scan_first_stop_pos": None,
            "alt_scan_stop_codon_count": 0,
            "alt_scan_stop_codons": [],
            "alt_scan_stop_codon_exons": [],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        marks=(Mark("alt", 40, 43, "s", "annotated stop codon"),),
        ruler=Ruler((0, 13, 25, 40, 43)),
    ),
    # SP-01
    Case(
        "snv_in_the_part_of_a_split_start_codon_in_the_second_exon_is_a_start_loss",
        (
            """
            tx      0             13   15                   40
            ref 5' [ccgccgccaccgc AT]|[G GCC GCC ..15.. GCC TAA c..28..ccc] 3'
            alt 5' [ccgccgccaccgc AT]|[C GCC GCC ..15.. GCC TAA c..28..ccc] 3'
                                       ^ G>C
            G>C at tx 15 turns the start codon ATG into ATC: a start loss. The start codon is split AT|G over the exon
            junction, and the SNV hits its part in exon 2.
            """
        ),
        SPLIT_START,
        Change("TTTCAG[G>C]GCCGCCGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(45, 69),
            "variant_end": per_strand(46, 70),
            "alt_cds_seq": "ATCGCCGCCGCCCTGACCGCCGCCGCCTAA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 2}, {"exon_number": 2, "length": 28}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 27,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 27, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATCGCCGCCGCCCTGACCGCCGCCGCCTAA" + REF_UTR3,
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": None,
            "alt_scan_start_codon_exon": None,
            "alt_scan_first_stop_codon": None,
            "alt_scan_first_stop_pos": None,
            "alt_scan_stop_codon_count": 0,
            "alt_scan_stop_codons": [],
            "alt_scan_stop_codon_exons": [],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 13, 15, 40)),
    ),
    # SP-02
    Case(
        "ensembl_start_codon_split_over_an_intron_is_an_annotated_start_codon",
        (
            """
            tx      0             13   15     20                         40
            ref 5' [ccgccgccaccgc AT]|[G GCC GCC GCC CTG ACC GCC GCC GCC TAA c..28..ccc] 3'
            alt 5' [ccgccgccaccgc AT]|[G GCC GAC GCC CTG ACC GCC GCC GCC TAA c..28..ccc] 3'
                                              ^ C>A
            C>A at tx 20 is a missense SNV: GCC>GAC. In the Ensembl flavor, the start codon AT|G is split over the
            intron, and it counts as the annotated start codon.
            """
        ),
        SPLIT_START_ENSEMBL,
        Change("TTTCAGGGCCG[C>A]CGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "variant_start": per_strand(50, 64),
            "variant_end": per_strand(51, 65),
            "alt_cds_seq": "ATGGCCGACGCCCTGACCGCCGCCGCCTAA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 2}, {"exon_number": 2, "length": 28}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 27,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 27, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGACGCCCTGACCGCCGCCGCCTAA" + REF_UTR3,
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 13, 15, 20, 40)),
    ),
    # SP-03, SP-04
    Case(
        "split_stop_codon_taa_to_tag_in_its_one_base_cds_row_is_neither_ptc_nor_stop_loss",
        (
            """
            tx      0             13                25                  40   42         52
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TA]|[A cccctgacctagc..16..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TA]|[G cccctgacctagc..16..ccc] 3'
                                                                             ^ A>G
            A>G at tx 42 turns the split stop codon TA|A into TAG. The A is the only CDS base of exon 3, and the stop
            codon stays a stop codon.
            """
        ),
        SPLIT_STOP,
        Change("TTTCAG[A>G]CCCCTGA"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(92, 42),
            "variant_end": per_strand(93, 43),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTAG",
            "alt_cds_length": 30,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 12},
                {"exon_number": 2, "length": 17},
                {"exon_number": 3, "length": 1},
            ],
            "alt_last_codon": "TAG",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 27,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 27, "codon": "TAG"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCTAG" + REF_UTR3,
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 13, 25, 40, 42, 52)),
    ),
    # SP-03, SP-04
    Case(
        "split_stop_codon_taa_to_caa_reads_through_into_the_next_exon_as_a_stop_loss",
        (
            """
            tx      0             13                25                  40   42         52
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TA]|[A cccctgacctagccc..14..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC CA]|[A cccctgacctagccc..14..ccc] 3'
                                                                        ^ T>C
                                                                        s annotated stop codon
                                                                                        fff first in-frame stop codon, in exon 3
                                                                        <--------------> annotated_stop_distance = -12
            T>C at tx 40 turns the split stop codon TA|A into CAA. The scan reads on into exon 3, to the in-frame TAG at
            alt tx 52.
            """
        ),
        SPLIT_STOP,
        Change("GCCGCCGCC[T>C]AGTAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(70, 64),
            "variant_end": per_strand(71, 65),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA",
            "alt_cds_length": 30,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 12},
                {"exon_number": 2, "length": 17},
                {"exon_number": 3, "length": 1},
            ],
            "alt_last_codon": "CAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA" + REF_UTR3,
            "alt_transcript_length": 75,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAG",
            "alt_scan_first_stop_pos": 52,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 52, "codon": "TAG"}],
            "alt_scan_stop_codon_exons": [3],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -12,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        marks=(
            Mark("alt", 40, 41, "s", "annotated stop codon"),
            Mark("alt", 52, 55, "f", "first in-frame stop codon, in exon 3"),
            Span("alt", 40, 52, "annotated_stop_distance = -12"),
        ),
        ruler=Ruler((0, 13, 25, 40, 42, 52)),
    ),
    # SP-05
    Case(
        "ptc_split_with_one_base_before_the_last_exon_junction_lies_in_the_exon_of_its_first_base",
        (
            """
            tx      0             13              25  26         34
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC C]|[AG GCC GCC TAA cccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC T]|[AG GCC GCC TAA cccc] 3'
                                                  ^ C>T
                                                  * PTC TAG
                                                                 sss annotated stop codon
                                                  <> ptc_to_exon_end = 1
            C>T at tx 25 turns CAG into TAG. This PTC is split T|AG over the last exon junction, and it lies in exon 1,
            the exon of its first base.
            """
        ),
        SPLIT_CODON_1,
        Change("GCCGCC[C>T]GTAAGT"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "variant_start": per_strand(35, 45),
            "variant_end": per_strand(36, 46),
            "alt_cds_seq": "ATGGCCGCCGCCTAGGCCGCCTAA",
            "alt_cds_length": 24,
            "alt_cds_exons": [{"exon_number": 1, "length": 13}, {"exon_number": 2, "length": 11}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 12,
            "alt_stop_codon_count": 2,
            "alt_stop_codons": [{"position": 12, "codon": "TAG"}, {"position": 21, "codon": "TAA"}],
            "alt_stop_codon_exons": [1, 2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCTAGGCCGCCTAA" + "CCCC",
            "alt_transcript_length": 41,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 0,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 12,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 26,
            "annotated_stop_distance": 9,
            "ptc_to_exon_end": 1,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": True,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Mark("alt", 25, 26, "*", "PTC TAG"),
            Mark("alt", 34, 37, "s", "annotated stop codon"),
            Span("alt", 25, 26, "ptc_to_exon_end = 1"),
        ),
        ruler=Ruler((0, 13, 25, 26, 34)),
    ),
    # SP-05, SP-06
    Case(
        "snv_at_an_exon_start_completes_a_ptc_with_two_bases_before_the_last_exon_junction",
        (
            """
            tx      0             13              25   27        34
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC TA]|[C GCC GCC TAA cccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC TA]|[G GCC GCC TAA cccc] 3'
                                                       ^ C>G
                                                  * PTC TAG
                                                                 sss annotated stop codon
                                                  <> ptc_to_exon_end = 2
            C>G at tx 27, the first base of exon 2, turns TAC into TAG. This PTC is split TA|G over the last exon
            junction, and it lies in exon 1, the exon of its first base.
            """
        ),
        SPLIT_CODON_2,
        Change("TTTCAG[C>G]GCCGCCTAA"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(57, 23),
            "variant_end": per_strand(58, 24),
            "alt_cds_seq": "ATGGCCGCCGCCTAGGCCGCCTAA",
            "alt_cds_length": 24,
            "alt_cds_exons": [{"exon_number": 1, "length": 14}, {"exon_number": 2, "length": 10}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 12,
            "alt_stop_codon_count": 2,
            "alt_stop_codons": [{"position": 12, "codon": "TAG"}, {"position": 21, "codon": "TAA"}],
            "alt_stop_codon_exons": [1, 2],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCTAGGCCGCCTAA" + "CCCC",
            "alt_transcript_length": 41,
            "alt_cds_start_in_transcript": 13,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 0,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 12,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 27,
            "annotated_stop_distance": 9,
            "ptc_to_exon_end": 2,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": True,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Mark("alt", 25, 26, "*", "PTC TAG"),
            Mark("alt", 34, 37, "s", "annotated stop codon"),
            Span("alt", 25, 27, "ptc_to_exon_end = 2"),
        ),
        ruler=Ruler((0, 13, 25, 27, 34)),
    ),
    # SP-07
    Case(
        "readthrough_stop_codon_split_over_an_exon_junction_lies_in_the_exon_of_its_first_base",
        (
            """
            tx      0             13                25                  40     46   48
            ref 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA cccta]|[gccc..14..ccc] 3'
            alt 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC CAA cccta]|[gccc..14..ccc] 3'
                                                                        ^ T>C
                                                                        sss annotated stop codon
                                                                               f first in-frame stop codon TA|G
                                                                        <-----> annotated_stop_distance = -6
            T>C turns the stop codon TAA into CAA. The scan reads on to the in-frame TAG at alt tx 46. It is split TA|G
            over the exon junction and lies in exon 2, the exon of its first base.
            """
        ),
        SPLIT_READTHROUGH,
        Change("GCC[T>C]AACCCTA"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(70, 58),
            "variant_end": per_strand(71, 59),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA",
            "alt_cds_length": 30,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}, {"exon_number": 2, "length": 18}],
            "alt_last_codon": "CAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": REF_UTR5 + "ATGGCCGCCGCCCTGACCGCCGCCGCCCAA" + "CCCTA" + "G" + "C" * 20,
            "alt_transcript_length": 69,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAG",
            "alt_scan_first_stop_pos": 46,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 46, "codon": "TAG"}],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -6,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        marks=(
            Mark("alt", 40, 43, "s", "annotated stop codon"),
            Mark("alt", 46, 47, "f", "first in-frame stop codon TA|G"),
            Span("alt", 40, 46, "annotated_stop_distance = -6"),
        ),
        ruler=Ruler((0, 13, 25, 40, 46, 48)),
    ),
    # ST-28. The deletion of the T of the start codon shortens alt exon 1 to 10 nt. The scan finds the ATG at tx 7,
    # and its PTC TAA at tx 10 is the first base of alt exon 2, the last exon. The exon numbers of the scan come from
    # the alt transcript: alt_scan_stop_codon_exons is [2], and the PTC features and rules are those of a PTC in
    # the last exon.
    Case(
        "ptc_of_the_scan_at_the_first_base_of_the_last_exon_after_a_start_codon_deletion",
        """
        alt tx     0    4    7      10           20
        ref    5' [gacc ATG GAT G]|[TA AGC TAA gc] 3'
        alt    5' [gacc A-G GAT G]|[TA AGC TAA gc] 3'
                         ^ T>-
                             a ATG of the scan
                                    * PTC TAA
                                    <-----------> ptc_to_exon_end = 10, to the transcript end
                             <--> ptc_to_start_codon = 3
        T>- deletes the T of the start codon: ATG>AG, a start loss. It shortens alt exon 1 to 10 nt.
        The scan finds the ATG at alt tx 7, in exon 1. Its PTC TAA at alt tx 10 is the first base of alt exon 2, the
        last exon.
        """,
        TWO_EXONS,
        Change("GACCA[T>]GGATG"),
        {
            "variant_id": "var1",
            "ref": per_strand("AT", "CA"),
            "alt": per_strand("A", "C"),
            "variant_start": per_strand(14, 44),
            "variant_end": per_strand(16, 46),
            "alt_cds_seq": "AGGATGTAAGCTAA",
            "alt_cds_length": 14,
            "alt_cds_exons": [{"exon_number": 1, "length": 6}, {"exon_number": 2, "length": 8}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 6,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 6, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": True,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GACCAGGATGTAAGCTAAGC",
            "alt_transcript_length": 20,
            "alt_cds_start_in_transcript": 4,
            "alt_scan_start_codon_pos": 7,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAA",
            "alt_scan_first_stop_pos": 10,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 10, "codon": "TAA"}],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 3,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 10,
            "annotated_stop_distance": 5,
            "ptc_to_exon_end": 10,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": [{"exon_number": 1, "length": 10}, {"exon_number": 2, "length": 10}],
            "nmd_model_status": "ok",
        },
        equivalent=(Change("GACC[ATG>AG]GATG"),),
        marks=(
            Mark("alt", 7, 8, "a", "ATG of the scan"),
            Mark("alt", 10, 11, "*", "PTC TAA"),
            Span("alt", 10, 20, "ptc_to_exon_end = 10, to the transcript end"),
            Span("alt", 7, 10, "ptc_to_start_codon = 3"),
        ),
        ruler=Ruler((0, 4, 7, 10, 20), "tx", "alt"),
    ),
]
