"""
Conformance cases of the NMD features: the values of the columns that add_nmd_features adds (utr5_length,
utr3_length, total_exon_count, upstream_exon_count, downstream_exon_count, ptc_to_start_codon,
ptc_less_than_150nt_to_start, ptc_exon_length, stop_codon_distance, ptc_to_intron and likely_misannotated).
The null cases of these columns are in cases_null.py.

The drawings show the features as the figures of "Technical Notes.md" do. They show the transcript 5' to 3', also
for the minus strand. The runner renders their layout block from the case. The line ref is the transcript and the
line alt is the transcript with the change. The two lines are aligned, and `-` fills a gap. `[...]` is an exon and
`|` is an exon junction. Upper case is the coding region, in codons of the annotated frame; lower case is UTR. An
intron that the change touches is drawn in lower case outside the brackets. `..N..` leaves out N bases. `^` marks the
change as ref>alt. A mark line puts `***` under the PTC, and `<-- label -->` spans a length, e.g. from the PTC to an
exon junction or to the transcript end. A ruler gives tx, CDS or alt CDS positions, as labelled. CDS positions count
from the first coding base, so the 5' UTR is negative.
"""

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
    Ruler,
    Span,
    Transcript,
    per_strand,
)

# PTC in exon 2 of 3: the figure "Features of a PTC in exon 2 of 3" of "Technical Notes.md"
FIGURE_CDS = "ATG" + "GCC" * 59 + "CAG" + "GCC" * 58 + "TAA"
FIGURE = Layout(
    Transcript(("c" * 50 + FIGURE_CDS[:100], FIGURE_CDS[100:300], FIGURE_CDS[300:] + "c" * 190)),
    {
        **IDS,
        "ref_cds_start": per_strand(60, 200),
        "ref_cds_stop": per_strand(460, 600),
        "ref_cds_seq": FIGURE_CDS,
        "ref_cds_len": 360,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 100), (2, 200), (3, 60)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 357,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(357, "TAA")],
        "ref_stop_codon_exons": [3],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 650,
        "transcript_seq": "C" * 50 + FIGURE_CDS + "C" * 190,
        "transcript_length": 600,
        "cds_start_in_transcript": 50,
        "cds_end_in_transcript": 410,
        "transcript_exon_info": [(1, 150), (2, 200), (3, 250)],
        "utr3_length": 190,
        "utr5_length": 50,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# The whole CDS in exon 2 of 3: exon 1 holds 5' UTR only, exon 3 holds 3' UTR only
INSIDE_CDS = "ATG" + "GCC" * 4 + "CAG" + "GCC" * 3 + "TAA"
CDS_INSIDE_EXON_2 = Layout(
    Transcript(("c" * 20, "c" * 15 + INSIDE_CDS + "c" * 25, "c" * 30)),
    {
        **IDS,
        "ref_cds_start": per_strand(65, 85),
        "ref_cds_stop": per_strand(95, 115),
        "ref_cds_seq": INSIDE_CDS,
        "ref_cds_len": 30,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(2, 30)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 2,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 27,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(27, "TAA")],
        "ref_stop_codon_exons": [2],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 170,
        "transcript_seq": "C" * 35 + INSIDE_CDS + "C" * 55,
        "transcript_length": 120,
        "cds_start_in_transcript": 35,
        "cds_end_in_transcript": 65,
        "transcript_exon_info": [(1, 20), (2, 70), (3, 30)],
        "utr3_length": 55,
        "utr5_length": 35,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# Three exons with a long last exon. The CDS repeats GCC, whose 3 frames hold no stop codon. CAG at CDS 30 can become
# the PTC TAG. CTA ACA at CDS 36 holds a TAA at 37, in the frame after a 1 nt deletion upstream. CCT AAC at CDS 45
# holds a TAA at 47, in the frame after a 1 nt insertion upstream. The 3' UTR holds a TAG in the frame of the CDS.
LAST_EXON_CDS = "ATGGCCGCCGCC" + "GCCAAGGCCGCCGCC" + "GCCCAGGCCCTAACAGCCCCTAACGCCTAA"
LAST_EXON_UTR3 = "gccgcctag" + "c" * 16
THREE_EXONS = Layout(
    Transcript(("c" * 8 + LAST_EXON_CDS[:12], LAST_EXON_CDS[12:27], LAST_EXON_CDS[27:] + LAST_EXON_UTR3)),
    {
        **IDS,
        "ref_cds_start": per_strand(18, 35),
        "ref_cds_stop": per_strand(115, 132),
        "ref_cds_seq": LAST_EXON_CDS,
        "ref_cds_len": 57,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 12), (2, 15), (3, 30)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 54,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(54, "TAA")],
        "ref_stop_codon_exons": [3],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 140,
        "transcript_seq": "C" * 8 + LAST_EXON_CDS + "GCCGCCTAG" + "C" * 16,
        "transcript_length": 90,
        "cds_start_in_transcript": 8,
        "cds_end_in_transcript": 65,
        "transcript_exon_info": [(1, 20), (2, 15), (3, 55)],
        "utr3_length": 25,
        "utr5_length": 8,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# The transcript of THREE_EXONS as a single exon
SINGLE_EXON = Layout(
    Transcript(("c" * 8 + LAST_EXON_CDS + LAST_EXON_UTR3,)),
    {
        **IDS,
        "ref_cds_start": per_strand(18, 35),
        "ref_cds_stop": per_strand(75, 92),
        "ref_cds_seq": LAST_EXON_CDS,
        "ref_cds_len": 57,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 57)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 54,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(54, "TAA")],
        "ref_stop_codon_exons": [1],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 100,
        "transcript_seq": "C" * 8 + LAST_EXON_CDS + "GCCGCCTAG" + "C" * 16,
        "transcript_length": 90,
        "cds_start_in_transcript": 8,
        "cds_end_in_transcript": 65,
        "transcript_exon_info": [(1, 90)],
        "utr3_length": 25,
        "utr5_length": 8,
        "total_exon_count": 1,
        "likely_misannotated": False,
    },
)

# Three exons. The coding exon 2 has 13 nt, so deleting it shifts the frame. In the frame after that deletion, exon 3
# starts with the stop codon TGA.
EXON_2_OF_13_NT = Layout(
    Transcript(("gccacATGGCCAAGCTGCTG", "CAGCAGCTGCTGC", "TGAAGCTGTAAgcccccccc")),
    {
        **IDS,
        "ref_cds_start": per_strand(15, 19),
        "ref_cds_stop": per_strand(94, 98),
        "ref_cds_seq": "ATGGCCAAGCTGCTG" + "CAGCAGCTGCTGC" + "TGAAGCTGTAA",
        "ref_cds_len": 39,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 15), (2, 13), (3, 11)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 36,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(36, "TAA")],
        "ref_stop_codon_exons": [3],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 103,
        "transcript_seq": "GCCACATGGCCAAGCTGCTG" + "CAGCAGCTGCTGC" + "TGAAGCTGTAAGCCCCCCCC",
        "transcript_length": 53,
        "cds_start_in_transcript": 5,
        "cds_end_in_transcript": 44,
        "transcript_exon_info": [(1, 20), (2, 13), (3, 20)],
        "utr3_length": 9,
        "utr5_length": 5,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# Exon 1 holds 1 nt of 5' UTR. TGG at CDS 6 can become the PTC TGA.
UTR_EXON_1_OF_1_NT = Layout(
    Transcript(("g", "accATGGCCTGGAAGTAAgcc")),
    {
        **IDS,
        "ref_cds_start": per_strand(34, 13),
        "ref_cds_stop": per_strand(49, 28),
        "ref_cds_seq": "ATGGCCTGGAAGTAA",
        "ref_cds_len": 15,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(2, 15)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 2,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 12,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(12, "TAA")],
        "ref_stop_codon_exons": [2],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 52,
        "transcript_seq": "G" + "ACCATGGCCTGGAAGTAAGCC",
        "transcript_length": 22,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 19,
        "transcript_exon_info": [(1, 1), (2, 21)],
        "utr3_length": 3,
        "utr5_length": 4,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)

# The CDS starts at the transcript start: no 5' UTR
NO_UTR5 = Layout(
    Transcript(("ATGGCCAAGCTG", "GGCTCCTAAgccgcc")),
    {
        **IDS,
        "ref_cds_start": per_strand(10, 16),
        "ref_cds_stop": per_strand(51, 57),
        "ref_cds_seq": "ATGGCCAAGCTGGGCTCCTAA",
        "ref_cds_len": 21,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 12), (2, 9)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 18,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(18, "TAA")],
        "ref_stop_codon_exons": [2],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 57,
        "transcript_seq": "ATGGCCAAGCTGGGCTCCTAAGCCGCC",
        "transcript_length": 27,
        "cds_start_in_transcript": 0,
        "cds_end_in_transcript": 21,
        "transcript_exon_info": [(1, 12), (2, 15)],
        "utr3_length": 6,
        "utr5_length": 0,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)

# The stop codon ends the transcript: no 3' UTR
NO_UTR3 = Layout(
    Transcript(("gaccATGGCCAAGCTG", "GGCTCCTAA")),
    {
        **IDS,
        "ref_cds_start": per_strand(14, 10),
        "ref_cds_stop": per_strand(55, 51),
        "ref_cds_seq": "ATGGCCAAGCTGGGCTCCTAA",
        "ref_cds_len": 21,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 12), (2, 9)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 18,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(18, "TAA")],
        "ref_stop_codon_exons": [2],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 55,
        "transcript_seq": "GACCATGGCCAAGCTGGGCTCCTAA",
        "transcript_length": 25,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 25,
        "transcript_exon_info": [(1, 16), (2, 9)],
        "utr3_length": 0,
        "utr5_length": 4,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)

# ATG GAT GCC CTA: the ATG at CDS 4 starts the ORF ATG CCC TAA, out of frame with the annotated stop codon
RESCUE_CDS = "ATGGATGCCCTA" + "AAA" * 50 + "GAC" + "TAA"
RESCUE = Layout(
    Transcript(("ggg" + RESCUE_CDS[:60], RESCUE_CDS[60:120], RESCUE_CDS[120:] + "ggggg")),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 15),
        "ref_cds_stop": per_strand(221, 223),
        "ref_cds_seq": RESCUE_CDS,
        "ref_cds_len": 168,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 60), (2, 60), (3, 48)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 165,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(165, "TAA")],
        "ref_stop_codon_exons": [3],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 226,
        "transcript_seq": "GGG" + RESCUE_CDS + "GGGGG",
        "transcript_length": 176,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 171,
        "transcript_exon_info": [(1, 63), (2, 60), (3, 53)],
        "utr3_length": 5,
        "utr5_length": 3,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# CAG at CDS 150 can become the PTC TAG, 150 nt from the start codon
PTC_150_CDS = "ATG" + "GCC" * 49 + "CAG" + "GCC" * 20 + "TAA"
PTC_150 = Layout(
    Transcript(("gacc" + PTC_150_CDS[:210], PTC_150_CDS[210:] + "gccgcc")),
    {
        **IDS,
        "ref_cds_start": per_strand(14, 16),
        "ref_cds_stop": per_strand(250, 252),
        "ref_cds_seq": PTC_150_CDS,
        "ref_cds_len": 216,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 210), (2, 6)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 213,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(213, "TAA")],
        "ref_stop_codon_exons": [2],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 256,
        "transcript_seq": "GACC" + PTC_150_CDS + "GCCGCC",
        "transcript_length": 226,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 220,
        "transcript_exon_info": [(1, 214), (2, 12)],
        "utr3_length": 6,
        "utr5_length": 4,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)


# The annotated start codon is CTG. The ATG at CDS 30 is an internal Met. TGG at CDS 159 can become the PTC TAG.
CTG_CDS = "CTG" + "AAA" * 9 + "ATG" + "AAA" * 42 + "TGG" + "GAC" + "TAA"
CTG_START = Layout(
    Transcript(("ggg" + CTG_CDS[:100], CTG_CDS[100:] + "ggggg")),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 15),
        "ref_cds_stop": per_strand(201, 203),
        "ref_cds_seq": CTG_CDS,
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
        "transcript_seq": "GGG" + CTG_CDS + "GGGGG",
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

# The transcript of CTG_START without start_codon rows, tagged cds_start_NF
CDS_START_NF = Layout(
    Transcript(("ggg" + CTG_CDS[:100], CTG_CDS[100:] + "ggggg"), start_codon=False, tags=("cds_start_NF",)),
    {
        **CTG_START.ref,
        "has_start_codon": False,
        "ref_start_codon_pos": None,
        "ref_start_codon_exon": None,
        "likely_misannotated": True,
    },
)

# No stop_codon rows, tagged cds_end_NF: the CDS ends in the sense codon TCC at the transcript end
CDS_END_NF = Layout(
    Transcript(("gaccATGGCCAAGCTG", "GGCTCC"), stop_codon=False, tags=("cds_end_NF",)),
    {
        **IDS,
        "ref_cds_start": per_strand(14, 10),
        "ref_cds_stop": per_strand(52, 48),
        "ref_cds_seq": "ATGGCCAAGCTGGGCTCC",
        "ref_cds_len": 18,
        "has_start_codon": True,
        "has_stop_codon": False,
        "cds_frame": 0,
        "ref_cds_info": [(1, 12), (2, 6)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TCC",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 52,
        "transcript_seq": "GACCATGGCCAAGCTGGGCTCC",
        "transcript_length": 22,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 22,
        "transcript_exon_info": [(1, 16), (2, 6)],
        "utr3_length": None,
        "utr5_length": 4,
        "total_exon_count": 2,
        "likely_misannotated": True,
    },
)

# The stop_codon rows lie on the sense codon TCC
STOP_CODON_ROWS_ON_A_SENSE_CODON = Layout(
    Transcript(("gaccATGGCCAAGCTG", "GGCTCCgccgcc")),
    {
        **IDS,
        "ref_cds_start": per_strand(14, 16),
        "ref_cds_stop": per_strand(52, 54),
        "ref_cds_seq": "ATGGCCAAGCTGGGCTCC",
        "ref_cds_len": 18,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 12), (2, 6)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TCC",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 58,
        "transcript_seq": "GACCATGGCCAAGCTGGGCTCCGCCGCC",
        "transcript_length": 28,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 22,
        "transcript_exon_info": [(1, 16), (2, 12)],
        "utr3_length": 6,
        "utr5_length": 4,
        "total_exon_count": 2,
        "likely_misannotated": True,
    },
)

# Ensembl GFF3 on chrM: the CDS ends in AGA, a mitochondrial stop codon that the codon scans do not know
MITOCHONDRIAL_AGA = Layout(
    Transcript(("ATGGCCAAGCTGGGCAGA",), contig="chrM", flavor="ensembl"),
    {
        **IDS,
        "chromosome": "chrM",
        "ref_cds_start": per_strand(10, 10),
        "ref_cds_stop": per_strand(28, 28),
        "ref_cds_seq": "ATGGCCAAGCTGGGCAGA",
        "ref_cds_len": 18,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 18)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "AGA",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 28,
        "transcript_seq": "ATGGCCAAGCTGGGCAGA",
        "transcript_length": 18,
        "cds_start_in_transcript": 0,
        "cds_end_in_transcript": 18,
        "transcript_exon_info": [(1, 18)],
        "utr3_length": 0,
        "utr5_length": 0,
        "total_exon_count": 1,
        "likely_misannotated": True,
    },
)


CASES = [
    # PF-07, PF-09, PF-10, PF-15 (positive), PF-16, PF-27 (internal exon)
    Case(
        "ptc_in_exon_2_of_3_has_the_features_of_the_figure",
        """
        C>T at CDS 180 makes CAG the PTC TAG. The values are those of the figure "Features of a PTC in exon 2 of 3".
        The PTC exon 2 is the smallest exon number in alt_stop_codon_exons [2, 3].
        exon 1: 150 nt, exon 2: 200 nt, exon 3: 250 nt
        upstream_exon_count = 1, downstream_exon_count = 1. No NMD rule is True.

        CDS     -50            0                      100                   180                           300                    357                550
        ref 5' [cccc..42..cccc ATG GCC ..90.. GCC G]|[CC GCC ..69.. GCC GCC CAG GCC GCC ..105.. GCC GCC]|[GCC GCC ..45.. GCC GCC TAA cccc..182..cccc] 3'
        alt 5' [cccc..42..cccc ATG GCC ..90.. GCC G]|[CC GCC ..69.. GCC GCC TAG GCC GCC ..105.. GCC GCC]|[GCC GCC ..45.. GCC GCC TAA cccc..182..cccc] 3'
                                                                            ^ C>T
                                                                            *** PTC
                <------------> utr5_length = 50
                               <-- ptc_to_start_codon = 180, not < 150 --->
                                                      <------------ ptc_exon_length = 200 ------------>
                                                                            <-- ptc_to_intron = 120 -->
                                                                            <----------- stop_codon_distance = 177 ------------>
                                                                                                                                     <-------------> utr3_length = 190
        """,
        FIGURE,
        Change("GCC[C>T]AGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "start_variant": per_strand(260, 399),
            "end_variant": per_strand(261, 400),
            "alt_cds_start": per_strand(60, 200),
            "alt_cds_stop": per_strand(460, 600),
            "alt_cds_seq": "ATG" + "GCC" * 59 + "TAG" + "GCC" * 58 + "TAA",
            "alt_cds_len": 360,
            "alt_cds_info": [(1, 100), (2, 200), (3, 60)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 180,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(180, "TAG"), (357, "TAA")],
            "alt_stop_codon_exons": [2, 3],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "C" * 50 + "ATG" + "GCC" * 59 + "TAG" + "GCC" * 58 + "TAA" + "C" * 190,
            "alt_transcript_length": 600,
            "alt_cds_start_in_transcript": 50,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 180,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 200,
            "stop_codon_distance": 177,
            "ptc_to_intron": 120,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Mark("alt", 230, 233, "*", "PTC"),
            Span("alt", 0, 50, "utr5_length = 50"),
            Span("alt", 50, 230, "ptc_to_start_codon = 180, not < 150"),
            Span("alt", 150, 350, "ptc_exon_length = 200"),
            Span("alt", 230, 350, "ptc_to_intron = 120"),
            Span("alt", 230, 407, "stop_codon_distance = 177"),
            Span("alt", 410, 600, "utr3_length = 190"),
        ),
        ruler=Ruler((-50, 0, 100, 180, 300, 357, 550), "CDS"),
    ),
    # PF-01, PF-02, PF-03, PF-14, PF-17, PF-27 (last CDS exon before a UTR-only exon)
    Case(
        "ptc_in_a_cds_inside_exon_2_of_3_measures_to_the_junction_in_the_3utr",
        """
        C>T at CDS 15 makes CAG the PTC TAG. The CDS lies inside exon 2. So the 5' UTR spans the intron after exon 1,
        and the 3' UTR spans the intron before exon 3. utr5_length is cds_start_in_transcript, not a genomic distance.
        The PTC exon 2 is the last CDS exon, and exon 3 holds 3' UTR only. So ptc_to_intron runs to the junction in
        the 3' UTR, not to the CDS end. ptc_exon_length counts the UTR parts of exon 2.
        exon 1: 20 nt, exon 2: 70 nt, exon 3: 30 nt
        upstream_exon_count = 1, downstream_exon_count = 1, stop_codon_distance = 27 - 15 = 12.

        tx      0               20              35                  50                  65               90            120
        ref 5' [cccc..13..ccc]|[ccccccccccccccc ATG GCC GCC GCC GCC CAG GCC GCC GCC TAA cccc..17..cccc]|[cccc..22..cccc] 3'
        alt 5' [cccc..13..ccc]|[ccccccccccccccc ATG GCC GCC GCC GCC TAG GCC GCC GCC TAA cccc..17..cccc]|[cccc..22..cccc] 3'
                                                                    ^ C>T
        CDS     -35             -15             0                   15                  30               55            85
                                                                    *** PTC
                <----- utr5_length = 35 ------>
                                                <-----------------> ptc_to_start_codon = 15
                                                                    <------ ptc_to_intron = 40 ------>
                                                                    <-----------------> not to the CDS end
                                                                                        <----- utr3_length = 55 ------>
                                <----------------------- ptc_exon_length = 70 ----------------------->
        """,
        CDS_INSIDE_EXON_2,
        Change("GCC[C>T]AGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "start_variant": per_strand(80, 99),
            "end_variant": per_strand(81, 100),
            "alt_cds_start": per_strand(65, 85),
            "alt_cds_stop": per_strand(95, 115),
            "alt_cds_seq": "ATG" + "GCC" * 4 + "TAG" + "GCC" * 3 + "TAA",
            "alt_cds_len": 30,
            "alt_cds_info": [(2, 30)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 2,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 15,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(15, "TAG"), (27, "TAA")],
            "alt_stop_codon_exons": [2, 2],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "C" * 35 + "ATG" + "GCC" * 4 + "TAG" + "GCC" * 3 + "TAA" + "C" * 55,
            "alt_transcript_length": 120,
            "alt_cds_start_in_transcript": 35,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 15,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 70,
            "stop_codon_distance": 12,
            "ptc_to_intron": 40,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": True,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Ruler((-35, -15, 0, 15, 30, 55, 85), "CDS"),
            Mark("alt", 50, 53, "*", "PTC"),
            Span("alt", 0, 35, "utr5_length = 35"),
            Span("alt", 35, 50, "ptc_to_start_codon = 15"),
            Span("alt", 50, 90, "ptc_to_intron = 40"),
            Span("alt", 50, 65, "not to the CDS end"),
            Span("alt", 65, 120, "utr3_length = 55"),
            Span("alt", 20, 90, "ptc_exon_length = 70"),
        ),
        ruler=Ruler((0, 20, 35, 50, 65, 90, 120)),
    ),
    # PF-18, PF-27 (last exon)
    Case(
        "ptc_in_the_last_exon_measures_to_the_transcript_end",
        """
        C>T at CDS 30 makes CAG the PTC TAG, in the last exon. ptc_to_intron runs to the transcript end: it is the
        length of the 3' UTR that the PTC creates.
        exon 1: 20 nt, exon 2: 15 nt, exon 3: 55 nt
        upstream_exon_count = 2, downstream_exon_count = 0.

        tx      0        8                 20                    35  38                                  65            90
        ref 5' [cccccccc ATG GCC GCC GCC]|[GCC AAG GCC GCC GCC]|[GCC CAG GCC CTA ACA GCC CCT AAC GCC TAA gccg..17..cccc] 3'
        alt 5' [cccccccc ATG GCC GCC GCC]|[GCC AAG GCC GCC GCC]|[GCC TAG GCC CTA ACA GCC CCT AAC GCC TAA gccg..17..cccc] 3'
                                                                     ^ C>T
        CDS     -8       0                 12                    27  30                                  57            82
                                                                     *** PTC
                         <-------- ptc_to_start_codon = 30 -------->
                                                                     <-------------- ptc_to_intron = 52 -------------->
                                                                     <-----------------------------> stop_codon_distance = 24 (the stop codon at CDS 54)
                                                                 <--------------- ptc_exon_length = 55 --------------->
        """,
        THREE_EXONS,
        Change("GCC[C>T]AGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "start_variant": per_strand(88, 61),
            "end_variant": per_strand(89, 62),
            "alt_cds_start": per_strand(18, 35),
            "alt_cds_stop": per_strand(115, 132),
            "alt_cds_seq": "ATGGCCGCCGCC" + "GCCAAGGCCGCCGCC" + "GCCTAGGCCCTAACAGCCCCTAACGCCTAA",
            "alt_cds_len": 57,
            "alt_cds_info": [(1, 12), (2, 15), (3, 30)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 30,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(30, "TAG"), (54, "TAA")],
            "alt_stop_codon_exons": [3, 3],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "C" * 8
            + "ATGGCCGCCGCC"
            + "GCCAAGGCCGCCGCC"
            + "GCCTAGGCCCTAACAGCCCCTAACGCCTAA"
            + "GCCGCCTAG"
            + "C" * 16,
            "alt_transcript_length": 90,
            "alt_cds_start_in_transcript": 8,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 2,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 30,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 55,
            "stop_codon_distance": 24,
            "ptc_to_intron": 52,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Ruler((-8, 0, 12, 27, 30, 57, 82), "CDS"),
            Mark("alt", 38, 41, "*", "PTC"),
            Span("alt", 8, 38, "ptc_to_start_codon = 30"),
            Span("alt", 38, 90, "ptc_to_intron = 52"),
            Span("alt", 38, 62, "stop_codon_distance = 24 (the stop codon at CDS 54)"),
            Span("alt", 35, 90, "ptc_exon_length = 55"),
        ),
        ruler=Ruler((0, 8, 20, 35, 38, 65, 90)),
    ),
    # PF-20 (deletion)
    Case(
        "deletion_upstream_of_a_ptc_in_the_last_exon_moves_the_exon_end_in_alt_cds_coordinates",
        """
        Deleting one A of AAG at CDS 15, in exon 2, shifts the frame. The TAA at ref CDS 37 becomes the PTC, at CDS 36
        in alt CDS coordinates. The transcript end moves with it, from CDS 82 to 81: ptc_to_intron = 81 - 36 = 45.
        ref exon 1: 20 nt, exon 2: 15 nt, exon 3: 55 nt
        upstream_exon_count = 2, downstream_exon_count = 0, ptc_to_start_codon = 36, ptc_exon_length = 55.

        CDS         -8       0                 12                    27                                      57            82
        ref     5' [cccccccc ATG GCC GCC GCC]|[GCC AAG GCC GCC GCC]|[GCC CAG GCC CTA ACA GCC CCT AAC GCC TAA gccg..17..cccc] 3'
        alt     5' [cccccccc ATG GCC GCC GCC]|[GCC -AG GCC GCC GCC]|[GCC CAG GCC CTA ACA GCC CCT AAC GCC TAA gccg..17..cccc] 3'
                                                   ^ A>-
        alt CDS     -8       0                 12                    26           36                         56            81
                                                                                  **** PTC
                                                                                  <--------- ptc_to_intron = 45 ---------->
                                                                                  <--------------------> stop_codon_distance = 17 (the stop codon at alt CDS 53)
        """,
        THREE_EXONS,
        Change("GCC[A>]AGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CA", "TT"),
            "alt": per_strand("C", "T"),
            "start_variant": per_strand(52, 95),
            "end_variant": per_strand(54, 97),
            "alt_cds_start": per_strand(18, 35),
            "alt_cds_stop": per_strand(115, 132),
            "alt_cds_seq": "ATGGCCGCCGCC" + "GCCAGGCCGCCGCC" + "GCCCAGGCCCTAACAGCCCCTAACGCCTAA",
            "alt_cds_len": 56,
            "alt_cds_info": [(1, 12), (2, 14), (3, 30)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 36,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(36, "TAA")],
            "alt_stop_codon_exons": [3],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "C" * 8
            + "ATGGCCGCCGCC"
            + "GCCAGGCCGCCGCC"
            + "GCCCAGGCCCTAACAGCCCCTAACGCCTAA"
            + "GCCGCCTAG"
            + "C" * 16,
            "alt_transcript_length": 89,
            "alt_cds_start_in_transcript": 8,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 2,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 36,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 55,
            "stop_codon_distance": 17,
            "ptc_to_intron": 45,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exon_info": [(1, 20), (2, 14), (3, 55)],
            "nmd_model_status": "ok",
        },
        equivalent=(Change("GCCA[A>]GGCC"),),
        marks=(
            Ruler((-8, 0, 12, 26, 36, 56, 81), "CDS", "alt"),
            Mark("alt", 44, 47, "*", "PTC"),
            Span("alt", 44, 89, "ptc_to_intron = 45"),
            Span("alt", 44, 61, "stop_codon_distance = 17 (the stop codon at alt CDS 53)"),
        ),
        ruler=Ruler((-8, 0, 12, 27, 57, 82), "CDS"),
    ),
    # PF-20 (insertion)
    Case(
        "insertion_upstream_of_a_ptc_in_the_last_exon_moves_the_exon_end_in_alt_cds_coordinates",
        """
        Inserting an A into AAG at CDS 15, in exon 2, shifts the frame. The TAA at ref CDS 47 becomes the PTC, at CDS
        48 in alt CDS coordinates. The transcript end moves with it, from CDS 82 to 83: ptc_to_intron = 83 - 48 = 35.
        ref exon 1: 20 nt, exon 2: 15 nt, exon 3: 55 nt
        upstream_exon_count = 2, downstream_exon_count = 0, ptc_to_start_codon = 48, ptc_exon_length = 55.

        CDS         -8       0                 12                     27                                      57            82
        ref     5' [cccccccc ATG GCC GCC GCC]|[GCC -AAG GCC GCC GCC]|[GCC CAG GCC CTA ACA GCC CCT AAC GCC TAA gccg..17..cccc] 3'
        alt     5' [cccccccc ATG GCC GCC GCC]|[GCC AAAG GCC GCC GCC]|[GCC CAG GCC CTA ACA GCC CCT AAC GCC TAA gccg..17..cccc] 3'
                                                   ^ ->A
        alt CDS     -8       0                 12                     28                        48            58            83
                                                                                                **** PTC
                                                                                                <--- ptc_to_intron = 35 --->
                                                                                                <-------> stop_codon_distance = 7 (the stop codon at alt CDS 55)
        """,
        THREE_EXONS,
        Change("GCC[>A]AAGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "T"),
            "alt": per_strand("CA", "TT"),
            "start_variant": per_strand(52, 96),
            "end_variant": per_strand(53, 97),
            "alt_cds_start": per_strand(18, 35),
            "alt_cds_stop": per_strand(115, 132),
            "alt_cds_seq": "ATGGCCGCCGCC" + "GCCAAAGGCCGCCGCC" + "GCCCAGGCCCTAACAGCCCCTAACGCCTAA",
            "alt_cds_len": 58,
            "alt_cds_info": [(1, 12), (2, 16), (3, 30)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 48,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(48, "TAA")],
            "alt_stop_codon_exons": [3],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "C" * 8
            + "ATGGCCGCCGCC"
            + "GCCAAAGGCCGCCGCC"
            + "GCCCAGGCCCTAACAGCCCCTAACGCCTAA"
            + "GCCGCCTAG"
            + "C" * 16,
            "alt_transcript_length": 91,
            "alt_cds_start_in_transcript": 8,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 2,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 48,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 55,
            "stop_codon_distance": 7,
            "ptc_to_intron": 35,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exon_info": [(1, 20), (2, 16), (3, 55)],
            "nmd_model_status": "ok",
        },
        equivalent=(Change("GCCA[>A]AGGCC"), Change("GCCAA[>A]GGCC")),
        marks=(
            Ruler((-8, 0, 12, 28, 48, 58, 83), "CDS", "alt"),
            Mark("alt", 56, 59, "*", "PTC"),
            Span("alt", 56, 91, "ptc_to_intron = 35"),
            Span("alt", 56, 63, "stop_codon_distance = 7 (the stop codon at alt CDS 55)"),
        ),
        ruler=Ruler((-8, 0, 12, 27, 57, 82), "CDS"),
    ),
    # Pins the rule of upstream_exon_count "an exon that the variant deletes is not in the mRNA and does not count".
    # The deleted exon 2 has length 0 in alt_transcript_exon_info and lies upstream of the PTC exon. The other case
    # with a deleted exon (coding_region_edges) deletes it in frame, so its row has no PTC and no exon counts.
    Case(
        "deletion_of_a_whole_exon_upstream_of_the_ptc_exon_is_no_upstream_exon",
        """
        The deletion removes all 13 nt of exon 2, and the acceptor AG and the donor GT stay. It has one placement:
        the G before it differs from its last base C, and the G after it differs from its first base C. So exon 2
        has length 0 in alt_transcript_exon_info. The equivalent description GCAGCAGCTGCTGC>G is a delins. Matched
        from the right, it deletes the G of the AG and loses the acceptor. So only the match from the left is valid.
        The deletion shifts the frame, and exon 3 starts with the PTC TGA, at alt CDS 15 and alt tx 20. An exon that
        the variant deletes is not in the mRNA and does not count. So upstream_exon_count = 1 (exon 1), and
        downstream_exon_count = 0. The PTC lies in the last exon, so ptc_to_intron runs to the transcript end.
        ref exon 1: 20 nt, exon 2: 13 nt, exon 3: 20 nt

        tx         0     5                                        20                  33                      53
        ref    5' [gccac ATG GCC AAG CTG CTG]gtaagtcccccccctttcag[CAG CAG CTG CTG C]|[TG AAG CTG TAA gcccccccc] 3'
        alt    5' [gccac ATG GCC AAG CTG CTG]gtaagtcccccccctttcag[--- --- --- --- -]|[TG AAG CTG TAA gcccccccc] 3'
                                                                  ^^^^^^^^^^^^^^^^^ 13 nt>-
        alt tx     0     5                                                            20         28           40
                                                                                      **** PTC
                                                                                      <--------> stop_codon_distance = 8
                                                                                      <----------------------> ptc_to_intron = 20
        """,
        EXON_2_OF_13_NT,
        Change("TTTCAG[CAGCAGCTGCTGC>]GTAAGT"),
        {
            "variant_id": "var1",
            "ref": per_strand("GCAGCAGCTGCTGC", "CGCAGCAGCTGCTG"),
            "alt": per_strand("G", "C"),
            "start_variant": per_strand(49, 49),
            "end_variant": per_strand(63, 63),
            "alt_cds_start": per_strand(15, 19),
            "alt_cds_stop": per_strand(94, 98),
            "alt_cds_seq": "ATGGCCAAGCTGCTG" + "TGAAGCTGTAA",
            "alt_cds_len": 26,
            "alt_cds_info": [(1, 15), (2, 0), (3, 11)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 15,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(15, "TGA")],
            "alt_stop_codon_exons": [3],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GCCACATGGCCAAGCTGCTG" + "TGAAGCTGTAAGCCCCCCCC",
            "alt_transcript_length": 40,
            "alt_cds_start_in_transcript": 5,
            "alt_transcript_exon_info": [(1, 20), (2, 0), (3, 20)],
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 15,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 20,
            "stop_codon_distance": 8,
            "ptc_to_intron": 20,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "nmd_model_status": "ok",
        },
        equivalent=(Change("TTTCA[GCAGCAGCTGCTGC>G]GTAAGT"),),
        marks=(
            Ruler((0, 5, 20, 28, 40), "tx", "alt"),
            Mark("alt", 20, 23, "*", "PTC"),
            Span("alt", 20, 28, "stop_codon_distance = 8"),
            Span("alt", 20, 40, "ptc_to_intron = 20"),
        ),
        ruler=Ruler((0, 5, 20, 33, 53)),
    ),
    # Pins the threshold of the same rule: an exon of 1 nt is in the mRNA, so it counts as an upstream exon. Only a
    # deleted exon, of length 0, does not count. No other case has an exon of 1 nt.
    Case(
        "exon_of_1_nt_upstream_of_the_ptc_exon_is_an_upstream_exon",
        """
        Exon 1 holds 1 nt of 5' UTR. G>A at CDS 8 makes TGG the PTC TGA, at CDS 6 and tx 10, in exon 2, the last
        exon. Exon 1 has 1 base in the mRNA, so it is an upstream exon: upstream_exon_count = 1, and
        downstream_exon_count = 0.
        exon 1: 1 nt, exon 2: 21 nt

        tx      0   1   4       10      16     22
        ref 5' [g]|[acc ATG GCC TGG AAG TAA gcc] 3'
        alt 5' [g]|[acc ATG GCC TGA AAG TAA gcc] 3'
                                  ^ G>A
                                *** PTC
                                <-----> stop_codon_distance = 6
                                <-------------> ptc_to_intron = 12
        """,
        UTR_EXON_1_OF_1_NT,
        Change("GCCTG[G>A]AAGTAA"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(42, 19),
            "end_variant": per_strand(43, 20),
            "alt_cds_start": per_strand(34, 13),
            "alt_cds_stop": per_strand(49, 28),
            "alt_cds_seq": "ATGGCCTGAAAGTAA",
            "alt_cds_len": 15,
            "alt_cds_info": [(2, 15)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 2,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 6,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(6, "TGA"), (12, "TAA")],
            "alt_stop_codon_exons": [2, 2],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "G" + "ACCATGGCCTGAAAGTAAGCC",
            "alt_transcript_length": 22,
            "alt_cds_start_in_transcript": 4,
            "alt_transcript_exon_info": SAME_EXONS,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 6,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 21,
            "stop_codon_distance": 6,
            "ptc_to_intron": 12,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "nmd_model_status": "ok",
        },
        marks=(
            Mark("alt", 10, 13, "*", "PTC"),
            Span("alt", 10, 16, "stop_codon_distance = 6"),
            Span("alt", 10, 22, "ptc_to_intron = 12"),
        ),
        ruler=Ruler((0, 1, 4, 10, 16, 22)),
    ),
    # PF-15 (negative)
    Case(
        "stop_loss_has_a_negative_stop_codon_distance",
        """
        T>C at CDS 54 turns the stop codon TAA into CAA. The alt transcript reads on into the 3' UTR, in the frame of
        the CDS, to the TAG at CDS 63. This stop codon lies 9 nt downstream of the annotated one, so
        stop_codon_distance = 54 - 63 = -9.
        exon 1: 20 nt, exon 2: 15 nt, exon 3: 55 nt

        tx      0        8                 20                    35                     62  65    71            90
        ref 5' [cccccccc ATG GCC GCC GCC]|[GCC AAG GCC GCC GCC]|[GCC CAG ..15.. AAC GCC TAA gccgcctagc..11..cccc] 3'
        alt 5' [cccccccc ATG GCC GCC GCC]|[GCC AAG GCC GCC GCC]|[GCC CAG ..15.. AAC GCC CAA gccgcctagc..11..cccc] 3'
                                                                                        ^ T>C
        CDS     -8       0                 12                    27                     54  57    63            82
                                                                                        <--------> stop_codon_distance = -9
        """,
        THREE_EXONS,
        Change("CGCC[T>C]AAGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "start_variant": per_strand(112, 37),
            "end_variant": per_strand(113, 38),
            "alt_cds_start": per_strand(18, 35),
            "alt_cds_stop": per_strand(115, 132),
            "alt_cds_seq": "ATGGCCGCCGCC" + "GCCAAGGCCGCCGCC" + "GCCCAGGCCCTAACAGCCCCTAACGCCCAA",
            "alt_cds_len": 57,
            "alt_cds_info": [(1, 12), (2, 15), (3, 30)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "CAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": 0,
            "alt_all_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": "C" * 8
            + "ATGGCCGCCGCC"
            + "GCCAAGGCCGCCGCC"
            + "GCCCAGGCCCTAACAGCCCCTAACGCCCAA"
            + "GCCGCCTAG"
            + "C" * 16,
            "alt_transcript_length": 90,
            "alt_cds_start_in_transcript": 8,
            "transcript_start_codon_pos": 8,
            "transcript_start_codon_exon": 1,
            "transcript_last_codon": "CCC",
            "transcript_valid_stop": False,
            "transcript_first_stop_codon": "TAG",
            "transcript_first_stop_pos": 71,
            "transcript_num_stop_codons": 1,
            "transcript_all_stop_codons": [(71, "TAG")],
            "transcript_stop_codon_exons": [3],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": -9,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        marks=(Ruler((-8, 0, 12, 27, 54, 57, 63, 82), "CDS"), Span("alt", 62, 71, "stop_codon_distance = -9")),
        ruler=Ruler((0, 8, 20, 35, 62, 65, 71, 90)),
    ),
    # PF-06, PF-26
    Case(
        "splice_site_snv_keeps_the_exon_count_and_the_utr_lengths_of_the_ref",
        """
        G>A destroys the donor GT after exon 1, so the alt transcript is unknown. The row keeps the features of the ref
        transcript: total_exon_count = 3, utr5_length = 8, utr3_length = 25 and likely_misannotated False. The features
        of the alt transcript are null.
        exon 1: 20 nt, exon 2: 15 nt, exon 3: 55 nt

        tx      0        8                                    20                    35                 65            90
        ref 5' [cccccccc ATG GCC GCC GCC]gtaagtcccccccctttcag[GCC AAG GCC GCC GCC]|[GCC CAG ..21.. TAA gccg..17..cccc] 3'
        alt 5' [cccccccc ATG GCC GCC GCC]ataagtcccccccctttcag[GCC AAG GCC GCC GCC]|[GCC CAG ..21.. TAA gccg..17..cccc] 3'
                                         ^ G>A
        """,
        THREE_EXONS,
        Change("ATGGCCGCCGCC[G>A]T"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(30, 119),
            "end_variant": per_strand(31, 120),
            **UNKNOWN_ALT,
            "unknown_reason": "splice_site_destroyed",
        },
        ruler=Ruler((0, 8, 20, 35, 65, 90)),
    ),
    # PF-19, PF-27 (single exon)
    Case(
        "ptc_in_a_single_exon_transcript_measures_to_the_transcript_end",
        """
        C>T at CDS 30 makes CAG the PTC TAG. The transcript has one exon. As in the last exon, ptc_to_intron runs to
        the transcript end.
        ptc_exon_length = 90. upstream_exon_count = 0, downstream_exon_count = 0.

        tx      0        8                      38                                  65            90
        ref 5' [cccccccc ATG GCC ..18.. GCC GCC CAG GCC CTA ACA GCC CCT AAC GCC TAA gccg..17..cccc] 3'
        alt 5' [cccccccc ATG GCC ..18.. GCC GCC TAG GCC CTA ACA GCC CCT AAC GCC TAA gccg..17..cccc] 3'
                                                ^ C>T
        CDS     -8       0                      30                                  57            82
                                                *** PTC
                         <--------------------> ptc_to_start_codon = 30
                                                <-------------- ptc_to_intron = 52 -------------->
                                                <-----------------------------> stop_codon_distance = 24 (the stop codon at CDS 54)
        """,
        SINGLE_EXON,
        Change("GCC[C>T]AGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "start_variant": per_strand(48, 61),
            "end_variant": per_strand(49, 62),
            "alt_cds_start": per_strand(18, 35),
            "alt_cds_stop": per_strand(75, 92),
            "alt_cds_seq": "ATGGCCGCCGCC" + "GCCAAGGCCGCCGCC" + "GCCTAGGCCCTAACAGCCCCTAACGCCTAA",
            "alt_cds_len": 57,
            "alt_cds_info": [(1, 57)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 30,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(30, "TAG"), (54, "TAA")],
            "alt_stop_codon_exons": [1, 1],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "C" * 8
            + "ATGGCCGCCGCC"
            + "GCCAAGGCCGCCGCC"
            + "GCCTAGGCCCTAACAGCCCCTAACGCCTAA"
            + "GCCGCCTAG"
            + "C" * 16,
            "alt_transcript_length": 90,
            "alt_cds_start_in_transcript": 8,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 0,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 30,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 90,
            "stop_codon_distance": 24,
            "ptc_to_intron": 52,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Ruler((-8, 0, 30, 57, 82), "CDS"),
            Mark("alt", 38, 41, "*", "PTC"),
            Span("alt", 8, 38, "ptc_to_start_codon = 30"),
            Span("alt", 38, 90, "ptc_to_intron = 52"),
            Span("alt", 38, 62, "stop_codon_distance = 24 (the stop codon at CDS 54)"),
        ),
        ruler=Ruler((0, 8, 38, 65, 90)),
    ),
    # PF-04, PF-15 (zero), PF-13 (no PTC row)
    Case(
        "missense_in_a_cds_at_the_transcript_start_has_no_5utr",
        """
        The CDS starts at the transcript start, so utr5_length is 0. The missense AAG>GAG keeps the annotated stop
        codon: stop_codon_distance is 0. Not a PTC row, so ptc_less_than_150nt_to_start is False.
        exon 1: 12 nt, exon 2: 15 nt

        tx      0                 12          21    27
        ref 5' [ATG GCC AAG CTG]|[GGC TCC TAA gccgcc] 3'
        alt 5' [ATG GCC GAG CTG]|[GGC TCC TAA gccgcc] 3'
                        ^ A>G
        """,
        NO_UTR5,
        Change("GCC[A>G]AGCTG"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("G", "C"),
            "start_variant": per_strand(16, 50),
            "end_variant": per_strand(17, 51),
            "alt_cds_start": per_strand(10, 16),
            "alt_cds_stop": per_strand(51, 57),
            "alt_cds_seq": "ATGGCCGAGCTGGGCTCCTAA",
            "alt_cds_len": 21,
            "alt_cds_info": [(1, 12), (2, 9)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 18,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(18, "TAA")],
            "alt_stop_codon_exons": [2],
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "ATGGCCGAGCTGGGCTCCTAAGCCGCC",
            "alt_transcript_length": 27,
            "alt_cds_start_in_transcript": 0,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 12, 21, 27)),
    ),
    # PF-05, PF-15 (zero)
    Case(
        "missense_in_a_cds_at_the_transcript_end_has_no_3utr",
        """
        The stop codon ends the transcript, so utr3_length is 0. The missense AAG>GAG keeps the annotated stop codon:
        stop_codon_distance is 0.
        exon 1: 16 nt, exon 2: 9 nt

        tx      0    4                 16         25
        ref 5' [gacc ATG GCC AAG CTG]|[GGC TCC TAA] 3'
        alt 5' [gacc ATG GCC GAG CTG]|[GGC TCC TAA] 3'
                             ^ A>G
        """,
        NO_UTR3,
        Change("GCC[A>G]AGCTG"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("G", "C"),
            "start_variant": per_strand(20, 44),
            "end_variant": per_strand(21, 45),
            "alt_cds_start": per_strand(14, 10),
            "alt_cds_stop": per_strand(55, 51),
            "alt_cds_seq": "ATGGCCGAGCTGGGCTCCTAA",
            "alt_cds_len": 21,
            "alt_cds_info": [(1, 12), (2, 9)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 18,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(18, "TAA")],
            "alt_stop_codon_exons": [2],
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATGGCCGAGCTGGGCTCCTAA",
            "alt_transcript_length": 25,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 4, 16, 25)),
    ),
    # PF-12 (150 nt)
    Case(
        "ptc_150_nt_after_the_start_codon_is_not_less_than_150_nt_to_start",
        """
        C>T at CDS 150 makes CAG the PTC TAG, 150 nt from the start codon ATG at CDS 0. 150 is not < 150, so
        ptc_less_than_150nt_to_start is False. No NMD rule is True.
        exon 1: 214 nt, exon 2: 12 nt
        ptc_exon_length = 214. upstream_exon_count = 0, downstream_exon_count = 1.

        tx      0    4                       154                          214     220   226
        ref 5' [gacc ATG GCC ..138.. GCC GCC CAG GCC GCC ..45.. GCC GCC]|[GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..138.. GCC GCC TAG GCC GCC ..45.. GCC GCC]|[GCC TAA gccgcc] 3'
                                             ^ C>T
        CDS     -4   0                       150                          210     216
                                             *** PTC
                     <---------------------> ptc_to_start_codon = 150
                                             <-- ptc_to_intron = 60 -->
                                             <------------------------------> stop_codon_distance = 63 (the stop codon at CDS 213)
        """,
        PTC_150,
        Change("GCC[C>T]AGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "start_variant": per_strand(164, 101),
            "end_variant": per_strand(165, 102),
            "alt_cds_start": per_strand(14, 16),
            "alt_cds_stop": per_strand(250, 252),
            "alt_cds_seq": "ATG" + "GCC" * 49 + "TAG" + "GCC" * 20 + "TAA",
            "alt_cds_len": 216,
            "alt_cds_info": [(1, 210), (2, 6)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 150,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(150, "TAG"), (213, "TAA")],
            "alt_stop_codon_exons": [1, 2],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACC" + "ATG" + "GCC" * 49 + "TAG" + "GCC" * 20 + "TAA" + "GCCGCC",
            "alt_transcript_length": 226,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 0,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 150,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 214,
            "stop_codon_distance": 63,
            "ptc_to_intron": 60,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Ruler((-4, 0, 150, 210, 216), "CDS"),
            Mark("alt", 154, 157, "*", "PTC"),
            Span("alt", 4, 154, "ptc_to_start_codon = 150"),
            Span("alt", 154, 214, "ptc_to_intron = 60"),
            Span("alt", 154, 217, "stop_codon_distance = 63 (the stop codon at CDS 213)"),
        ),
        ruler=Ruler((0, 4, 154, 214, 220, 226)),
    ),
    # PF-10, PF-22 (CTG start)
    Case(
        "ptc_after_a_ctg_start_codon_is_measured_from_the_ctg_not_from_an_internal_atg",
        """
        The annotated start codon is CTG at CDS 0. The ATG at CDS 30 is an internal Met. G>A at CDS 160 makes TGG the
        PTC TAG, at CDS 159: 159 nt from the start codon CTG, but only 129 nt from the ATG. ptc_to_start_codon is
        measured from the annotated start codon, so it is not < 150. A CTG start codon is no misannotation:
        likely_misannotated is False.
        exon 1: 103 nt, exon 2: 73 nt
        ptc_exon_length = 73. upstream_exon_count = 1, downstream_exon_count = 0.

        tx      0   3                  33                     103                   162     168      176
        ref 5' [ggg CTG AAA ..21.. AAA ATG AAA ..60.. AAA A]|[AA AAA ..48.. AAA AAA TGG GAC TAA ggggg] 3'
        alt 5' [ggg CTG AAA ..21.. AAA ATG AAA ..60.. AAA A]|[AA AAA ..48.. AAA AAA TAG GAC TAA ggggg] 3'
                                                                                     ^ G>A
        CDS     -3  0                  30                     100                   159     165      173
                                                                                    *** PTC
                    <----------------- ptc_to_start_codon = 159 ------------------>
                                                                                    <---------------> ptc_to_intron = 14
                                                                                    <-----> stop_codon_distance = 6
        """,
        CTG_START,
        Change("AAAT[G>A]GGACTAA"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(193, 22),
            "end_variant": per_strand(194, 23),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(201, 203),
            "alt_cds_seq": "CTG" + "AAA" * 9 + "ATG" + "AAA" * 42 + "TAG" + "GAC" + "TAA",
            "alt_cds_len": 168,
            "alt_cds_info": [(1, 100), (2, 68)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 159,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(159, "TAG"), (165, "TAA")],
            "alt_stop_codon_exons": [2, 2],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GGG" + "CTG" + "AAA" * 9 + "ATG" + "AAA" * 42 + "TAG" + "GAC" + "TAA" + "GGGGG",
            "alt_transcript_length": 176,
            "alt_cds_start_in_transcript": 3,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 159,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 73,
            "stop_codon_distance": 6,
            "ptc_to_intron": 14,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Ruler((-3, 0, 30, 100, 159, 165, 173), "CDS"),
            Mark("alt", 162, 165, "*", "PTC"),
            Span("alt", 3, 162, "ptc_to_start_codon = 159"),
            Span("alt", 162, 176, "ptc_to_intron = 14"),
            Span("alt", 162, 168, "stop_codon_distance = 6"),
        ),
        ruler=Ruler((0, 3, 33, 103, 162, 168, 176)),
    ),
    # PF-13 (no annotated start codon), PF-23; the same instance as NU-22
    Case(
        "ptc_in_a_cds_start_nf_transcript_has_no_distance_to_the_start_codon",
        """
        The transcript of the CTG case, without start_codon rows and tagged cds_start_NF. Its start codon lies
        upstream of the CDS, at an unknown distance. So ptc_to_start_codon is null, and ptc_less_than_150nt_to_start
        is False. Without an annotated start codon, likely_misannotated is True.
        exon 1: 103 nt, exon 2: 73 nt
        ptc_exon_length = 73. upstream_exon_count = 1, downstream_exon_count = 0.

        tx      0   3                      103                   162     168      176
        ref 5' [ggg CTG AAA ..90.. AAA A]|[AA AAA ..48.. AAA AAA TGG GAC TAA ggggg] 3'
        alt 5' [ggg CTG AAA ..90.. AAA A]|[AA AAA ..48.. AAA AAA TAG GAC TAA ggggg] 3'
                                                                  ^ G>A
        CDS     -3  0                      100                   159     165      173
                                                                 *** PTC
                                                                 <---------------> ptc_to_intron = 14
        """,
        CDS_START_NF,
        Change("AAAT[G>A]GGACTAA"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(193, 22),
            "end_variant": per_strand(194, 23),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(201, 203),
            "alt_cds_seq": "CTG" + "AAA" * 9 + "ATG" + "AAA" * 42 + "TAG" + "GAC" + "TAA",
            "alt_cds_len": 168,
            "alt_cds_info": [(1, 100), (2, 68)],
            "alt_start_codon_pos": None,
            "alt_start_codon_exon": None,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 159,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(159, "TAG"), (165, "TAA")],
            "alt_stop_codon_exons": [2, 2],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GGG" + "CTG" + "AAA" * 9 + "ATG" + "AAA" * 42 + "TAG" + "GAC" + "TAA" + "GGGGG",
            "alt_transcript_length": 176,
            "alt_cds_start_in_transcript": 3,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": None,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 73,
            "stop_codon_distance": 6,
            "ptc_to_intron": 14,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_annotated_start",
        },
        marks=(
            Ruler((-3, 0, 100, 159, 165, 173), "CDS"),
            Mark("alt", 162, 165, "*", "PTC"),
            Span("alt", 162, 176, "ptc_to_intron = 14"),
        ),
        ruler=Ruler((0, 3, 103, 162, 168, 176)),
    ),
    # PF-24 (no stop_codon rows)
    Case(
        "missense_without_stop_codon_rows_is_likely_misannotated",
        """
        The transcript has no stop_codon rows and is tagged cds_end_NF. The CDS ends in the sense codon TCC at the
        transcript end. has_stop_codon is False, so ref_valid_stop is False, and likely_misannotated is True. Without
        an annotated stop codon, utr3_length and stop_codon_distance are null.
        exon 1: 16 nt, exon 2: 6 nt

        tx      0    4                 16     22
        ref 5' [gacc ATG GCC AAG CTG]|[GGC TCC] 3'
        alt 5' [gacc ATG GCC GAG CTG]|[GGC TCC] 3'
                             ^ A>G
        """,
        CDS_END_NF,
        Change("GCC[A>G]AGCTG"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("G", "C"),
            "start_variant": per_strand(20, 41),
            "end_variant": per_strand(21, 42),
            "alt_cds_start": per_strand(14, 10),
            "alt_cds_stop": per_strand(52, 48),
            "alt_cds_seq": "ATGGCCGAGCTGGGCTCC",
            "alt_cds_len": 18,
            "alt_cds_info": [(1, 12), (2, 6)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TCC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": 0,
            "alt_all_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATGGCCGAGCTGGGCTCC",
            "alt_transcript_length": 22,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 4, 16, 22)),
    ),
    # PF-24 (stop_codon rows on a sense codon)
    Case(
        "missense_with_stop_codon_rows_on_a_sense_codon_is_likely_misannotated",
        """
        The stop_codon rows lie on the sense codon TCC. has_stop_codon is True, but ref_valid_stop is False, so
        likely_misannotated is True. The row keeps the flags from the CDS, whose alt CDS has no in-frame stop codon,
        so stop_codon_distance is null.
        exon 1: 16 nt, exon 2: 12 nt

        tx      0    4                 16      22    28
        ref 5' [gacc ATG GCC AAG CTG]|[GGC TCC gccgcc] 3'
        alt 5' [gacc ATG GCC GAG CTG]|[GGC TCC gccgcc] 3'
                             ^ A>G
        """,
        STOP_CODON_ROWS_ON_A_SENSE_CODON,
        Change("GCC[A>G]AGCTG"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("G", "C"),
            "start_variant": per_strand(20, 47),
            "end_variant": per_strand(21, 48),
            "alt_cds_start": per_strand(14, 16),
            "alt_cds_stop": per_strand(52, 54),
            "alt_cds_seq": "ATGGCCGAGCTGGGCTCC",
            "alt_cds_len": 18,
            "alt_cds_info": [(1, 12), (2, 6)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TCC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": 0,
            "alt_all_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATGGCCGAGCTGGGCTCCGCCGCC",
            "alt_transcript_length": 28,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 4, 16, 22, 28)),
    ),
    # PF-24 (Ensembl, chrM, AGA)
    Case(
        "missense_in_a_mitochondrial_cds_ending_in_aga_is_likely_misannotated",
        """
        Ensembl GFF3 on chrM. The CDS ends in AGA, a mitochondrial stop codon, so has_stop_codon is True. The codon
        scans know only TAA, TAG and TGA. So ref_valid_stop is False, and likely_misannotated is True. The row keeps
        the flags from the CDS, whose alt CDS has no in-frame stop codon, so stop_codon_distance is null.
        exon 1: 18 nt

        tx      0                      18
        ref 5' [ATG GCC AAG CTG GGC AGA] 3'
        alt 5' [ATG GCC GAG CTG GGC AGA] 3'
                        ^ A>G
        """,
        MITOCHONDRIAL_AGA,
        Change("GCC[A>G]AGCTG"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("G", "C"),
            "start_variant": per_strand(16, 21),
            "end_variant": per_strand(17, 22),
            "alt_cds_start": per_strand(10, 10),
            "alt_cds_stop": per_strand(28, 28),
            "alt_cds_seq": "ATGGCCGAGCTGGGCAGA",
            "alt_cds_len": 18,
            "alt_cds_info": [(1, 18)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "AGA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": 0,
            "alt_all_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "ATGGCCGAGCTGGGCAGA",
            "alt_transcript_length": 18,
            "alt_cds_start_in_transcript": 0,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 18)),
    ),
]
