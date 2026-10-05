"""
Conformance cases of the NMD rules: each rule True and False, nmd_escape, and the thresholds of the 50 nt rule
(50/51 nt, and 0/1 nt at the last exon junction), the long exon rule (407/408 nt) and the start-proximal rule
(147/150 nt).

The drawings show the transcript 5' to 3', also for the minus strand. The runner renders their layout block from
the case. The line ref is the transcript and the line alt is the transcript with the change. The two lines are
aligned, and `-` fills a gap. `[...]` is an exon and `|` is an exon junction. Upper case is the coding region, in
codons of the annotated frame; lower case is UTR. An intron that the change touches is drawn in lower case outside
the brackets. `..N..` leaves out N bases. `^` marks the change as ref>alt. A mark line puts `***` under the PTC, and
`<-- label -->` spans a length, e.g. from the PTC to an exon junction or to the transcript end. The PTC mark stands
under the first base of the codon only where the drawing splits the codon: by an exon junction or the codon spacing
of another frame. The ruler CDS gives CDS positions, as alt_first_stop_pos, and the ruler alt CDS gives them for the
alt line. They count from the first coding base, so the 5' UTR has negative positions. Most CDS are built from GCC
repeats, so the codons that matter stand out.
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

# Layouts A: 3 exons, a CDS of 261 nt. Only the junction of exon 2 and exon 3 moves between them.
A_HEAD = "ATG" + "GCC" * 9 + "AAA" + "GCC" * 9 + "GAG"  # CDS 0 to 62, the CDS part of exon 1
A_MID = "GCC" * 28  # CDS 63 to 146
A_TAIL = "GCC" * 23 + "TGG" + "GCC" * 10 + "TAA"  # CDS 156 to 260, after CAG TTG AGC at 147 to 155
A_CDS = A_HEAD + A_MID + "CAGTTGAGC" + A_TAIL
A_REF = {
    **IDS,
    "ref_cds_start": per_strand(14, 16),
    "ref_cds_stop": per_strand(315, 317),
    "ref_cds_seq": A_CDS,
    "ref_cds_len": 261,
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "cds_in_transcript": True,
    "ref_start_codon_pos": 0,
    "ref_start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_first_stop_pos": 258,
    "ref_num_stop_codons": 1,
    "ref_all_stop_codons": [(258, "TAA")],
    "ref_stop_codon_exons": [3],
    "ref_is_premature": False,
    "transcript_start": 10,
    "transcript_end": 321,
    "transcript_seq": "GACC" + A_CDS + "GCCGCC",
    "transcript_length": 271,
    "cds_start_in_transcript": 4,
    "cds_end_in_transcript": 265,
    "utr3_length": 6,
    "utr5_length": 4,
    "total_exon_count": 3,
    "likely_misannotated": False,
}
# The last exon junction at CDS position 201
A201 = Layout(
    Transcript(("gacc" + A_HEAD, A_MID + "CAGTTGAGC" + "GCC" * 15, "GCC" * 8 + "TGG" + "GCC" * 10 + "TAA" + "gccgcc")),
    {**A_REF, "ref_cds_info": [(1, 63), (2, 138), (3, 60)], "transcript_exon_info": [(1, 67), (2, 138), (3, 66)]},
)
# The last exon junction at CDS position 200
A200 = Layout(
    Transcript(
        (
            "gacc" + A_HEAD,
            A_MID + "CAGTTGAGC" + "GCC" * 14 + "GC",
            "C" + "GCC" * 8 + "TGG" + "GCC" * 10 + "TAA" + "gccgcc",
        )
    ),
    {**A_REF, "ref_cds_info": [(1, 63), (2, 137), (3, 61)], "transcript_exon_info": [(1, 67), (2, 137), (3, 67)]},
)
# The last exon junction at CDS position 150, before TTG
A150 = Layout(
    Transcript(("gacc" + A_HEAD, A_MID + "CAG", "TTGAGC" + A_TAIL + "gccgcc")),
    {**A_REF, "ref_cds_info": [(1, 63), (2, 87), (3, 111)], "transcript_exon_info": [(1, 67), (2, 87), (3, 117)]},
)
# The last exon junction at CDS position 151, inside TTG
A151 = Layout(
    Transcript(("gacc" + A_HEAD, A_MID + "CAGT", "TGAGC" + A_TAIL + "gccgcc")),
    {**A_REF, "ref_cds_info": [(1, 63), (2, 88), (3, 110)], "transcript_exon_info": [(1, 67), (2, 88), (3, 116)]},
)
# The alt columns of an SNV in a layout A that keeps the start codon and the stop codon
A_SNV = {
    "variant_id": "var1",
    "alt_cds_start": per_strand(14, 16),
    "alt_cds_stop": per_strand(315, 317),
    "alt_cds_len": 261,
    "alt_start_codon_pos": 0,
    "alt_start_codon_exon": 1,
    "alt_last_codon": "TAA",
    "alt_valid_stop": True,
    "start_loss": False,
    "stop_loss": False,
    "alt_transcript_length": 271,
    "alt_cds_start_in_transcript": 4,
    **NOT_SCANNED,
    "unknown_reason": None,
}
# T>A at CDS position 151: TTG becomes the PTC TAG at 150
A_PTC_150 = Change("CAGT[T>A]GAGC")
A_PTC_150_SEQ = A_HEAD + A_MID + "CAGTAGAGC" + A_TAIL

# Layout SINGLE: the CDS of the layouts A in a single exon
SINGLE = Layout(
    Transcript(("gacc" + A_CDS + "gccgcc",)),
    {
        **A_REF,
        "ref_cds_start": per_strand(14, 16),
        "ref_cds_stop": per_strand(275, 277),
        "ref_stop_codon_exons": [1],
        "transcript_end": 281,
        "ref_cds_info": [(1, 261)],
        "transcript_exon_info": [(1, 271)],
        "total_exon_count": 1,
    },
)
SINGLE_SNV = {**A_SNV, "alt_cds_stop": per_strand(275, 277)}

# Layouts U: 3 exons, exon 3 holds only 3' UTR. Exon 2 ends at CDS position 240 (U51) or 239 (U50): the last exon
# junction lies in the 3' UTR.
U_CDS = "ATG" + "GCC" * 49 + "AAG" + "GCC" * 12 + "CAG" + "GCC" * 5 + "TAA"
U_EXON_1 = "gacc" + "ATG" + "GCC" * 49 + "AAG" + "GCC" * 9
U_REF = {
    **IDS,
    "ref_cds_start": per_strand(14, 72),
    "ref_cds_stop": per_strand(244, 302),
    "ref_cds_seq": U_CDS,
    "ref_cds_len": 210,
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "ref_cds_info": [(1, 180), (2, 30)],
    "cds_in_transcript": True,
    "ref_start_codon_pos": 0,
    "ref_start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_first_stop_pos": 207,
    "ref_num_stop_codons": 1,
    "ref_all_stop_codons": [(207, "TAA")],
    "ref_stop_codon_exons": [2],
    "ref_is_premature": False,
    "transcript_start": 10,
    "transcript_end": 306,
    "transcript_seq": "GACC" + U_CDS + "GCC" * 10 + "GCC" * 4,
    "transcript_length": 256,
    "cds_start_in_transcript": 4,
    "cds_end_in_transcript": 214,
    "transcript_exon_info": [(1, 184), (2, 60), (3, 12)],
    "utr3_length": 42,
    "utr5_length": 4,
    "total_exon_count": 3,
    "likely_misannotated": False,
}
U51 = Layout(Transcript((U_EXON_1, "GCC" * 3 + "CAG" + "GCC" * 5 + "TAA" + "gcc" * 10, "gcc" * 4)), U_REF)
U50 = Layout(
    Transcript((U_EXON_1, "GCC" * 3 + "CAG" + "GCC" * 5 + "TAA" + "gcc" * 9 + "gc", "gcc" * 4)),
    {
        **U_REF,
        "ref_cds_start": per_strand(14, 71),
        "ref_cds_stop": per_strand(244, 301),
        "transcript_end": 305,
        "transcript_seq": "GACC" + U_CDS + "GCC" * 9 + "GC" + "GCC" * 4,
        "transcript_length": 255,
        "transcript_exon_info": [(1, 184), (2, 59), (3, 12)],
        "utr3_length": 41,
    },
)
U_SNV = {
    "variant_id": "var1",
    "alt_cds_start": per_strand(14, 72),
    "alt_cds_stop": per_strand(244, 302),
    "alt_cds_len": 210,
    "alt_cds_info": [(1, 180), (2, 30)],
    "alt_start_codon_pos": 0,
    "alt_start_codon_exon": 1,
    "alt_last_codon": "TAA",
    "alt_valid_stop": True,
    "alt_is_premature": True,
    "start_loss": False,
    "stop_loss": False,
    "alt_transcript_length": 256,
    "alt_cds_start_in_transcript": 4,
    **NOT_SCANNED,
    "unknown_reason": None,
}
# C>T at CDS position 189: CAG becomes the PTC TAG
U_PTC_189_SEQ = "ATG" + "GCC" * 49 + "AAG" + "GCC" * 12 + "TAG" + "GCC" * 5 + "TAA"

# Layouts C: 2 exons, the whole CDS lies in exon 1, and exon 2 holds only 3' UTR. In C15, the stop codon ends at the
# end of exon 1; in C35, 20 nt of 3' UTR follow it in exon 1.
C_CDS = "ATG" + "GCC" * 49 + "CAG" + "GCC" * 3 + "TAA"
C_REF = {
    **IDS,
    "ref_cds_start": per_strand(14, 42),
    "ref_cds_stop": per_strand(179, 207),
    "ref_cds_seq": C_CDS,
    "ref_cds_len": 165,
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "ref_cds_info": [(1, 165)],
    "cds_in_transcript": True,
    "ref_start_codon_pos": 0,
    "ref_start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_first_stop_pos": 162,
    "ref_num_stop_codons": 1,
    "ref_all_stop_codons": [(162, "TAA")],
    "ref_stop_codon_exons": [1],
    "ref_is_premature": False,
    "transcript_start": 10,
    "transcript_end": 211,
    "transcript_seq": "GACC" + C_CDS + "GCC" * 4,
    "transcript_length": 181,
    "cds_start_in_transcript": 4,
    "cds_end_in_transcript": 169,
    "transcript_exon_info": [(1, 169), (2, 12)],
    "utr3_length": 12,
    "utr5_length": 4,
    "total_exon_count": 2,
    "likely_misannotated": False,
}
C15 = Layout(Transcript(("gacc" + C_CDS, "gcc" * 4)), C_REF)
C35 = Layout(
    Transcript(("gacc" + C_CDS + "gcc" * 6 + "gc", "gcc" * 4)),
    {
        **C_REF,
        "ref_cds_start": per_strand(14, 62),
        "ref_cds_stop": per_strand(179, 227),
        "transcript_end": 231,
        "transcript_seq": "GACC" + C_CDS + "GCC" * 6 + "GC" + "GCC" * 4,
        "transcript_length": 201,
        "transcript_exon_info": [(1, 189), (2, 12)],
        "utr3_length": 32,
    },
)
C_PTC = {
    "variant_id": "var1",
    "ref": per_strand("C", "G"),
    "alt": per_strand("T", "A"),
    "alt_cds_seq": "ATG" + "GCC" * 49 + "TAG" + "GCC" * 3 + "TAA",
    "alt_cds_len": 165,
    "alt_cds_info": [(1, 165)],
    "alt_start_codon_pos": 0,
    "alt_start_codon_exon": 1,
    "alt_last_codon": "TAA",
    "alt_valid_stop": True,
    "alt_first_stop_codon": "TAG",
    "alt_first_stop_pos": 150,
    "alt_num_stop_codons": 2,
    "alt_all_stop_codons": [(150, "TAG"), (162, "TAA")],
    "alt_stop_codon_exons": [1, 1],
    "alt_is_premature": True,
    "start_loss": False,
    "stop_loss": False,
    "alt_cds_start_in_transcript": 4,
    **NOT_SCANNED,
    "unknown_reason": None,
    "upstream_exon_count": 0,
    "downstream_exon_count": 1,
    "ptc_to_start_codon": 150,
    "ptc_less_than_150nt_to_start": False,
    "stop_codon_distance": 12,
    **NO_RULE,
    "nmd_50nt_penultimate_rule": True,
    "nmd_escape": True,
}

# Layouts L: 2 exons, the PTC exon 1 has 407 nt (L407) or 408 nt (L408). Its CDS part has 403 nt in both, and only
# its 5' UTR differs.
L_CDS = "ATG" + "GCC" * 49 + "CAG" + "GCC" * 89 + "TAA"
L_EXON_1 = "ATG" + "GCC" * 49 + "CAG" + "GCC" * 83 + "G"
L_EXON_2 = "CC" + "GCC" * 5 + "TAA" + "gccgcc"
L_REF = {
    **IDS,
    "ref_cds_start": per_strand(14, 16),
    "ref_cds_stop": per_strand(457, 459),
    "ref_cds_seq": L_CDS,
    "ref_cds_len": 423,
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "ref_cds_info": [(1, 403), (2, 20)],
    "cds_in_transcript": True,
    "ref_start_codon_pos": 0,
    "ref_start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_first_stop_pos": 420,
    "ref_num_stop_codons": 1,
    "ref_all_stop_codons": [(420, "TAA")],
    "ref_stop_codon_exons": [2],
    "ref_is_premature": False,
    "transcript_start": 10,
    "transcript_end": 463,
    "transcript_seq": "GACC" + L_CDS + "GCCGCC",
    "transcript_length": 433,
    "cds_start_in_transcript": 4,
    "cds_end_in_transcript": 427,
    "transcript_exon_info": [(1, 407), (2, 26)],
    "utr3_length": 6,
    "utr5_length": 4,
    "total_exon_count": 2,
    "likely_misannotated": False,
}
L407 = Layout(Transcript(("gacc" + L_EXON_1, L_EXON_2)), L_REF)
L408 = Layout(
    Transcript(("cgacc" + L_EXON_1, L_EXON_2)),
    {
        **L_REF,
        "ref_cds_start": per_strand(15, 16),
        "ref_cds_stop": per_strand(458, 459),
        "transcript_end": 464,
        "transcript_seq": "CGACC" + L_CDS + "GCCGCC",
        "transcript_length": 434,
        "cds_start_in_transcript": 5,
        "cds_end_in_transcript": 428,
        "transcript_exon_info": [(1, 408), (2, 26)],
        "utr5_length": 5,
    },
)
L_PTC_SEQ = "ATG" + "GCC" * 49 + "TAG" + "GCC" * 89 + "TAA"
L_PTC = {
    "variant_id": "var1",
    "ref": per_strand("C", "G"),
    "alt": per_strand("T", "A"),
    "alt_cds_seq": L_PTC_SEQ,
    "alt_cds_len": 423,
    "alt_cds_info": [(1, 403), (2, 20)],
    "alt_start_codon_pos": 0,
    "alt_start_codon_exon": 1,
    "alt_last_codon": "TAA",
    "alt_valid_stop": True,
    "alt_first_stop_codon": "TAG",
    "alt_first_stop_pos": 150,
    "alt_num_stop_codons": 2,
    "alt_all_stop_codons": [(150, "TAG"), (420, "TAA")],
    "alt_stop_codon_exons": [1, 2],
    "alt_is_premature": True,
    "start_loss": False,
    "stop_loss": False,
    **NOT_SCANNED,
    "unknown_reason": None,
    "upstream_exon_count": 0,
    "downstream_exon_count": 1,
    "ptc_to_start_codon": 150,
    "ptc_less_than_150nt_to_start": False,
    "stop_codon_distance": 270,
    "ptc_to_intron": 253,
}

# Layouts M: 2 exons, exon 2 is the last exon and has 432 nt
M_CDS = "ATG" + "GCC" * 30 + "CAG" + "GCC" * 68 + "TGG" + "GCC" * 61 + "TAA"
M_EXONS = ("gacc" + "ATG" + "GCC" * 20, "GCC" * 10 + "CAG" + "GCC" * 68 + "TGG" + "GCC" * 61 + "TAA" + "gccgcc")
M_REF = {
    **IDS,
    "ref_cds_start": per_strand(14, 16),
    "ref_cds_stop": per_strand(523, 525),
    "ref_cds_seq": M_CDS,
    "ref_cds_len": 489,
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "ref_cds_info": [(1, 63), (2, 426)],
    "cds_in_transcript": True,
    "ref_start_codon_pos": 0,
    "ref_start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_first_stop_pos": 486,
    "ref_num_stop_codons": 1,
    "ref_all_stop_codons": [(486, "TAA")],
    "ref_stop_codon_exons": [2],
    "ref_is_premature": False,
    "transcript_start": 10,
    "transcript_end": 529,
    "transcript_seq": "GACC" + M_CDS + "GCCGCC",
    "transcript_length": 499,
    "cds_start_in_transcript": 4,
    "cds_end_in_transcript": 493,
    "transcript_exon_info": [(1, 67), (2, 432)],
    "utr3_length": 6,
    "utr5_length": 4,
    "total_exon_count": 2,
    "likely_misannotated": False,
}
M = Layout(Transcript(M_EXONS), M_REF)


def number_the_cds_of_exon_2_as_exon_9(rows, strand):
    """Give the CDS and stop_codon rows of exon 2 the exon_number 9, which no exon row has."""
    return [
        [*row[:8], row[8].replace("exon_number=2", "exon_number=9")] if row[2] in ("CDS", "stop_codon") else row
        for row in rows
    ]


# Layout M9: layout M, but the CDS rows of exon 2 have the exon_number 9
M9 = Layout(
    Transcript(M_EXONS, edit_gff3=number_the_cds_of_exon_2_as_exon_9),
    {**M_REF, "ref_cds_info": [(1, 63), (9, 426)], "ref_stop_codon_exons": [9]},
)
M_SNV = {
    "variant_id": "var1",
    "ref": per_strand("C", "G"),
    "alt": per_strand("T", "A"),
    "alt_cds_start": per_strand(14, 16),
    "alt_cds_stop": per_strand(523, 525),
    "alt_cds_len": 489,
    "alt_start_codon_pos": 0,
    "alt_start_codon_exon": 1,
    "alt_last_codon": "TAA",
    "alt_valid_stop": True,
    "alt_num_stop_codons": 2,
    "alt_is_premature": True,
    "start_loss": False,
    "stop_loss": False,
    "alt_transcript_length": 499,
    "alt_cds_start_in_transcript": 4,
    **NOT_SCANNED,
    "unknown_reason": None,
}

# Layout SL: 3 exons. CAT GCC holds an ATG in the frame +1 at CDS position 100, and GTG AGC a TGA in the same frame
# at CDS position 160.
SL_CDS = "ATG" + "GCC" * 32 + "CAT" + "GCC" * 19 + "GTGAGC" + "GCC" * 21 + "TAA"
SL = Layout(
    Transcript(
        (
            "gacc" + "ATG" + "GCC" * 20,
            "GCC" * 12 + "CAT" + "GCC" * 19 + "GTGAGC" + "GCC" * 12,
            "GCC" * 9 + "TAA" + "gccgcc",
        )
    ),
    {
        **IDS,
        "ref_cds_start": per_strand(14, 16),
        "ref_cds_stop": per_strand(285, 287),
        "ref_cds_seq": SL_CDS,
        "ref_cds_len": 231,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 63), (2, 138), (3, 30)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 228,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(228, "TAA")],
        "ref_stop_codon_exons": [3],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 291,
        "transcript_seq": "GACC" + SL_CDS + "GCCGCC",
        "transcript_length": 241,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 235,
        "transcript_exon_info": [(1, 67), (2, 138), (3, 36)],
        "utr3_length": 6,
        "utr5_length": 4,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# Layout NF: 3 exons, a cds_start_NF transcript without a start codon and without a 5' UTR
NF_CDS = "GCC" * 3 + "CAG" + "GCC" * 52 + "TAA"
NF = Layout(
    Transcript(
        ("GCC" * 3 + "CAG" + "GCC" * 17, "GCC" * 30, "GCC" * 5 + "TAA" + "gccgcc"),
        start_codon=False,
        tags=("cds_start_NF",),
    ),
    {
        **IDS,
        "ref_cds_start": per_strand(10, 16),
        "ref_cds_stop": per_strand(221, 227),
        "ref_cds_seq": NF_CDS,
        "ref_cds_len": 171,
        "has_start_codon": False,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 63), (2, 90), (3, 18)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": None,
        "ref_start_codon_exon": None,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 168,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(168, "TAA")],
        "ref_stop_codon_exons": [3],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 227,
        "transcript_seq": NF_CDS + "GCCGCC",
        "transcript_length": 177,
        "cds_start_in_transcript": 0,
        "cds_end_in_transcript": 171,
        "transcript_exon_info": [(1, 63), (2, 90), (3, 24)],
        "utr3_length": 6,
        "utr5_length": 0,
        "total_exon_count": 3,
        "likely_misannotated": True,
    },
)

CASES = [
    # NR-04 (51 nt, False), NR-12 and PF-12 (150 nt, False), NR-18, NR-03 (PTC in an earlier exon), NR-16
    Case(
        "ptc_51_nt_before_the_last_exon_junction_and_150_nt_after_the_start_codon_escapes_by_no_rule",
        """
        CDS     -4   0                    63                 147 150                          201                258       267
        ref 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..75.. GCC CAG TTG AGC GCC ..36.. GCC GCC]|[GCC GCC ..48.. GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..75.. GCC CAG TAG AGC GCC ..36.. GCC GCC]|[GCC GCC ..48.. GCC TAA gccgcc] 3'
                                                                  ^ T>A
                                                                 *** PTC
                                                                 <-- ptc_to_intron = 51 -->
                     <------- ptc_to_start_codon = 150 -------->
        exon 1: 67 nt, exon 2: 138 nt, exon 3: 66 nt
        """,
        A201,
        A_PTC_150,
        {
            **A_SNV,
            "ref": per_strand("T", "A"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(185, 145),
            "end_variant": per_strand(186, 146),
            "alt_cds_seq": A_PTC_150_SEQ,
            "alt_cds_info": [(1, 63), (2, 138), (3, 60)],
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 150,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(150, "TAG"), (258, "TAA")],
            "alt_stop_codon_exons": [2, 3],
            "alt_is_premature": True,
            "alt_transcript_seq": "GACC" + A_PTC_150_SEQ + "GCCGCC",
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 150,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 138,
            "stop_codon_distance": 108,
            "ptc_to_intron": 51,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Mark("alt", 154, 157, "*", "PTC"),
            Span("alt", 154, 205, "ptc_to_intron = 51"),
            Span("alt", 4, 154, "ptc_to_start_codon = 150"),
        ),
        ruler=Ruler((-4, 0, 63, 147, 150, 201, 258, 267), "CDS"),
    ),
    # NR-04 (50 nt, True)
    Case(
        "ptc_50_nt_before_the_last_exon_junction_escapes_by_the_50nt_rule",
        """
        CDS     -4   0                    63                 147 150                         200              258       267
        ref 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..75.. GCC CAG TTG AGC GCC ..36.. GCC GC]|[C GCC ..51.. GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..75.. GCC CAG TAG AGC GCC ..36.. GCC GC]|[C GCC ..51.. GCC TAA gccgcc] 3'
                                                                  ^ T>A
                                                                 *** PTC
                                                                 <-----------------------> ptc_to_intron = 50
        exon 1: 67 nt, exon 2: 137 nt, exon 3: 67 nt
        """,
        A200,
        A_PTC_150,
        {
            **A_SNV,
            "ref": per_strand("T", "A"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(185, 145),
            "end_variant": per_strand(186, 146),
            "alt_cds_seq": A_PTC_150_SEQ,
            "alt_cds_info": [(1, 63), (2, 137), (3, 61)],
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 150,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(150, "TAG"), (258, "TAA")],
            "alt_stop_codon_exons": [2, 3],
            "alt_is_premature": True,
            "alt_transcript_seq": "GACC" + A_PTC_150_SEQ + "GCCGCC",
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 150,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 137,
            "stop_codon_distance": 108,
            "ptc_to_intron": 50,
            **NO_RULE,
            "nmd_50nt_penultimate_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Mark("alt", 154, 157, "*", "PTC"), Span("alt", 154, 204, "ptc_to_intron = 50")),
        ruler=Ruler((-4, 0, 63, 147, 150, 200, 258, 267), "CDS"),
    ),
    # NR-05 (0 nt, False): the PTC starts at the first base of the last exon, so the last exon rule fires instead
    Case(
        "ptc_at_the_first_base_of_the_last_exon_is_0_nt_before_the_junction_and_escapes_by_the_last_exon_rule",
        """
        CDS     -4   0                    63                 147   150                    258       267
        ref 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..75.. GCC CAG]|[TTG AGC GCC ..96.. GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..75.. GCC CAG]|[TAG AGC GCC ..96.. GCC TAA gccgcc] 3'
                                                                    ^ T>A
                                                                   *** PTC
                                                                   <----- ptc_to_intron = 117 ----->
        exon 1: 67 nt, exon 2: 87 nt, exon 3: 117 nt
        The PTC TAG is the first codon of exon 3, the last exon, so ptc_to_intron counts to the transcript end.
        """,
        A150,
        A_PTC_150,
        {
            **A_SNV,
            "ref": per_strand("T", "A"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(205, 125),
            "end_variant": per_strand(206, 126),
            "alt_cds_seq": A_PTC_150_SEQ,
            "alt_cds_info": [(1, 63), (2, 87), (3, 111)],
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 150,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(150, "TAG"), (258, "TAA")],
            "alt_stop_codon_exons": [3, 3],
            "alt_is_premature": True,
            "alt_transcript_seq": "GACC" + A_PTC_150_SEQ + "GCCGCC",
            "upstream_exon_count": 2,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 150,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 117,
            "stop_codon_distance": 108,
            "ptc_to_intron": 117,
            **NO_RULE,
            "nmd_last_exon_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Mark("alt", 154, 157, "*", "PTC"), Span("alt", 154, 271, "ptc_to_intron = 117")),
        ruler=Ruler((-4, 0, 63, 147, 150, 258, 267), "CDS"),
    ),
    # NR-05 (1 nt, True): the PTC codon is split over the last exon junction, its first base lies in exon 2
    Case(
        "ptc_codon_split_1_nt_before_the_last_exon_junction_escapes_by_the_50nt_rule",
        """
        CDS     -4   0                    63                 147 150 151                   258       267
        ref 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..75.. GCC CAG T]|[TG AGC GCC ..96.. GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..75.. GCC CAG T]|[AG AGC GCC ..96.. GCC TAA gccgcc] 3'
                                                                     ^ T>A
                                                                 * PTC TAG
                                                                 <> ptc_to_intron = 1
        exon 1: 67 nt, exon 2: 88 nt, exon 3: 116 nt
        The PTC TAG is split T|AG over the last exon junction. It lies in exon 2, the exon of its first base, 1 nt
        before the junction.
        """,
        A151,
        Change("TTTCAG[T>A]GAGC"),
        {
            **A_SNV,
            "ref": per_strand("T", "A"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(205, 125),
            "end_variant": per_strand(206, 126),
            "alt_cds_seq": A_PTC_150_SEQ,
            "alt_cds_info": [(1, 63), (2, 88), (3, 110)],
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 150,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(150, "TAG"), (258, "TAA")],
            "alt_stop_codon_exons": [2, 3],
            "alt_is_premature": True,
            "alt_transcript_seq": "GACC" + A_PTC_150_SEQ + "GCCGCC",
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 150,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 88,
            "stop_codon_distance": 108,
            "ptc_to_intron": 1,
            **NO_RULE,
            "nmd_50nt_penultimate_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Mark("alt", 154, 155, "*", "PTC TAG"), Span("alt", 154, 155, "ptc_to_intron = 1")),
        ruler=Ruler((-4, 0, 63, 147, 150, 151, 258, 267), "CDS"),
    ),
    # NR-12 and PF-12 (147 nt, True), NR-17 (nmd_escape by the start-proximal rule only)
    Case(
        "ptc_147_nt_after_the_start_codon_escapes_by_the_start_proximal_rule_only",
        """
        CDS     -4   0                    63                     147                          201                258       267
        ref 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..72.. GCC GCC CAG TTG AGC ..39.. GCC GCC]|[GCC GCC ..48.. GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..72.. GCC GCC TAG TTG AGC ..39.. GCC GCC]|[GCC GCC ..48.. GCC TAA gccgcc] 3'
                                                                 ^ C>T
                                                                 *** PTC
                     <------- ptc_to_start_codon = 147 -------->
                                                                 <-- ptc_to_intron = 54 -->
        exon 1: 67 nt, exon 2: 138 nt, exon 3: 66 nt
        """,
        A201,
        Change("GCC[C>T]AGTTG"),
        {
            **A_SNV,
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "start_variant": per_strand(181, 149),
            "end_variant": per_strand(182, 150),
            "alt_cds_seq": A_HEAD + A_MID + "TAGTTGAGC" + A_TAIL,
            "alt_cds_info": [(1, 63), (2, 138), (3, 60)],
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 147,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(147, "TAG"), (258, "TAA")],
            "alt_stop_codon_exons": [2, 3],
            "alt_is_premature": True,
            "alt_transcript_seq": "GACC" + A_HEAD + A_MID + "TAGTTGAGC" + A_TAIL + "GCCGCC",
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 147,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 138,
            "stop_codon_distance": 111,
            "ptc_to_intron": 54,
            **NO_RULE,
            "nmd_start_proximal_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Mark("alt", 151, 154, "*", "PTC"),
            Span("alt", 4, 151, "ptc_to_start_codon = 147"),
            Span("alt", 151, 205, "ptc_to_intron = 54"),
        ),
        ruler=Ruler((-4, 0, 63, 147, 201, 258, 267), "CDS"),
    ),
    # NR-01
    Case(
        "nonsense_snv_in_the_last_exon_escapes_by_the_last_exon_rule",
        """
        CDS     -4   0                    63                    201                    225                    258       267
        ref 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..129.. GCC]|[GCC GCC ..12.. GCC GCC TGG GCC GCC ..21.. GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..54.. GAG]|[GCC GCC ..129.. GCC]|[GCC GCC ..12.. GCC GCC TGA GCC GCC ..21.. GCC TAA gccgcc] 3'
                                                                                         ^ G>A
                                                                                       *** PTC
                                                                                       <----- ptc_to_intron = 42 ------>
        exon 1: 67 nt, exon 2: 138 nt, exon 3: 66 nt
        """,
        A201,
        Change("GCCTG[G>A]GCC"),
        {
            **A_SNV,
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(281, 49),
            "end_variant": per_strand(282, 50),
            "alt_cds_seq": A_HEAD + A_MID + "CAGTTGAGC" + "GCC" * 23 + "TGA" + "GCC" * 10 + "TAA",
            "alt_cds_info": [(1, 63), (2, 138), (3, 60)],
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 225,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(225, "TGA"), (258, "TAA")],
            "alt_stop_codon_exons": [3, 3],
            "alt_is_premature": True,
            "alt_transcript_seq": "GACC" + A_HEAD + A_MID + "CAGTTGAGC" + "GCC" * 23 + "TGA" + "GCC" * 10 + "TAAGCCGCC",
            "upstream_exon_count": 2,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 225,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 66,
            "stop_codon_distance": 33,
            "ptc_to_intron": 42,
            **NO_RULE,
            "nmd_last_exon_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Mark("alt", 229, 232, "*", "PTC"), Span("alt", 229, 271, "ptc_to_intron = 42")),
        ruler=Ruler((-4, 0, 63, 201, 225, 258, 267), "CDS"),
    ),
    # NR-08: the deletion moves the last exon junction from 201 to 200 in alt CDS coordinates. Measured to the
    # junction in ref CDS coordinates, the PTC would lie 51 nt before it.
    Case(
        "frameshift_deletion_upstream_of_the_ptc_moves_the_last_exon_junction_into_the_50nt_rule",
        """
        CDS         -4   0                      30                       63                                          201                258       267
        ref     5' [gacc ATG GCC ..18.. GCC GCC AAA GCC GCC ..21.. GAG]|[GCC GCC ..78.. CAG TTG AGC ..39.. GCC GCC]|[GCC GCC ..48.. GCC TAA gccgcc] 3'
        alt     5' [gacc ATG GCC ..18.. GCC GCC -AA GCC GCC ..21.. GAG]|[GCC GCC ..78.. CAG TTG AGC ..39.. GCC GCC]|[GCC GCC ..48.. GCC TAA gccgcc] 3'
                                                ^ A>-
        alt CDS                                                          62                  150                     200                257       266
                                                                                             * PTC TGA, in the alt frame
                                                                                             <-------------------> ptc_to_intron = 50
        exon 1: 67 nt, exon 2: 138 nt, exon 3: 66 nt
        A>- deletes one A of AAA at CDS 30 to 32 and moves the frame by 1. The alt frame reads the PTC TGA at alt
        CDS 150. The deletion moves the last exon junction from CDS 201 to alt CDS 200.
        """,
        A201,
        Change("GCC[A>]AAGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CA", "TT"),
            "alt": per_strand("C", "T"),
            "start_variant": per_strand(43, 285),
            "end_variant": per_strand(45, 287),
            "alt_cds_start": per_strand(14, 16),
            "alt_cds_stop": per_strand(315, 317),
            "alt_cds_seq": "ATG" + "GCC" * 9 + "AA" + "GCC" * 9 + "GAG" + A_MID + "CAGTTGAGC" + A_TAIL,
            "alt_cds_len": 260,
            "alt_cds_info": [(1, 62), (2, 138), (3, 60)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 150,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(150, "TGA")],
            "alt_stop_codon_exons": [2],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATG"
            + "GCC" * 9
            + "AA"
            + "GCC" * 9
            + "GAG"
            + A_MID
            + "CAGTTGAGC"
            + A_TAIL
            + "GCCGCC",
            "alt_transcript_length": 270,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 150,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 138,
            "stop_codon_distance": 107,
            "ptc_to_intron": 50,
            **NO_RULE,
            "nmd_50nt_penultimate_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": [(1, 66), (2, 138), (3, 66)],
            "nmd_model_status": "ok",
        },
        equivalent=(Change("GCCAA[A>]GCC"), Change("GCC[AAA>AA]GCC")),
        marks=(
            Ruler((62, 150, 200, 257, 266), "CDS", "alt"),
            Mark("alt", 154, 155, "*", "PTC TGA, in the alt frame"),
            Span("alt", 154, 204, "ptc_to_intron = 50"),
        ),
        ruler=Ruler((-4, 0, 30, 63, 201, 258, 267), "CDS"),
    ),
    # NR-21
    Case(
        "snv_in_the_splice_donor_gives_an_unknown_row_with_null_rules",
        """
        CDS     -4   0                                           63
        ref 5' [gacc ATG GCC ..51.. GCC GAG]gtaagtcccccccctttcag[GCC GCC ..129.. GCC]|[GCC ..54.. TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..51.. GCC GAG]ataagtcccccccctttcag[GCC GCC ..129.. GCC]|[GCC ..54.. TAA gccgcc] 3'
                                            ^ G>A
        exon 1: 67 nt, exon 2: 138 nt, exon 3: 66 nt
        G>A at the first base of intron 1 turns the donor GT into AT.
        """,
        A201,
        Change("GAG[G>A]TAAGT"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(77, 253),
            "end_variant": per_strand(78, 254),
            **UNKNOWN_ALT,
            "unknown_reason": "splice_site_destroyed",
        },
        ruler=Ruler((-4, 0, 63), "CDS"),
    ),
    # NR-02, NR-09 (no last exon junction: 50 nt rule False), NR-15
    Case(
        "ptc_in_a_single_exon_transcript_escapes_by_the_single_exon_and_last_exon_rules",
        """
        CDS     -4   0                   147 150                    258       267
        ref 5' [gacc ATG GCC ..138.. GCC CAG TTG AGC GCC ..96.. GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..138.. GCC CAG TAG AGC GCC ..96.. GCC TAA gccgcc] 3'
                                              ^ T>A
                                             *** PTC
                                             <----- ptc_to_intron = 117 ----->
        exon 1: 271 nt, the only exon
        """,
        SINGLE,
        A_PTC_150,
        {
            **SINGLE_SNV,
            "ref": per_strand("T", "A"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(165, 125),
            "end_variant": per_strand(166, 126),
            "alt_cds_seq": A_PTC_150_SEQ,
            "alt_cds_info": [(1, 261)],
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 150,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(150, "TAG"), (258, "TAA")],
            "alt_stop_codon_exons": [1, 1],
            "alt_is_premature": True,
            "alt_transcript_seq": "GACC" + A_PTC_150_SEQ + "GCCGCC",
            "upstream_exon_count": 0,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 150,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 271,
            "stop_codon_distance": 108,
            "ptc_to_intron": 117,
            **NO_RULE,
            "nmd_last_exon_rule": True,
            "nmd_single_exon_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Mark("alt", 154, 157, "*", "PTC"), Span("alt", 154, 271, "ptc_to_intron = 117")),
        ruler=Ruler((-4, 0, 147, 150, 258, 267), "CDS"),
    ),
    # NR-20: a row that is not a PTC row has every rule False, also the single exon rule of a single exon transcript
    Case(
        "missense_snv_in_a_single_exon_transcript_is_no_ptc_row_and_every_rule_is_false",
        """
        CDS     -4   0                       147                    258       267
        ref 5' [gacc ATG GCC ..135.. GCC GCC CAG TTG AGC ..99.. GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..135.. GCC GCC GAG TTG AGC ..99.. GCC TAA gccgcc] 3'
                                             ^ C>G
        exon 1: 271 nt, the only exon
        """,
        SINGLE,
        Change("GCC[C>G]AGTTG"),
        {
            **SINGLE_SNV,
            "ref": per_strand("C", "G"),
            "alt": per_strand("G", "C"),
            "start_variant": per_strand(161, 129),
            "end_variant": per_strand(162, 130),
            "alt_cds_seq": A_HEAD + A_MID + "GAGTTGAGC" + A_TAIL,
            "alt_cds_info": [(1, 261)],
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 258,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(258, "TAA")],
            "alt_stop_codon_exons": [1],
            "alt_is_premature": False,
            "alt_transcript_seq": "GACC" + A_HEAD + A_MID + "GAGTTGAGC" + A_TAIL + "GCCGCC",
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((-4, 0, 147, 258, 267), "CDS"),
    ),
    # NR-06 (90 nt, False), NR-03: the PTC lies 30 nt before the end of its exon, but the last exon junction lies in
    # the 3' UTR, 90 nt after the PTC
    Case(
        "ptc_30_nt_before_its_exon_end_but_90_nt_before_the_last_exon_junction_in_the_3utr_escapes_by_no_rule",
        """
        CDS     -4   0                       150                          180                207               240         252
        ref 5' [gacc ATG GCC ..138.. GCC GCC AAG GCC GCC ..15.. GCC GCC]|[GCC GCC ..18.. GCC TAA g..25..cgcc]|[gccgccgccgcc] 3'
        alt 5' [gacc ATG GCC ..138.. GCC GCC TAG GCC GCC ..15.. GCC GCC]|[GCC GCC ..18.. GCC TAA g..25..cgcc]|[gccgccgccgcc] 3'
                                             ^ A>T
                                             *** PTC
                                             <-- ptc_to_intron = 30 -->
                                             <-------------- 90 nt to the last exon junction -------------->
        exon 1: 184 nt, exon 2: 60 nt, exon 3: 12 nt
        """,
        U51,
        Change("GCC[A>T]AGGCC"),
        {
            **U_SNV,
            "ref": per_strand("A", "T"),
            "alt": per_strand("T", "A"),
            "start_variant": per_strand(164, 151),
            "end_variant": per_strand(165, 152),
            "alt_cds_seq": "ATG" + "GCC" * 49 + "TAG" + "GCC" * 12 + "CAG" + "GCC" * 5 + "TAA",
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 150,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(150, "TAG"), (207, "TAA")],
            "alt_stop_codon_exons": [1, 2],
            "alt_transcript_seq": "GACC"
            + "ATG"
            + "GCC" * 49
            + "TAG"
            + "GCC" * 12
            + "CAG"
            + "GCC" * 5
            + "TAA"
            + "GCC" * 14,
            "upstream_exon_count": 0,
            "downstream_exon_count": 2,
            "ptc_to_start_codon": 150,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 184,
            "stop_codon_distance": 57,
            "ptc_to_intron": 30,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Mark("alt", 154, 157, "*", "PTC"),
            Span("alt", 154, 184, "ptc_to_intron = 30"),
            Span("alt", 154, 244, "90 nt to the last exon junction"),
        ),
        ruler=Ruler((-4, 0, 150, 180, 207, 240, 252), "CDS"),
    ),
    # NR-06 (51 nt, False), NR-03 (PTC in the last CDS exon before an exon with only 3' UTR)
    Case(
        "ptc_in_the_last_cds_exon_51_nt_before_the_last_exon_junction_in_the_3utr_escapes_by_no_rule",
        """
        CDS     -4   0                     180         189                     207               240         252
        ref 5' [gacc ATG GCC ..171.. GCC]|[GCC GCC GCC CAG GCC GCC GCC GCC GCC TAA g..25..cgcc]|[gccgccgccgcc] 3'
        alt 5' [gacc ATG GCC ..171.. GCC]|[GCC GCC GCC TAG GCC GCC GCC GCC GCC TAA g..25..cgcc]|[gccgccgccgcc] 3'
                                                       ^ C>T
                                                       *** PTC
                                                       <-------- ptc_to_intron = 51 --------->
        exon 1: 184 nt, exon 2: 60 nt, exon 3: 12 nt
        """,
        U51,
        Change("GCC[C>T]AGGCC"),
        {
            **U_SNV,
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "start_variant": per_strand(223, 92),
            "end_variant": per_strand(224, 93),
            "alt_cds_seq": U_PTC_189_SEQ,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 189,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(189, "TAG"), (207, "TAA")],
            "alt_stop_codon_exons": [2, 2],
            "alt_transcript_seq": "GACC" + U_PTC_189_SEQ + "GCC" * 14,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 189,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 60,
            "stop_codon_distance": 18,
            "ptc_to_intron": 51,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Mark("alt", 193, 196, "*", "PTC"), Span("alt", 193, 244, "ptc_to_intron = 51")),
        ruler=Ruler((-4, 0, 180, 189, 207, 240, 252), "CDS"),
    ),
    # NR-06 (50 nt, True)
    Case(
        "ptc_in_the_last_cds_exon_50_nt_before_the_last_exon_junction_in_the_3utr_escapes_by_the_50nt_rule",
        """
        CDS     -4   0                     180         189                     207               239         251
        ref 5' [gacc ATG GCC ..171.. GCC]|[GCC GCC GCC CAG GCC GCC GCC GCC GCC TAA g..24..ccgc]|[gccgccgccgcc] 3'
        alt 5' [gacc ATG GCC ..171.. GCC]|[GCC GCC GCC TAG GCC GCC GCC GCC GCC TAA g..24..ccgc]|[gccgccgccgcc] 3'
                                                       ^ C>T
                                                       *** PTC
                                                       <-------- ptc_to_intron = 50 --------->
        exon 1: 184 nt, exon 2: 59 nt, exon 3: 12 nt
        """,
        U50,
        Change("GCC[C>T]AGGCC"),
        {
            **U_SNV,
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "start_variant": per_strand(223, 91),
            "end_variant": per_strand(224, 92),
            "alt_cds_start": per_strand(14, 71),
            "alt_cds_stop": per_strand(244, 301),
            "alt_cds_seq": U_PTC_189_SEQ,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 189,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(189, "TAG"), (207, "TAA")],
            "alt_stop_codon_exons": [2, 2],
            "alt_transcript_seq": "GACC" + U_PTC_189_SEQ + "GCC" * 9 + "GC" + "GCC" * 4,
            "alt_transcript_length": 255,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 189,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 59,
            "stop_codon_distance": 18,
            "ptc_to_intron": 50,
            **NO_RULE,
            "nmd_50nt_penultimate_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Mark("alt", 193, 196, "*", "PTC"), Span("alt", 193, 243, "ptc_to_intron = 50")),
        ruler=Ruler((-4, 0, 180, 189, 207, 239, 251), "CDS"),
    ),
    # NR-07: the whole CDS in exon 1 of 2, followed by 3' UTR in exon 1
    Case(
        "whole_cds_in_exon_1_of_2_ptc_35_nt_before_the_junction_in_the_3utr_escapes_by_the_50nt_rule",
        """
        CDS     -4   0                       150             162               185         197
        ref 5' [gacc ATG GCC ..138.. GCC GCC CAG GCC GCC GCC TAA g..15..ccgc]|[gccgccgccgcc] 3'
        alt 5' [gacc ATG GCC ..138.. GCC GCC TAG GCC GCC GCC TAA g..15..ccgc]|[gccgccgccgcc] 3'
                                             ^ C>T
                                             *** PTC
                                             <---- ptc_to_intron = 35 ----->
        exon 1: 189 nt, exon 2: 12 nt
        """,
        C35,
        Change("GCC[C>T]AGGCC"),
        {
            **C_PTC,
            "start_variant": per_strand(164, 76),
            "end_variant": per_strand(165, 77),
            "alt_cds_start": per_strand(14, 62),
            "alt_cds_stop": per_strand(179, 227),
            "alt_transcript_seq": "GACC"
            + "ATG"
            + "GCC" * 49
            + "TAG"
            + "GCC" * 3
            + "TAA"
            + "GCC" * 6
            + "GC"
            + "GCC" * 4,
            "alt_transcript_length": 201,
            "ptc_exon_length": 189,
            "ptc_to_intron": 35,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Mark("alt", 154, 157, "*", "PTC"), Span("alt", 154, 189, "ptc_to_intron = 35")),
        ruler=Ruler((-4, 0, 150, 162, 185, 197), "CDS"),
    ),
    # NR-07: the whole CDS in exon 1 of 2, the stop codon ends at the end of exon 1
    Case(
        "whole_cds_ending_at_the_end_of_exon_1_of_2_ptc_15_nt_before_the_junction_escapes_by_the_50nt_rule",
        """
        CDS     -4   0                       150             162   165         177
        ref 5' [gacc ATG GCC ..138.. GCC GCC CAG GCC GCC GCC TAA]|[gccgccgccgcc] 3'
        alt 5' [gacc ATG GCC ..138.. GCC GCC TAG GCC GCC GCC TAA]|[gccgccgccgcc] 3'
                                             ^ C>T
                                             *** PTC
                                             <-----------------> ptc_to_intron = 15
        exon 1: 169 nt, exon 2: 12 nt
        """,
        C15,
        Change("GCC[C>T]AGGCC"),
        {
            **C_PTC,
            "start_variant": per_strand(164, 56),
            "end_variant": per_strand(165, 57),
            "alt_cds_start": per_strand(14, 42),
            "alt_cds_stop": per_strand(179, 207),
            "alt_transcript_seq": "GACC" + "ATG" + "GCC" * 49 + "TAG" + "GCC" * 3 + "TAA" + "GCC" * 4,
            "alt_transcript_length": 181,
            "ptc_exon_length": 169,
            "ptc_to_intron": 15,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Mark("alt", 154, 157, "*", "PTC"), Span("alt", 154, 169, "ptc_to_intron = 15")),
        ruler=Ruler((-4, 0, 150, 162, 165, 177), "CDS"),
    ),
    # NR-11 (407 nt, False)
    Case(
        "ptc_in_an_exon_of_407_nt_escapes_by_no_rule",
        """
        CDS     -4   0                       150                         403                    420       429
        ref 5' [gacc ATG GCC ..138.. GCC GCC CAG GCC GCC ..240.. GCC G]|[CC GCC GCC GCC GCC GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..138.. GCC GCC TAG GCC GCC ..240.. GCC G]|[CC GCC GCC GCC GCC GCC TAA gccgcc] 3'
                                             ^ C>T
                <-------------- ptc_exon_length = 407 --------------->
        exon 1: 407 nt, exon 2: 26 nt
        """,
        L407,
        Change("GCC[C>T]AGGCC"),
        {
            **L_PTC,
            "start_variant": per_strand(164, 308),
            "end_variant": per_strand(165, 309),
            "alt_cds_start": per_strand(14, 16),
            "alt_cds_stop": per_strand(457, 459),
            "alt_transcript_seq": "GACC" + L_PTC_SEQ + "GCCGCC",
            "alt_transcript_length": 433,
            "alt_cds_start_in_transcript": 4,
            "ptc_exon_length": 407,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Span("alt", 0, 407, "ptc_exon_length = 407"),),
        ruler=Ruler((-4, 0, 150, 403, 420, 429), "CDS"),
    ),
    # NR-11 (408 nt, True): the CDS part of the exon has 403 nt, the 5' UTR counts too
    Case(
        "ptc_in_an_exon_of_408_nt_with_its_5utr_escapes_by_the_long_exon_rule",
        """
        CDS     -5    0                       150                         403                    420       429
        ref 5' [cgacc ATG GCC ..138.. GCC GCC CAG GCC GCC ..240.. GCC G]|[CC GCC GCC GCC GCC GCC TAA gccgcc] 3'
        alt 5' [cgacc ATG GCC ..138.. GCC GCC TAG GCC GCC ..240.. GCC G]|[CC GCC GCC GCC GCC GCC TAA gccgcc] 3'
                                              ^ C>T
                <--------------- ptc_exon_length = 408 --------------->
        exon 1: 408 nt, exon 2: 26 nt
        ptc_exon_length counts the 5 nt of 5' UTR in exon 1.
        """,
        L408,
        Change("GCC[C>T]AGGCC"),
        {
            **L_PTC,
            "start_variant": per_strand(165, 308),
            "end_variant": per_strand(166, 309),
            "alt_cds_start": per_strand(15, 16),
            "alt_cds_stop": per_strand(458, 459),
            "alt_transcript_seq": "CGACC" + L_PTC_SEQ + "GCCGCC",
            "alt_transcript_length": 434,
            "alt_cds_start_in_transcript": 5,
            "ptc_exon_length": 408,
            **NO_RULE,
            "nmd_long_exon_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Span("alt", 0, 408, "ptc_exon_length = 408"),),
        ruler=Ruler((-5, 0, 150, 403, 420, 429), "CDS"),
    ),
    # NR-19: the last exon, long exon and start-proximal rules at once
    Case(
        "ptc_in_a_long_last_exon_93_nt_after_the_start_codon_escapes_by_three_rules",
        """
        CDS     -4   0                    63                     93                      486       495
        ref 5' [gacc ATG GCC ..54.. GCC]|[GCC GCC ..18.. GCC GCC CAG GCC GCC ..381.. GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..54.. GCC]|[GCC GCC ..18.. GCC GCC TAG GCC GCC ..381.. GCC TAA gccgcc] 3'
                                                                 ^ C>T
                                                                 *** PTC
                     <-------- ptc_to_start_codon = 93 -------->
                                                                 <----- ptc_to_intron = 402 ------>
                                          <---------------- ptc_exon_length = 432 ---------------->
        exon 1: 67 nt, exon 2: 432 nt
        """,
        M,
        Change("GCC[C>T]AGGCC"),
        {
            **M_SNV,
            "start_variant": per_strand(127, 411),
            "end_variant": per_strand(128, 412),
            "alt_cds_seq": "ATG" + "GCC" * 30 + "TAG" + "GCC" * 68 + "TGG" + "GCC" * 61 + "TAA",
            "alt_cds_info": [(1, 63), (2, 426)],
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 93,
            "alt_all_stop_codons": [(93, "TAG"), (486, "TAA")],
            "alt_stop_codon_exons": [2, 2],
            "alt_transcript_seq": "GACC" + "ATG" + "GCC" * 30 + "TAG" + "GCC" * 68 + "TGG" + "GCC" * 61 + "TAAGCCGCC",
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 93,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 432,
            "stop_codon_distance": 393,
            "ptc_to_intron": 402,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": True,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Mark("alt", 97, 100, "*", "PTC"),
            Span("alt", 4, 97, "ptc_to_start_codon = 93"),
            Span("alt", 97, 499, "ptc_to_intron = 402"),
            Span("alt", 67, 499, "ptc_exon_length = 432"),
        ),
        ruler=Ruler((-4, 0, 63, 93, 486, 495), "CDS"),
    ),
    # NR-22: the CDS rows of exon 2 say exon_number 9, which no exon row has. The exon features find the PTC exon by
    # its position in alt_transcript_exon_info: exon 2, the last exon, of 432 nt. So the last exon rule and the long
    # exon rule fire.
    Case(
        "ptc_in_an_exon_whose_cds_rows_have_another_exon_number_fires_the_last_exon_and_long_exon_rules",
        """
        CDS     -4   0                    63                      300                     486       495
        ref 5' [gacc ATG GCC ..54.. GCC]|[GCC GCC ..225.. GCC GCC TGG GCC GCC ..174.. GCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC ..54.. GCC]|[GCC GCC ..225.. GCC GCC TGA GCC GCC ..174.. GCC TAA gccgcc] 3'
                                                                    ^ G>A
                                                                  *** PTC
                                                                  <----- ptc_to_intron = 195 ------>
                                          <---------------- ptc_exon_length = 432 ----------------->
        exon 1: 67 nt, exon 2: 432 nt
        The CDS rows of exon 2 say exon_number 9. The PTC TGA at CDS 300 lies in exon 2 of alt_transcript_exon_info,
        the last exon.
        """,
        M9,
        Change("GCCTG[G>A]GCC"),
        {
            **M_SNV,
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(336, 202),
            "end_variant": per_strand(337, 203),
            "alt_cds_seq": "ATG" + "GCC" * 30 + "CAG" + "GCC" * 68 + "TGA" + "GCC" * 61 + "TAA",
            "alt_cds_info": [(1, 63), (9, 426)],
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 300,
            "alt_all_stop_codons": [(300, "TGA"), (486, "TAA")],
            "alt_stop_codon_exons": [9, 9],
            "alt_transcript_seq": "GACC" + "ATG" + "GCC" * 30 + "CAG" + "GCC" * 68 + "TGA" + "GCC" * 61 + "TAAGCCGCC",
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 300,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 432,
            "stop_codon_distance": 186,
            "ptc_to_intron": 195,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": True,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Mark("alt", 304, 307, "*", "PTC"),
            Span("alt", 304, 499, "ptc_to_intron = 195"),
            Span("alt", 67, 499, "ptc_exon_length = 432"),
        ),
        ruler=Ruler((-4, 0, 63, 300, 486, 495), "CDS"),
    ),
    # NR-10, NR-14: after the start loss, the PTC is the TGA of the rescued ORF. Measured from the annotated start
    # codon, it would lie 160 nt downstream; the annotated stop codon in exon 3 would fire the last exon rule.
    Case(
        "start_loss_rescued_orf_ptc_escapes_by_the_50nt_and_start_proximal_rules_measured_from_the_rescued_atg",
        """
        CDS     -4   0                        63                  100                160                     201                228       237
        ref 5' [gacc ATG GCC GCC ..51.. GCC]|[GCC GCC ..27.. GCC CAT GCC ..51.. GCC GTG AGC ..30.. GCC GCC]|[GCC GCC ..18.. GCC TAA gccgcc] 3'
        alt 5' [gacc ACG GCC GCC ..51.. GCC]|[GCC GCC ..27.. GCC CAT GCC ..51.. GCC GTG AGC ..30.. GCC GCC]|[GCC GCC ..18.. GCC TAA gccgcc] 3'
                      ^ T>C
                                                                                     * PTC TGA, in the frame +1
                                                                  <-----------------> ptc_to_start_codon = 60, in the frame +1
                                                                                     <-------------------> ptc_to_intron = 41
        exon 1: 67 nt, exon 2: 138 nt, exon 3: 36 nt
        T>C turns the start codon ATG into ACG. The scan finds the ATG at CDS 100, in the frame +1, and reads on to
        the PTC TGA at CDS 160.
        """,
        SL,
        Change("GACCA[T>C]GGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "start_variant": per_strand(15, 285),
            "end_variant": per_strand(16, 286),
            "alt_cds_start": per_strand(14, 16),
            "alt_cds_stop": per_strand(285, 287),
            "alt_cds_seq": "ACG" + "GCC" * 32 + "CAT" + "GCC" * 19 + "GTGAGC" + "GCC" * 21 + "TAA",
            "alt_cds_len": 231,
            "alt_cds_info": [(1, 63), (2, 138), (3, 30)],
            "alt_start_codon_pos": None,
            "alt_start_codon_exon": None,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 228,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(228, "TAA")],
            "alt_stop_codon_exons": [3],
            "alt_is_premature": True,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GACC"
            + "ACG"
            + "GCC" * 32
            + "CAT"
            + "GCC" * 19
            + "GTGAGC"
            + "GCC" * 21
            + "TAAGCCGCC",
            "alt_transcript_length": 241,
            "alt_cds_start_in_transcript": 4,
            "transcript_start_codon_pos": 104,
            "transcript_start_codon_exon": 2,
            "transcript_last_codon": "GCC",
            "transcript_valid_stop": False,
            "transcript_first_stop_codon": "TGA",
            "transcript_first_stop_pos": 164,
            "transcript_num_stop_codons": 1,
            "transcript_all_stop_codons": [(164, "TGA")],
            "transcript_stop_codon_exons": [2],
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 60,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 138,
            "stop_codon_distance": 68,
            "ptc_to_intron": 41,
            **NO_RULE,
            "nmd_50nt_penultimate_rule": True,
            "nmd_start_proximal_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(
            Mark("alt", 164, 165, "*", "PTC TGA, in the frame +1"),
            Span("alt", 104, 164, "ptc_to_start_codon = 60, in the frame +1"),
            Span("alt", 164, 205, "ptc_to_intron = 41"),
        ),
        ruler=Ruler((-4, 0, 63, 100, 160, 201, 228, 237), "CDS"),
    ),
    # NR-13: without an annotated start codon, the start-proximal rule is False, also 9 nt after the CDS start
    Case(
        "cds_start_nf_ptc_9_nt_after_the_cds_start_escapes_by_no_rule",
        """
        CDS     0           9                            63                   153                 168       177
        ref 5' [GCC GCC GCC CAG GCC GCC ..39.. GCC GCC]|[GCC GCC ..81.. GCC]|[GCC GCC GCC GCC GCC TAA gccgcc] 3'
        alt 5' [GCC GCC GCC TAG GCC GCC ..39.. GCC GCC]|[GCC GCC ..81.. GCC]|[GCC GCC GCC GCC GCC TAA gccgcc] 3'
                            ^ C>T
                            *** PTC
                            <-- ptc_to_intron = 54 -->
        exon 1: 63 nt, exon 2: 90 nt, exon 3: 24 nt
        cds_start_NF: the transcript has no start codon and no 5' UTR.
        """,
        NF,
        Change("GCC[C>T]AGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("T", "A"),
            "start_variant": per_strand(19, 217),
            "end_variant": per_strand(20, 218),
            "alt_cds_start": per_strand(10, 16),
            "alt_cds_stop": per_strand(221, 227),
            "alt_cds_seq": "GCC" * 3 + "TAG" + "GCC" * 52 + "TAA",
            "alt_cds_len": 171,
            "alt_cds_info": [(1, 63), (2, 90), (3, 18)],
            "alt_start_codon_pos": None,
            "alt_start_codon_exon": None,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 9,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(9, "TAG"), (168, "TAA")],
            "alt_stop_codon_exons": [1, 3],
            "alt_is_premature": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GCC" * 3 + "TAG" + "GCC" * 52 + "TAAGCCGCC",
            "alt_transcript_length": 177,
            "alt_cds_start_in_transcript": 0,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 0,
            "downstream_exon_count": 2,
            "ptc_to_start_codon": None,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 63,
            "stop_codon_distance": 159,
            "ptc_to_intron": 54,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_annotated_start",
        },
        marks=(Mark("alt", 9, 12, "*", "PTC"), Span("alt", 9, 63, "ptc_to_intron = 54")),
        ruler=Ruler((0, 9, 63, 153, 168, 177), "CDS"),
    ),
]
