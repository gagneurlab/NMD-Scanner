"""
Conformance cases of the null cases: each case pins a constellation where the column tables of "Technical Notes.md"
or "Input Defects.md" say "Null when", and some misannotations that make a column null. A CDS row outside the exon rows of its
transcript is an error.

The drawings show the transcript 5' to 3', also for the minus strand. The runner renders their layout block from the
case. The line ref is the transcript and the line alt is the transcript with the change. The two lines are aligned,
and `-` fills a gap. `[...]` is an exon and `|` is an exon junction. Upper case is the coding region, in codons of the
annotated frame; lower case is UTR. An intron that the change touches is drawn in lower case outside the brackets.
`^` marks the change as ref>alt. A mark line puts a character under bases, e.g. `***` under the PTC, and
`<-- label -->` spans a length. A ruler gives tx or CDS positions, as labelled.
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
    Raises,
    Ruler,
    Span,
    Transcript,
    names,
    per_strand,
)


def _exon_1_starts_at_the_second_cds_base(rows, strand):
    """Move the 5' end of exon 1 by 5 nt, past its 5' UTR and the first base of the start codon."""
    for row in rows:
        if row[2] == "exon" and row[8].endswith("exon_number=1"):
            if strand == "+":
                row[3] = str(int(row[3]) + 5)
            else:
                row[4] = str(int(row[4]) - 5)
    return rows


def _cds_row_of_exon_2_starts_in_the_acceptor(rows, strand):
    """Move the 5' end of the CDS row of exon 2 by 3 nt into the intron, over the CAG of the acceptor."""
    for row in rows:
        if row[2] == "CDS" and row[8].endswith("exon_number=2"):
            if strand == "+":
                row[3] = str(int(row[3]) - 3)
            else:
                row[4] = str(int(row[4]) + 3)
    return rows


def _exon_row_twice(rows, strand):
    """Write the exon row twice."""
    return [copy for row in rows for copy in ([row, list(row)] if row[2] == "exon" else [row])]


def _cds_row_of_exon_3_has_exon_number_4(rows, strand):
    """Give the CDS row of exon 3 the exon_number 4, which no exon row has."""
    for row in rows:
        if row[2] == "CDS" and row[8].endswith("exon_number=3"):
            row[8] = row[8].replace("exon_number=3", "exon_number=4")
    return rows


# 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
THREE_EXONS_REF = {
    **IDS,
    "cds_start": per_strand(14, 16),
    "cds_end": per_strand(75, 77),
    "ref_cds_seq": "ATGGCCAAGTGGGGCTCCTAA",
    "ref_cds_length": 21,
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "ref_cds_exons": [(1, 6), (2, 9), (3, 6)],
    "cds_in_transcript": True,
    "start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_first_stop_pos": 18,
    "ref_stop_codon_count": 1,
    "ref_stop_codons": [(18, "TAA")],
    "ref_stop_codon_exons": [3],
    "ref_has_ptc": False,
    "transcript_start": 10,
    "transcript_end": 81,
    "transcript_seq": "GACCATGGCCAAGTGGGGCTCCTAAGCCGCC",
    "transcript_length": 31,
    "cds_start_in_transcript": 4,
    "cds_end_in_transcript": 25,
    "transcript_exons": [(1, 10), (2, 9), (3, 12)],
    "utr3_length": 6,
    "utr5_length": 4,
    "total_exon_count": 3,
    "likely_misannotated": False,
}
THREE_EXONS = Layout(Transcript(("gaccATGGCC", "AAGTGGGGC", "TCCTAAgccgcc")), THREE_EXONS_REF)
# The same transcript without start_codon rows
NO_START_CODON = Layout(
    Transcript(("gaccATGGCC", "AAGTGGGGC", "TCCTAAgccgcc"), start_codon=False),
    {
        **THREE_EXONS_REF,
        "has_start_codon": False,
        "start_codon_exon": None,
        "likely_misannotated": True,
    },
)
# The CDS ends in the sense codon TCA, and there are no stop_codon rows
NO_STOP_CODON = Layout(
    Transcript(("gaccATGGCC", "AAGTGGGGC", "TCCTCAgccgcc"), stop_codon=False),
    {
        **THREE_EXONS_REF,
        "ref_cds_seq": "ATGGCCAAGTGGGGCTCCTCA",
        "has_stop_codon": False,
        "ref_last_codon": "TCA",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_stop_codon_count": 0,
        "ref_stop_codons": [],
        "ref_stop_codon_exons": [],
        "transcript_seq": "GACCATGGCCAAGTGGGGCTCCTCAGCCGCC",
        "utr3_length": None,
        "likely_misannotated": True,
    },
)
# Exon 2 starts with an in-frame ATG
SECOND_ATG = Layout(
    Transcript(("gaccATGGCC", "ATGTGGGGC", "TCCTAAgccgcc")),
    {
        **THREE_EXONS_REF,
        "ref_cds_seq": "ATGGCCATGTGGGGCTCCTAA",
        "transcript_seq": "GACCATGGCCATGTGGGGCTCCTAAGCCGCC",
    },
)
# The 3' UTR holds an ATG that starts on the last base of the stop codon TGA
ATG_ON_THE_STOP_CODON = Layout(
    Transcript(("gaccATGGCC", "AAGTGGGGC", "TCCTGAtgccgcc")),
    {
        **THREE_EXONS_REF,
        "cds_start": per_strand(14, 17),
        "cds_end": per_strand(75, 78),
        "ref_cds_seq": "ATGGCCAAGTGGGGCTCCTGA",
        "ref_last_codon": "TGA",
        "ref_first_stop_codon": "TGA",
        "ref_stop_codons": [(18, "TGA")],
        "transcript_end": 82,
        "transcript_seq": "GACCATGGCCAAGTGGGGCTCCTGATGCCGCC",
        "transcript_length": 32,
        "transcript_exons": [(1, 10), (2, 9), (3, 13)],
        "utr3_length": 7,
    },
)
# The 3' UTR holds a TGA in the frame of the CDS shifted by 1 nt
STOP_IN_THE_SHIFTED_FRAME = Layout(
    Transcript(("gaccATGGCC", "AAGTGGGGC", "TCCTAAgtgaccgcc")),
    {
        **THREE_EXONS_REF,
        "cds_start": per_strand(14, 19),
        "cds_end": per_strand(75, 80),
        "transcript_end": 84,
        "transcript_seq": "GACCATGGCCAAGTGGGGCTCCTAAGTGACCGCC",
        "transcript_length": 34,
        "transcript_exons": [(1, 10), (2, 9), (3, 15)],
        "utr3_length": 9,
    },
)
# No start_codon rows, and the 3' UTR holds a TGA in frame: gcc TGA cc
NO_START_CODON_STOP_IN_THE_UTR = Layout(
    Transcript(("gaccATGGCC", "AAGTGGGGC", "TCCTAAgcctgacc"), start_codon=False),
    {
        **THREE_EXONS_REF,
        "cds_start": per_strand(14, 18),
        "cds_end": per_strand(75, 79),
        "has_start_codon": False,
        "start_codon_exon": None,
        "transcript_end": 83,
        "transcript_seq": "GACCATGGCCAAGTGGGGCTCCTAAGCCTGACC",
        "transcript_length": 33,
        "transcript_exons": [(1, 10), (2, 9), (3, 14)],
        "utr3_length": 8,
        "likely_misannotated": True,
    },
)
# The CDS row of exon 3 has exon_number 4
CDS_ROW_WITH_ANOTHER_EXON_NUMBER = Layout(
    Transcript(("gaccATGGCC", "AAGTGGGGC", "TCGTAAgccgcc"), edit_gff3=_cds_row_of_exon_3_has_exon_number_4),
    {
        **THREE_EXONS_REF,
        "ref_cds_seq": "ATGGCCAAGTGGGGCTCGTAA",
        "ref_cds_exons": [(1, 6), (2, 9), (4, 6)],
        "ref_stop_codon_exons": [4],
        "transcript_seq": "GACCATGGCCAAGTGGGGCTCGTAAGCCGCC",
    },
)
# 5' [gacc AT gcc] 3': a CDS of 2 nt, with start_codon and stop_codon rows on its 2 nt
TWO_NT_CDS = Layout(
    Transcript(("gaccATgcc",)),
    {
        **IDS,
        "cds_start": per_strand(14, 13),
        "cds_end": per_strand(16, 15),
        "ref_cds_seq": "AT",
        "ref_cds_length": 2,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [(1, 2)],
        "cds_in_transcript": True,
        "start_codon_exon": None,
        "ref_last_codon": None,
        "ref_valid_stop": None,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_stop_codon_count": None,
        "ref_stop_codons": None,
        "ref_stop_codon_exons": None,
        "ref_has_ptc": None,
        "transcript_start": 10,
        "transcript_end": 19,
        "transcript_seq": "GACCATGCC",
        "transcript_length": 9,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 6,
        "transcript_exons": [(1, 9)],
        "utr3_length": 3,
        "utr5_length": 4,
        "total_exon_count": 1,
        "likely_misannotated": True,
    },
)
# 5' [gacc ATG CCC GGG TGA gcc] 3'
SINGLE_EXON = Layout(
    Transcript(("gaccATGCCCGGGTGAgcc",)),
    {
        **IDS,
        "cds_start": per_strand(14, 13),
        "cds_end": per_strand(26, 25),
        "ref_cds_seq": "ATGCCCGGGTGA",
        "ref_cds_length": 12,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [(1, 12)],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TGA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TGA",
        "ref_first_stop_pos": 9,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [(9, "TGA")],
        "ref_stop_codon_exons": [1],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 29,
        "transcript_seq": "GACCATGCCCGGGTGAGCC",
        "transcript_length": 19,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 16,
        "transcript_exons": [(1, 19)],
        "utr3_length": 3,
        "utr5_length": 4,
        "total_exon_count": 1,
        "likely_misannotated": False,
    },
)
# 5' [gacc GCC CGG GAA ATG A gcc] 3': 13 nt CDS without start_codon rows; its stop_codon rows are on TGA, out of frame
STOP_CODON_OUT_OF_FRAME = Layout(
    Transcript(("gaccGCCCGGGAAATGAgcc",), start_codon=False),
    {
        **IDS,
        "cds_start": per_strand(14, 13),
        "cds_end": per_strand(27, 26),
        "ref_cds_seq": "GCCCGGGAAATGA",
        "ref_cds_length": 13,
        "has_start_codon": False,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [(1, 13)],
        "cds_in_transcript": True,
        "start_codon_exon": None,
        "ref_last_codon": "TGA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_stop_codon_count": 0,
        "ref_stop_codons": [],
        "ref_stop_codon_exons": [],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 30,
        "transcript_seq": "GACCGCCCGGGAAATGAGCC",
        "transcript_length": 20,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 17,
        "transcript_exons": [(1, 20)],
        "utr3_length": 3,
        "utr5_length": 4,
        "total_exon_count": 1,
        "likely_misannotated": True,
    },
)
# 5' [ggac ATG GCC AAG TAA gcc] 3', with the exon row twice in the GFF3
EXON_ROW_TWICE = Layout(
    Transcript(("ggacATGGCCAAGTAAgcc",), edit_gff3=_exon_row_twice),
    {
        **IDS,
        "cds_start": per_strand(14, 13),
        "cds_end": per_strand(26, 25),
        "ref_cds_seq": "ATGGCCAAGTAA",
        "ref_cds_length": 12,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [(1, 12)],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 9,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [(9, "TAA")],
        "ref_stop_codon_exons": [1],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 29,
        "transcript_seq": "GGACATGGCCAAGTAAGCC" + "GGACATGGCCAAGTAAGCC",
        "transcript_length": 38,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 16,
        "transcript_exons": [(1, 19), (1, 19)],
        "utr3_length": 22,
        "utr5_length": 4,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)

# The alt CDS of a variant that changes neither the length of the CDS nor its stop codon, on THREE_EXONS
THREE_EXONS_ALT_CDS = {
    "alt_cds_length": 21,
    "alt_cds_exons": [(1, 6), (2, 9), (3, 6)],
}
TGG_TO_TTG = {
    "variant_id": "var1",
    "ref": per_strand("G", "C"),
    "alt": per_strand("T", "A"),
    "variant_start": per_strand(44, 46),
    "variant_end": per_strand(45, 47),
}
TGG_TO_TGA = {
    "variant_id": "var1",
    "ref": per_strand("G", "C"),
    "alt": per_strand("A", "T"),
    "variant_start": per_strand(45, 45),
    "variant_end": per_strand(46, 46),
}
TAA_TO_TAC = {
    "variant_id": "var1",
    "ref": per_strand("A", "T"),
    "alt": per_strand("C", "G"),
    "variant_start": per_strand(74, 16),
    "variant_end": per_strand(75, 17),
}
ATG_TO_CTG = {
    "variant_id": "var1",
    "ref": per_strand("A", "T"),
    "alt": per_strand("C", "G"),
    "variant_start": per_strand(14, 76),
    "variant_end": per_strand(15, 77),
}

# tx      0    4           11           21
# ref 5' [gacc ATG GAT G]|[TA AGC TAA gc] 3'
# CDS          0    4      7      12
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
    "ref_cds_exons": [(1, 7), (2, 8)],
    "start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_first_stop_pos": 12,
    "ref_stop_codon_count": 1,
    "ref_stop_codons": [(12, "TAA")],
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
        "transcript_exons": [(1, 11), (2, 10)],
        "utr3_length": 2,
        "utr5_length": 4,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)

# The same transcript without exon rows in the GFF3: its CDS, start_codon and stop_codon rows stay. "Input
# Defects.md" makes every column of the transcript null then, and cds_in_transcript False.
TWO_EXONS_WITHOUT_EXON_ROWS = Layout(
    Transcript(TWO_EXONS_SEQUENCES, exon_rows=False),
    {
        **TWO_EXONS_CDS,
        "cds_in_transcript": False,
        "transcript_start": None,
        "transcript_end": None,
        "transcript_seq": None,
        "transcript_length": None,
        "cds_start_in_transcript": None,
        "cds_end_in_transcript": None,
        "transcript_exons": None,
        "utr3_length": None,
        "utr5_length": None,
        "total_exon_count": None,
        "likely_misannotated": True,
    },
)


# THREE_EXONS without exon rows in the GFF3, as TWO_EXONS_WITHOUT_EXON_ROWS
THREE_EXONS_WITHOUT_EXON_ROWS = Layout(
    Transcript(THREE_EXONS.transcript.exons, exon_rows=False),
    {
        **THREE_EXONS_REF,
        "cds_in_transcript": False,
        "transcript_start": None,
        "transcript_end": None,
        "transcript_seq": None,
        "transcript_length": None,
        "cds_start_in_transcript": None,
        "cds_end_in_transcript": None,
        "transcript_exons": None,
        "utr3_length": None,
        "utr5_length": None,
        "total_exon_count": None,
        "likely_misannotated": True,
    },
)


# The missense AGC>AGA in exon 2 on TWO_EXONS, with or without exon rows: the columns of the alt CDS
MISSENSE_ALT_CDS = {
    "variant_id": "var1",
    "ref": per_strand("C", "G"),
    "alt": per_strand("A", "T"),
    "variant_start": per_strand(45, 15),
    "variant_end": per_strand(46, 16),
    "alt_cds_seq": "ATGGATGTAAGATAA",
    "alt_cds_length": 15,
    "alt_cds_exons": [(1, 7), (2, 8)],
    "alt_last_codon": "TAA",
    "alt_valid_stop": True,
    "alt_first_stop_codon": "TAA",
    "alt_first_stop_pos": 12,
    "alt_stop_codon_count": 1,
    "alt_stop_codons": [(12, "TAA")],
    "alt_stop_codon_exons": [2],
    "alt_has_ptc": False,
    "start_loss": False,
    "stop_loss": False,
    **NOT_SCANNED,
    "unknown_reason": None,
    **NO_PTC_FEATURES,
    "annotated_stop_distance": 0,
    **NO_RULE,
}

# One exon of 16 nt: the 5'UTR gacc, the CDS ATGGCCTAA and the 3'UTR tcc
SHORT_TRANSCRIPT = Layout(
    Transcript(("gaccATGGCCTAAtcc",)),
    {
        **IDS,
        "cds_start": per_strand(14, 13),
        "cds_end": per_strand(23, 22),
        "ref_cds_seq": "ATGGCCTAA",
        "ref_cds_length": 9,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [(1, 9)],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 6,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [(6, "TAA")],
        "ref_stop_codon_exons": [1],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 26,
        "transcript_seq": "GACCATGGCCTAATCC",
        "transcript_length": 16,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 13,
        "transcript_exons": [(1, 16)],
        "utr3_length": 3,
        "utr5_length": 4,
        "total_exon_count": 1,
        "likely_misannotated": False,
    },
)
# One exon whose start_codon rows mark the stop codon TAG, a misannotation
TAG_START_CODON = Layout(
    Transcript(("gccTAGGCCAAGCTGTAAgcc",)),
    {
        **IDS,
        "cds_start": 13,
        "cds_end": 28,
        "ref_cds_seq": "TAGGCCAAGCTGTAA",
        "ref_cds_length": 15,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [(1, 15)],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAG",
        "ref_first_stop_pos": 0,
        "ref_stop_codon_count": 2,
        "ref_stop_codons": [(0, "TAG"), (12, "TAA")],
        "ref_stop_codon_exons": [1, 1],
        "ref_has_ptc": True,
        "transcript_start": 10,
        "transcript_end": 31,
        "transcript_seq": "GCCTAGGCCAAGCTGTAAGCC",
        "transcript_length": 21,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 18,
        "transcript_exons": [(1, 21)],
        "utr3_length": 3,
        "utr5_length": 3,
        "total_exon_count": 1,
        "likely_misannotated": False,
    },
)

CASES = [
    # NU-12, NU-17, NU-20
    Case(
        "missense_without_start_or_stop_loss_is_not_scanned_and_has_no_ptc_features",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG TTG GGC]|[TCC TAA gccgcc] 3'
                                    ^ G>T
        G>T changes TGG to TTG: no start loss, no stop loss, no PTC.
        """,
        THREE_EXONS,
        Change("AAGT[G>T]GGG"),
        {
            **TGG_TO_TTG,
            **THREE_EXONS_ALT_CDS,
            "alt_cds_seq": "ATGGCCAAGTTGGGCTCCTAA",
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 18,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [(18, "TAA")],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATGGCCAAGTTGGGCTCCTAAGCCGCC",
            "alt_transcript_length": 31,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
    ),
    # NU-01, NU-24
    Case(
        "snv_in_the_acceptor_of_the_last_exon_makes_the_alt_columns_features_and_rules_null",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]gtaagtcccccccctttcag[TCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG TGG GGC]gtaagtcccccccctttcgg[TCC TAA gccgcc] 3'
                                                             ^ A>G
        A>G changes the acceptor AG of the last exon to GG, so unknown_reason is splice_site_destroyed.
        """,
        THREE_EXONS,
        Change("TTTC[A>G]GTCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(67, 23),
            "variant_end": per_strand(68, 24),
            **UNKNOWN_ALT,
            "unknown_reason": "splice_site_destroyed",
        },
    ),
    # NU-02
    Case(
        "without_start_codon_rows_start_codon_exon_is_null",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG TTG GGC]|[TCC TAA gccgcc] 3'
                                    ^ G>T
        The GFF3 has no start_codon rows, so has_start_codon is False. G>T changes TGG to TTG.
        """,
        NO_START_CODON,
        Change("AAGT[G>T]GGG"),
        {
            **TGG_TO_TTG,
            **THREE_EXONS_ALT_CDS,
            "alt_cds_seq": "ATGGCCAAGTTGGGCTCCTAA",
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 18,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [(18, "TAA")],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATGGCCAAGTTGGGCTCCTAAGCCGCC",
            "alt_transcript_length": 31,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
    ),
    # NU-22
    Case(
        "ptc_without_an_annotated_start_codon_has_no_ptc_to_start_codon",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG TGA GGC]|[TCC TAA gccgcc] 3'
                                     ^ G>A
                                   *** PTC
                                   <-----> ptc_to_exon_end = 6
        The GFF3 has no start_codon rows, so has_start_codon is False. G>A changes TGG to TGA, the PTC at CDS 9.
        """,
        NO_START_CODON,
        Change("AAGTG[G>A]GGC"),
        {
            **TGG_TO_TGA,
            **THREE_EXONS_ALT_CDS,
            "alt_cds_seq": "ATGGCCAAGTGAGGCTCCTAA",
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 9,
            "alt_stop_codon_count": 2,
            "alt_stop_codons": [(9, "TGA"), (18, "TAA")],
            "alt_stop_codon_exons": [2, 3],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATGGCCAAGTGAGGCTCCTAAGCCGCC",
            "alt_transcript_length": 31,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": None,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 9,
            "annotated_stop_distance": 9,
            "ptc_to_exon_end": 6,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": True,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_annotated_start",
        },
        marks=(Mark("alt", 13, 16, "*", "PTC"), Span("alt", 13, 19, "ptc_to_exon_end = 6")),
    ),
    # Pins the null clause of ptc_to_start_codon "the annotated start codon is a stop codon, such as TAG" ("Input
    # Defects.md": "A variant that leaves it unchanged gives a PTC row whose PTC is this start codon"). So
    # ptc_less_than_150nt_to_start is False ("False if `ptc_to_start_codon` is null"), and so is
    # nmd_start_proximal_rule. The row keeps the flags from the CDS: the ref transcript, read in frame, stops at the
    # TAG and not at the annotated stop codon. likely_misannotated is False ("Input Defects.md": "`likely_misannotated`
    # does not flag it, because its start codon check only asks for an annotated start codon at CDS position 0"). NU-22
    # pins the other clause, a PTC row without an annotated start codon.
    Case(
        "missense_in_a_cds_whose_annotated_start_codon_is_tag_is_a_ptc_row_without_ptc_to_start_codon",
        """
        tx      0   3               15
        ref 5' [gcc TAG GCC AAG CTG TAA gcc] 3'
        alt 5' [gcc TAG GAC AAG CTG TAA gcc] 3'
                         ^ C>A
                    *** PTC
                    <---------------------> ptc_to_exon_end = 18
        C>A at CDS 4: GCC>GAC. The start_codon rows mark TAG at CDS 0, a stop codon. It is the first in-frame stop
        codon of the alt CDS, so it is the PTC, and ptc_to_start_codon is null. The ref transcript, read in frame
        from tx 3, stops at this TAG, so the row keeps the flags from the CDS: alt_has_ptc is True.
        annotated_stop_distance = 15 - 3 = 12. likely_misannotated is False: its start codon check only asks for an
        annotated start codon at CDS position 0.
        """,
        TAG_START_CODON,
        Change("TAGG[C>A]CAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "variant_start": per_strand(17, 23),
            "variant_end": per_strand(18, 24),
            "alt_cds_seq": "TAGGACAAGCTGTAA",
            "alt_cds_length": 15,
            "alt_cds_exons": [(1, 15)],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 0,
            "alt_stop_codon_count": 2,
            "alt_stop_codons": [(0, "TAG"), (12, "TAA")],
            "alt_stop_codon_exons": [1, 1],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GCCTAGGACAAGCTGTAAGCC",
            "alt_transcript_length": 21,
            "alt_cds_start_in_transcript": 3,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 0,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": None,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 21,
            "annotated_stop_distance": 12,
            "ptc_to_exon_end": 18,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": True,
            "nmd_escape": True,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "ref_ptc",
        },
        marks=(Mark("alt", 3, 6, "*", "PTC"), Span("alt", 3, 21, "ptc_to_exon_end = 18")),
        ruler=Ruler((0, 3, 15)),
    ),
    # NU-03
    Case(
        "cds_of_two_nt_makes_every_codon_column_null",
        """
        ref 5' [gacc AT gcc] 3'
        alt 5' [gacc GT gcc] 3'
                     ^ A>G
        The GFF3 has start_codon and stop_codon rows on the 2 nt CDS.
        """,
        TWO_NT_CDS,
        Change("acc[A>G]Tgc"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(14, 14),
            "variant_end": per_strand(15, 15),
            "alt_cds_seq": "GT",
            "alt_cds_length": 2,
            "alt_cds_exons": [(1, 2)],
            "alt_last_codon": None,
            "alt_valid_stop": None,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": None,
            "alt_stop_codons": None,
            "alt_stop_codon_exons": None,
            "alt_has_ptc": None,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCGTGCC",
            "alt_transcript_length": 9,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
    ),
    # NU-04
    Case(
        "deletion_that_leaves_two_nt_of_the_cds_makes_the_alt_codon_columns_null",
        """
        ref 5' [gacc ATG CCC GGG TGA gcc] 3'
        alt 5' [gacc A-- --- --- --A gcc] 3'
                      ^^^^^^^^^^^^^ TGCCCGGGTG>-
        Deleting TGCCCGGGTG leaves the alt CDS AA. It is a start loss, and the alt transcript has no ATG.
        """,
        SINGLE_EXON,
        Change("gaccA[TGCCCGGGTG>]Agcc"),
        {
            "variant_id": "var1",
            "ref": per_strand("ATGCCCGGGTG", "TCACCCGGGCA"),
            "alt": per_strand("A", "T"),
            "variant_start": per_strand(14, 13),
            "variant_end": per_strand(25, 24),
            "alt_cds_seq": "AA",
            "alt_cds_length": 2,
            "alt_cds_exons": [(1, 2)],
            "alt_last_codon": None,
            "alt_valid_stop": None,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": None,
            "alt_stop_codons": None,
            "alt_stop_codon_exons": None,
            "alt_has_ptc": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GACCAAGCC",
            "alt_transcript_length": 9,
            "alt_cds_start_in_transcript": 4,
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
            "alt_transcript_exons": [(1, 9)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("gacc[ATGCCCGGGTG>A]Agcc"),),
    ),
    # Pins the threshold of the null clause "the CDS has fewer than 3 nt" of the alt codon columns, alt_last_codon to
    # alt_stop_codon_exons ("as ref_last_codon", "as ref_valid_stop" and so on). The alt CDS has exactly 3 nt, so these
    # columns have values. The closest case, NU-04, deletes one base more and leaves 2 nt, which makes them null.
    Case(
        "deletion_that_leaves_three_nt_of_the_cds_gives_the_alt_codon_columns",
        """
        ref 5' [gacc ATG CCC GGG TGA gcc] 3'
        alt 5' [gacc A-- --- --- -GA gcc] 3'
                      ^^^^^^^^^^^^ TGCCCGGGT>-
                     <-------------> alt_cds_length = 3
        Deleting TGCCCGGGT leaves the alt CDS AGA of 3 nt, so the alt codon columns have values: alt_last_codon is
        AGA, alt_valid_stop is False, and alt_stop_codon_count is 0. It is a start loss, and the alt transcript
        GACCAGAGCC has no ATG.
        """,
        SINGLE_EXON,
        Change("gaccA[TGCCCGGGT>]GAgcc"),
        {
            "variant_id": "var1",
            "ref": per_strand("ATGCCCGGGT", "CACCCGGGCA"),
            "alt": per_strand("A", "C"),
            "variant_start": per_strand(14, 14),
            "variant_end": per_strand(24, 24),
            "alt_cds_seq": "AGA",
            "alt_cds_length": 3,
            "alt_cds_exons": [(1, 3)],
            "alt_last_codon": "AGA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GACCAGAGCC",
            "alt_transcript_length": 10,
            "alt_cds_start_in_transcript": 4,
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
            "alt_transcript_exons": [(1, 10)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("gacc[ATGCCCGGGT>A]GAgcc"),),
        marks=(Span("alt", 4, 7, "alt_cds_length = 3"),),
    ),
    # Pins the rows that keep the flags from the CDS because their alt_transcript_seq has fewer than 3 nt ("Stop codon
    # classification"): "Such a sequence holds no codon, so it is neither classified nor scanned, also after a start
    # loss. Its alt CDS has fewer than 3 nt too. So `alt_has_ptc` is null, and `stop_loss` is True if
    # `ref_valid_stop` is True." The scan "runs only if `start_loss` or `stop_loss` is True and `alt_transcript_seq`
    # has at least 3 nt". Every other start loss with alt_transcript_seq is scanned. The next case leaves 3 nt.
    Case(
        "deletion_that_leaves_two_nt_of_the_transcript_is_not_scanned_and_keeps_the_flags_from_the_cds",
        """
        tx      0    4           13
        ref 5' [gacc ATG GCC TAA tcc] 3'
        alt 5' [g--- --- --- --- --c] 3'
                 ^^^^^^^^^^^^^^^^^^ 14 nt>-
        Deleting tx 1 to 14 leaves the alt transcript GC of 2 nt and an empty alt CDS at alt tx 1. It is a start
        loss. GC holds no codon, so the row is not scanned and keeps the flags from the CDS: alt_has_ptc is
        null, and stop_loss is True, because ref_valid_stop is True and alt_valid_stop is null.
        """,
        SHORT_TRANSCRIPT,
        Change("g[accATGGCCTAAtc>]c"),
        {
            "variant_id": "var1",
            "ref": per_strand("GACCATGGCCTAATC", "GGATTAGGCCATGGT"),
            "alt": "G",
            "variant_start": 10,
            "variant_end": 25,
            "alt_cds_seq": "",
            "alt_cds_length": 0,
            "alt_cds_exons": [(1, 0)],
            "alt_last_codon": None,
            "alt_valid_stop": None,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": None,
            "alt_stop_codons": None,
            "alt_stop_codon_exons": None,
            "alt_has_ptc": None,
            "start_loss": True,
            "stop_loss": True,
            "alt_transcript_seq": "GC",
            "alt_transcript_length": 2,
            "alt_cds_start_in_transcript": 1,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": [(1, 2)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("[gaccATGGCCTAAtc>g]c"),),
        ruler=Ruler((0, 4, 13)),
    ),
    # Pins the threshold of the scan ("PTC columns"): "An `alt_transcript_seq` of exactly 3 nt is scanned too." The
    # row is a start loss, and the scan finds no ATG: "Without an ATG, ...
    # Both flags are then False, and `annotated_stop_distance` is null." The case before leaves 2 nt, which are not
    # scanned.
    Case(
        "deletion_that_leaves_three_nt_of_the_transcript_is_scanned",
        """
        tx      0    4           13
        ref 5' [gacc ATG GCC TAA tcc] 3'
        alt 5' [g--- --- --- --- -cc] 3'
                 ^^^^^^^^^^^^^^^^^ 13 nt>-
        Deleting tx 1 to 13 leaves the alt transcript GCC of 3 nt and an empty alt CDS at alt tx 1. It is a start
        loss, and 3 nt are scanned. The scan finds no ATG, so both flags are False and annotated_stop_distance is null.
        """,
        SHORT_TRANSCRIPT,
        Change("g[accATGGCCTAAt>]cc"),
        {
            "variant_id": "var1",
            "ref": per_strand("GACCATGGCCTAAT", "GATTAGGCCATGGT"),
            "alt": "G",
            "variant_start": per_strand(10, 11),
            "variant_end": per_strand(24, 25),
            "alt_cds_seq": "",
            "alt_cds_length": 0,
            "alt_cds_exons": [(1, 0)],
            "alt_last_codon": None,
            "alt_valid_stop": None,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": None,
            "alt_stop_codons": None,
            "alt_stop_codon_exons": None,
            "alt_has_ptc": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GCC",
            "alt_transcript_length": 3,
            "alt_cds_start_in_transcript": 1,
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
            "alt_transcript_exons": [(1, 3)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("[gaccATGGCCTAAt>g]cc"),),
        ruler=Ruler((0, 4, 13)),
    ),
    # NU-05, NU-18, NU-25
    Case(
        "cds_without_stop_codon_rows_and_without_in_frame_stop_has_no_utr3_length_and_no_annotated_stop_distance",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TCA gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG TTG GGC]|[TCC TCA gccgcc] 3'
                                    ^ G>T
        The GFF3 has no stop_codon rows, and the last codon is TCA. G>T changes TGG to TTG.
        """,
        NO_STOP_CODON,
        Change("AAGT[G>T]GGG"),
        {
            **TGG_TO_TTG,
            **THREE_EXONS_ALT_CDS,
            "alt_cds_seq": "ATGGCCAAGTTGGGCTCCTCA",
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
            "alt_transcript_seq": "GACCATGGCCAAGTTGGGCTCCTCAGCCGCC",
            "alt_transcript_length": 31,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
    ),
    # NU-06
    Case(
        "frameshift_without_stop_in_the_alt_cds_has_no_alt_first_stop_codon",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gtgaccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG T-G GGC]|[TCC TAA gtgaccgcc] 3'
                                    ^ G>-
                                                      *** alt_scan_first_stop_pos = 25
        G>- deletes one G of the run GGGG and shifts the frame. In the alt transcript, the first in-frame stop codon
        is the TGA at tx 25 (`*`). It lies in the 3' UTR: a stop loss.
        """,
        STOP_IN_THE_SHIFTED_FRAME,
        Change("AAGT[G>]GGGC"),
        {
            "variant_id": "var1",
            "ref": per_strand("TG", "CC"),
            "alt": per_strand("T", "C"),
            "variant_start": per_strand(43, 48),
            "variant_end": per_strand(45, 50),
            "alt_cds_seq": "ATGGCCAAGTGGGCTCCTAA",
            "alt_cds_length": 20,
            "alt_cds_exons": [(1, 6), (2, 8), (3, 6)],
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
            "alt_transcript_seq": "GACCATGGCCAAGTGGGCTCCTAAGTGACCGCC",
            "alt_transcript_length": 33,
            "alt_cds_start_in_transcript": 4,
            "alt_scan_start_codon_pos": 4,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TGA",
            "alt_scan_first_stop_pos": 25,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [(25, "TGA")],
            "alt_scan_stop_codon_exons": [3],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -4,
            **NO_RULE,
            "alt_transcript_exons": [(1, 10), (2, 8), (3, 15)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("AAGTGG[G>]GC"), Change("AAGTGGG[G>]C")),
        marks=(Mark("alt", 25, 28, "*", "alt_scan_first_stop_pos = 25"),),
    ),
    # NU-07
    Case(
        "snv_in_the_start_codon_makes_the_alt_start_codon_columns_null",
        """
        tx      0    4         10                22        31
        ref 5' [gacc ATG GCC]|[ATG TGG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc ATA GCC]|[ATG TGG GGC]|[TCC TAA gccgcc] 3'
                       ^ G>A
        G>A changes ATG to ATA: a start loss. The scan finds the ATG at tx 10, in frame with the TAA at tx 22.
        """,
        SECOND_ATG,
        Change("ccAT[G>A]GCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "variant_start": per_strand(16, 74),
            "variant_end": per_strand(17, 75),
            **THREE_EXONS_ALT_CDS,
            "alt_cds_seq": "ATAGCCATGTGGGGCTCCTAA",
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 18,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [(18, "TAA")],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATAGCCATGTGGGGCTCCTAAGCCGCC",
            "alt_transcript_length": 31,
            "alt_cds_start_in_transcript": 4,
            "alt_scan_start_codon_pos": 10,
            "alt_scan_start_codon_exon": 2,
            "alt_scan_first_stop_codon": "TAA",
            "alt_scan_first_stop_pos": 22,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [(22, "TAA")],
            "alt_scan_stop_codon_exons": [3],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 4, 10, 22, 31)),
    ),
    # NU-14, NU-27 (no ATG)
    Case(
        "start_loss_without_atg_in_the_alt_transcript_has_no_scan_start_codon",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc CTG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
                     ^ A>C
        A>C changes ATG to CTG: a start loss. The alt transcript has no ATG.
        """,
        THREE_EXONS,
        Change("acc[A>C]TGG"),
        {
            **ATG_TO_CTG,
            **THREE_EXONS_ALT_CDS,
            "alt_cds_seq": "CTGGCCAAGTGGGGCTCCTAA",
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 18,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [(18, "TAA")],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GACCCTGGCCAAGTGGGGCTCCTAAGCCGCC",
            "alt_transcript_length": 31,
            "alt_cds_start_in_transcript": 4,
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
    ),
    # NU-27 (ATG downstream of the first base of the annotated stop codon)
    Case(
        "start_loss_with_the_next_atg_on_the_last_base_of_the_stop_codon_has_no_annotated_stop_distance",
        """
        tx      0    4                           22         32
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TGA tgccgcc] 3'
        alt 5' [gacc CTG GCC]|[AAG TGG GGC]|[TCC TGA tgccgcc] 3'
                     ^ A>C
                                                   aaaa ATG at tx 24
        A>C changes ATG to CTG: a start loss. The next ATG is at tx 24 (`a`), on the A of the TGA at tx 22. So no
        ORF overlaps the CDS.
        """,
        ATG_ON_THE_STOP_CODON,
        Change("acc[A>C]TGG"),
        {
            **ATG_TO_CTG,
            "variant_start": per_strand(14, 77),
            "variant_end": per_strand(15, 78),
            **THREE_EXONS_ALT_CDS,
            "alt_cds_seq": "CTGGCCAAGTGGGGCTCCTGA",
            "alt_last_codon": "TGA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 18,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [(18, "TGA")],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GACCCTGGCCAAGTGGGGCTCCTGATGCCGCC",
            "alt_transcript_length": 32,
            "alt_cds_start_in_transcript": 4,
            "alt_scan_start_codon_pos": 24,
            "alt_scan_start_codon_exon": 3,
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
        marks=(Mark("alt", 24, 27, "a", "ATG at tx 24"),),
        ruler=Ruler((0, 4, 22, 32)),
    ),
    # NU-09, NU-13, NU-19, NU-23, PF-25 (5' CDS base outside the exons), ST-21
    Case(
        "cds_row_that_starts_before_exon_1_is_an_error_that_names_the_transcript_and_the_cds_row",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG TGA GGC]|[TCC TAA gccgcc] 3'
                                     ^ G>A
                     <-----> the CDS row of exon 1
                      <----> the exon row of exon 1 in the GFF3
        The exon row of exon 1 starts at the T of ATG, so the CDS row of exon 1 starts 1 nt before its exon row. The
        block draws exon 1 before this edit of the GFF3. annotate() rejects the GFF3 with a ValueError that names
        the transcript tx1 and the CDS row.
        """,
        Layout(
            Transcript(THREE_EXONS.transcript.exons, edit_gff3=_exon_1_starts_at_the_second_cds_base), THREE_EXONS.ref
        ),
        Change("AAGTG[G>A]GGC"),
        Raises(
            ValueError,
            match=per_strand(
                names("CDS row of transcript tx1 at chr1:15-20"), names("CDS row of transcript tx1 at chr1:72-77")
            ),
        ),
        marks=(Span("ref", 4, 10, "the CDS row of exon 1"), Span("ref", 5, 10, "the exon row of exon 1 in the GFF3")),
    ),
    # NU-10
    Case(
        "cds_row_that_starts_in_the_acceptor_is_an_error_that_names_the_transcript_and_the_cds_row",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG TTG GGC]|[TCC TAA gccgcc] 3'
                                    ^ G>T
        The CDS row of exon 2 starts 3 nt before its exon row, on the CAG of the acceptor. The block does not draw
        this intron, since the change does not touch it. annotate() rejects the GFF3 with a ValueError that names
        the transcript tx1 and the CDS row.
        """,
        Layout(
            Transcript(THREE_EXONS.transcript.exons, edit_gff3=_cds_row_of_exon_2_starts_in_the_acceptor),
            THREE_EXONS.ref,
        ),
        Change("AAGT[G>T]GGG"),
        Raises(
            ValueError,
            match=per_strand(
                names("CDS row of transcript tx1 at chr1:38-49"), names("CDS row of transcript tx1 at chr1:43-54")
            ),
        ),
    ),
    # NU-11
    Case(
        "deleted_5utr_bases_that_the_transcript_does_not_hold_leave_the_alt_transcript_null",
        """
        ref 5' [ggac ATG GCC AAG TAA gcc] 3'
        alt 5' [gga- -TG GCC AAG TAA gcc] 3'
                   ^^^ cA>-
        The GFF3 has the exon row twice. The alt line draws the deletion cA>-. Its other placement ac>- lies in the
        5' UTR: the alt keeps ATG GCC AAG TAA, and the 5' UTR loses AC.
        transcript_seq holds the exon twice. The deleted 5' UTR bases come from both exon rows, ACAC, and
        transcript_seq does not hold ACAC before the CDS.
        """,
        EXON_ROW_TWICE,
        Change("gga[cA>]TGG"),
        {
            "variant_id": "var1",
            "ref": per_strand("ACA", "ATG"),
            "alt": per_strand("A", "A"),
            "variant_start": per_strand(12, 23),
            "variant_end": per_strand(15, 26),
            "alt_cds_seq": "ATGGCCAAGTAA",
            "alt_cds_length": 12,
            "alt_cds_exons": [(1, 12)],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 9,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [(9, "TAA")],
            "alt_stop_codon_exons": [1],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": None,
            "alt_transcript_length": None,
            "alt_cds_start_in_transcript": None,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": None,
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("gg[aC>]ATGG"),),
    ),
    # NU-15
    Case(
        "stop_loss_without_an_annotated_start_codon_has_no_scan_start_codon",
        """
        tx      0    4                           22     28   33
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gcctgacc] 3'
        alt 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAC gcctgacc] 3'
                                                   ^ A>C
        The GFF3 has no start_codon rows, so has_start_codon is False. A>C changes TAA to TAC: a stop loss. The
        next in-frame stop codon is the TGA at tx 28.
        """,
        NO_START_CODON_STOP_IN_THE_UTR,
        Change("TCCTA[A>C]GCC"),
        {
            **TAA_TO_TAC,
            "variant_start": per_strand(74, 18),
            "variant_end": per_strand(75, 19),
            **THREE_EXONS_ALT_CDS,
            "alt_cds_seq": "ATGGCCAAGTGGGGCTCCTAC",
            "alt_last_codon": "TAC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": "GACCATGGCCAAGTGGGGCTCCTACGCCTGACC",
            "alt_transcript_length": 33,
            "alt_cds_start_in_transcript": 4,
            "alt_scan_start_codon_pos": None,
            "alt_scan_start_codon_exon": None,
            "alt_scan_first_stop_codon": "TGA",
            "alt_scan_first_stop_pos": 28,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [(28, "TGA")],
            "alt_scan_stop_codon_exons": [3],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -6,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 4, 22, 28, 33)),
    ),
    # NU-16, NU-26
    Case(
        "nonstop_stop_loss_has_no_first_stop_codon_and_no_annotated_stop_distance",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAC gccgcc] 3'
                                                   ^ A>C
        A>C changes TAA to TAC: a stop loss, with no in-frame stop codon up to the transcript end.
        """,
        THREE_EXONS,
        Change("TCCTA[A>C]GCC"),
        {
            **TAA_TO_TAC,
            **THREE_EXONS_ALT_CDS,
            "alt_cds_seq": "ATGGCCAAGTGGGGCTCCTAC",
            "alt_last_codon": "TAC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": "GACCATGGCCAAGTGGGGCTCCTACGCCGCC",
            "alt_transcript_length": 31,
            "alt_cds_start_in_transcript": 4,
            "alt_scan_start_codon_pos": 4,
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
    ),
    # NU-21
    Case(
        "ptc_in_an_exon_whose_cds_row_has_another_exon_number_takes_the_exon_features_from_the_alt_transcript",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCG TAA gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG TGG GGC]|[TAG TAA gccgcc] 3'
                                              ^ C>A
                                             *** PTC
                                             <------------> ptc_to_exon_end = 12
        The CDS row of exon 3 has exon_number 4. C>A changes TCG to TAG, the PTC at CDS 15, in "exon 4".
        The exon features read the exons of the alt transcript: the PTC at tx 19 lies in exon 3, the last exon.
        So upstream_exon_count = 2, downstream_exon_count = 0 and ptc_exon_length = 12.
        """,
        CDS_ROW_WITH_ANOTHER_EXON_NUMBER,
        Change("CAGT[C>A]GTAA"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "variant_start": per_strand(70, 20),
            "variant_end": per_strand(71, 21),
            **THREE_EXONS_ALT_CDS,
            "alt_cds_exons": [(1, 6), (2, 9), (4, 6)],
            "alt_cds_seq": "ATGGCCAAGTGGGGCTAGTAA",
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 15,
            "alt_stop_codon_count": 2,
            "alt_stop_codons": [(15, "TAG"), (18, "TAA")],
            "alt_stop_codon_exons": [4, 4],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATGGCCAAGTGGGGCTAGTAAGCCGCC",
            "alt_transcript_length": 31,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": 2,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 15,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 12,
            "annotated_stop_distance": 3,
            "ptc_to_exon_end": 12,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "ok",
        },
        marks=(Mark("alt", 19, 22, "*", "PTC"), Span("alt", 19, 31, "ptc_to_exon_end = 12")),
    ),
    # NU-28
    Case(
        "missense_on_a_cds_whose_stop_codon_is_out_of_frame_has_no_annotated_stop_distance",
        """
        CDS          0            10
        ref 5' [gacc GCC CGG GAA ATG A gcc] 3'
        alt 5' [gacc GAC CGG GAA ATG A gcc] 3'
                      ^ C>A
        The GFF3 has stop_codon rows on the TGA at CDS 10, out of frame. C>A changes GCC to GAC. The row keeps the
        flags from the CDS, and the alt CDS has no in-frame stop codon.
        """,
        STOP_CODON_OUT_OF_FRAME,
        Change("ACCG[C>A]CCGG"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "variant_start": per_strand(15, 24),
            "variant_end": per_strand(16, 25),
            "alt_cds_seq": "GACCGGGAAATGA",
            "alt_cds_length": 13,
            "alt_cds_exons": [(1, 13)],
            "alt_last_codon": "TGA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCGACCGGGAAATGAGCC",
            "alt_transcript_length": 20,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 10), "CDS"),
    ),
    # NU-29
    Case(
        "deletion_that_leaves_two_nt_of_a_cds_with_an_out_of_frame_stop_codon_has_a_null_alt_has_ptc",
        """
        ref 5' [gacc GCC CGG GAA ATG A gcc] 3'
        alt 5' [gacc G-- --- --- --- A gcc] 3'
                      ^^^^^^^^^^^^^^ 11 nt>-
        The GFF3 has stop_codon rows on TGA, out of frame: the row keeps the flags from the CDS. Deleting
        CCCGGGAAATG leaves the alt CDS GA.
        """,
        STOP_CODON_OUT_OF_FRAME,
        Change("gaccG[CCCGGGAAATG>]Agcc"),
        {
            "variant_id": "var1",
            "ref": per_strand("GCCCGGGAAATG", "TCATTTCCCGGG"),
            "alt": per_strand("G", "T"),
            "variant_start": per_strand(14, 13),
            "variant_end": per_strand(26, 25),
            "alt_cds_seq": "GA",
            "alt_cds_length": 2,
            "alt_cds_exons": [(1, 2)],
            "alt_last_codon": None,
            "alt_valid_stop": None,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": None,
            "alt_stop_codons": None,
            "alt_stop_codon_exons": None,
            "alt_has_ptc": None,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": "GACCGAGCC",
            "alt_transcript_length": 9,
            "alt_cds_start_in_transcript": 4,
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
            "alt_transcript_exons": [(1, 9)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("gacc[GCCCGGGAAAT>]GAgcc"),),
    ),
    # NU-08, MI-17. The only transcript has no exon rows. The column table gives null transcript columns,
    # cds_in_transcript False, likely_misannotated True, and the flags from the CDS (alt_transcript_seq is null).
    Case(
        "missense_in_a_transcript_without_exon_rows_gives_null_transcript_columns",
        """
        ref 5' [gacc ATG GAT G]|[TA AGC TAA gc] 3'
        alt 5' [gacc ATG GAT G]|[TA AGA TAA gc] 3'
                                      ^ C>A
        The GFF3 has CDS rows only, no exon rows. C>A changes AGC to AGA.
        """,
        TWO_EXONS_WITHOUT_EXON_ROWS,
        Change("AAG[C>A]TAA"),
        {
            **MISSENSE_ALT_CDS,
            "alt_transcript_seq": None,
            "alt_transcript_length": None,
            "alt_cds_start_in_transcript": None,
            "alt_transcript_exons": None,
            "nmd_model_status": "no_ptc",
        },
    ),
    # A PTC row on a transcript without exon rows. alt_transcript_exons is null, so the PTC has no exon features,
    # and the rules that read them are False. ptc_to_start_codon reads CDS positions only: 9 - 0.
    Case(
        "nonsense_in_a_transcript_without_exon_rows_has_no_exon_features",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG TGA GGC]|[TCC TAA gccgcc] 3'
                                     ^ G>A
                                   *** PTC
        The GFF3 has CDS rows only, no exon rows. G>A changes TGG to TGA, the PTC at CDS 9. alt_transcript_exons
        is null, so the PTC has no exon features.
        """,
        THREE_EXONS_WITHOUT_EXON_ROWS,
        Change("AAGTG[G>A]GGC"),
        {
            **TGG_TO_TGA,
            **THREE_EXONS_ALT_CDS,
            "alt_cds_seq": "ATGGCCAAGTGAGGCTCCTAA",
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 9,
            "alt_stop_codon_count": 2,
            "alt_stop_codons": [(9, "TGA"), (18, "TAA")],
            "alt_stop_codon_exons": [2, 3],
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": None,
            "alt_transcript_length": None,
            "alt_cds_start_in_transcript": None,
            "alt_transcript_exons": None,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": None,
            "downstream_exon_count": None,
            "ptc_to_start_codon": 9,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": None,
            "annotated_stop_distance": 9,
            "ptc_to_exon_end": None,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
            "nmd_model_status": "missing_input",
        },
        marks=(Mark("alt", 13, 16, "*", "PTC"),),
    ),
    # Pins the null clause of annotated_stop_distance "on a row that keeps the flags from the CDS, the alt CDS has no
    # in-frame stop codon" on a row without alt_transcript_seq. The other cases of the clause have alt_transcript_seq.
    # The transcript has no exon rows, so the row keeps the flags from the CDS: stop_loss is True, because
    # ref_valid_stop is True and alt_valid_stop is False.
    Case(
        "stop_loss_in_a_transcript_without_exon_rows_has_no_annotated_stop_distance",
        """
        ref 5' [gacc ATG GAT G]|[TA AGC TAA gc] 3'
        alt 5' [gacc ATG GAT G]|[TA AGC CAA gc] 3'
                                        ^ T>C
        The GFF3 has CDS rows only, no exon rows. T>C changes the stop codon TAA to CAA, and the alt CDS has no
        in-frame stop codon. alt_transcript_seq is null, so the row keeps the flags from the CDS: stop_loss is
        True, and annotated_stop_distance is null.
        """,
        TWO_EXONS_WITHOUT_EXON_ROWS,
        Change("AGC[T>C]AAGC"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(46, 14),
            "variant_end": per_strand(47, 15),
            "alt_cds_seq": "ATGGATGTAAGCCAA",
            "alt_cds_length": 15,
            "alt_cds_exons": [(1, 7), (2, 8)],
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
            "alt_transcript_seq": None,
            "alt_transcript_length": None,
            "alt_cds_start_in_transcript": None,
            "alt_transcript_exons": None,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "nmd_model_status": "no_ptc",
        },
    ),
    # A start loss on a transcript without exon rows: the row has no alt_transcript_seq, so the scan does not run, and
    # the row keeps the PTC of the alt CDS. ptc_to_start_codon is null, so the model cannot score the row: start_lost.
    Case(
        "start_loss_in_a_transcript_without_exon_rows_keeps_the_ptc_of_the_alt_cds_and_is_start_lost",
        """
        ref 5' [gacc ATG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc TAG GCC]|[AAG TGG GGC]|[TCC TAA gccgcc] 3'
                     ^^ AT>TA
                     *** PTC
        The GFF3 has CDS rows only, no exon rows. AT>TA changes the start codon ATG to TAG: a start loss. The row has
        no alt_transcript_seq, so the scan of the alt transcript does not run, and the row keeps the flags from the
        alt CDS. Its first in-frame stop codon is the TAG at CDS 0, upstream of the annotated stop codon at CDS 18: a
        PTC. No ATG of a scan starts the ORF, so ptc_to_start_codon is null, and nmd_model_status is start_lost.
        annotated_stop_distance is 18 - 0 = 18, in alt CDS coordinates.
        """,
        THREE_EXONS_WITHOUT_EXON_ROWS,
        Change("GACC[AT>TA]GGCC"),
        {
            "variant_id": "var1",
            "ref": "AT",
            "alt": "TA",
            "variant_start": per_strand(14, 75),
            "variant_end": per_strand(16, 77),
            **THREE_EXONS_ALT_CDS,
            "alt_cds_seq": "TAGGCCAAGTGGGGCTCCTAA",
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 0,
            "alt_stop_codon_count": 2,
            "alt_stop_codons": [(0, "TAG"), (18, "TAA")],
            "alt_stop_codon_exons": [1, 3],
            "alt_has_ptc": True,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": None,
            "alt_transcript_length": None,
            "alt_cds_start_in_transcript": None,
            "alt_transcript_exons": None,
            **NOT_SCANNED,
            "unknown_reason": None,
            "upstream_exon_count": None,
            "downstream_exon_count": None,
            "ptc_to_start_codon": None,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": None,
            "annotated_stop_distance": 18,
            "ptc_to_exon_end": None,
            **NO_RULE,
            "nmd_model_status": "start_lost",
        },
        marks=(Mark("alt", 4, 7, "*", "PTC"),),
    ),
    # The GFF3 has the exon row twice, and the deletion lies in both copies: each exon of the alt transcript has 18 nt,
    # 36 nt in all. alt_transcript_seq holds the deletion once and has 37 nt. The exon lengths do not add up, so
    # alt_transcript_exons is null, and so are the exon numbers of the scan. The frameshift reads on into the second
    # copy of the exon, to the TAA at tx 31: a stop loss.
    Case(
        "frameshift_in_an_exon_with_two_exon_rows_has_no_alt_transcript_exons",
        """
        ref 5' [ggac ATG GCC AAG TAA gcc] 3'
        alt 5' [ggac ATG -CC AAG TAA gcc] 3'
                         ^ G>-
        The GFF3 has the exon row twice, so transcript_seq holds the exon twice. Deleting one G of GG shifts the
        frame. The deletion lies in both exon rows: the exons of the alt transcript have 18 + 18 = 36 nt, but
        alt_transcript_seq has 37 nt. So alt_transcript_exons is null. The scan reads on into the second copy
        of the exon, to the TAA at tx 31: a stop loss.
        """,
        EXON_ROW_TWICE,
        Change("ATG[G>]CCAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("GG", "GC"),
            "alt": per_strand("G", "G"),
            "variant_start": per_strand(16, 20),
            "variant_end": per_strand(18, 22),
            "alt_cds_seq": "ATGCCAAGTAA",
            "alt_cds_length": 11,
            "alt_cds_exons": [(1, 11)],
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
            "alt_transcript_seq": "GGACATGCCAAGTAAGCC" + "GGACATGGCCAAGTAAGCC",
            "alt_transcript_length": 37,
            "alt_cds_start_in_transcript": 4,
            "alt_transcript_exons": None,
            "alt_scan_start_codon_pos": 4,
            "alt_scan_start_codon_exon": None,
            "alt_scan_first_stop_codon": "TAA",
            "alt_scan_first_stop_pos": 31,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [(31, "TAA")],
            "alt_scan_stop_codon_exons": None,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -19,
            **NO_RULE,
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("AT[G>]GCCAAG"),),
    ),
]
