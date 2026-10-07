"""
Conformance cases of variant placement at the edges of the coding region and at exon boundaries
("Technical Notes.md", section "Variants at exon boundaries"). The cases cover the deletion of a whole exon or intron,
splice sites of exons whose edge is also an edge of the coding region or of no coding row, intron variants that
touch nothing, insertions and delins at the start codon and at the stop codon, UTR changes and
alt_cds_start_in_transcript.

The drawings show the transcript 5' to 3', also for the minus strand. The runner renders their layout block from the
case. The line ref is the transcript and the line alt is the transcript with the change. The two lines are aligned,
and `-` fills a gap. `[...]` is an exon and `|` is an exon junction. Upper case is the coding region, in codons of the
annotated frame; lower case is UTR. An intron or a flank that the change touches is drawn in lower case outside the
brackets. `..N..` leaves out N bases. `^` marks the change as ref>alt, with `-` for no bases. `<-- 12 nt -->` spans a
length. A ruler gives tx positions, 0-based as in alt_transcript_seq; "alt tx" numbers are positions of the alt line.
A "placement" is one of the equivalent positions of an indel, or one of the two matchings of a delins. The alt line
draws one placement; the prose says which placement the tool takes, and "alt [...]" writes out its alt transcript.
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
    NoRow,
    Ruler,
    Span,
    Transcript,
    per_strand,
)

THREE_EXONS = Layout(
    Transcript(("gaccATGGCC", "AAGCTGGGC", "TCCTAAgccgcc")),
    {
        **IDS,
        "cds_start": per_strand(14, 16),
        "cds_end": per_strand(75, 77),
        "ref_cds_seq": "ATGGCCAAGCTGGGCTCCTAA",
        "ref_cds_length": 21,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [
            {"exon_number": 1, "length": 6},
            {"exon_number": 2, "length": 9},
            {"exon_number": 3, "length": 6},
        ],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 18,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [{"position": 18, "codon": "TAA"}],
        "ref_stop_codon_exons": [3],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 81,
        "transcript_seq": "GACCATGGCCAAGCTGGGCTCCTAAGCCGCC",
        "transcript_length": 31,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 25,
        "transcript_exons": [
            {"exon_number": 1, "length": 10},
            {"exon_number": 2, "length": 9},
            {"exon_number": 3, "length": 12},
        ],
        "utr3_length": 6,
        "utr5_length": 4,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# THREE_EXONS with a 3'UTR that starts with the run aaa, so the stop codon TAA ends in a run of 5 A
UTR3_A_RUN = Layout(
    Transcript(("gaccATGGCC", "AAGCTGGGC", "TCCTAAaaacctagcc")),
    {
        **IDS,
        "cds_start": per_strand(14, 20),
        "cds_end": per_strand(75, 81),
        "ref_cds_seq": "ATGGCCAAGCTGGGCTCCTAA",
        "ref_cds_length": 21,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [
            {"exon_number": 1, "length": 6},
            {"exon_number": 2, "length": 9},
            {"exon_number": 3, "length": 6},
        ],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 18,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [{"position": 18, "codon": "TAA"}],
        "ref_stop_codon_exons": [3],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 85,
        "transcript_seq": "GACCATGGCCAAGCTGGGCTCCTAAAAACCTAGCC",
        "transcript_length": 35,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 25,
        "transcript_exons": [
            {"exon_number": 1, "length": 10},
            {"exon_number": 2, "length": 9},
            {"exon_number": 3, "length": 16},
        ],
        "utr3_length": 10,
        "utr5_length": 4,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# THREE_EXONS with a 3'UTR that starts with cctag: the last CDS base C repeats as the second 3'UTR base
UTR3_CCTAG = Layout(
    Transcript(("gaccATGGCC", "AAGCTGGGC", "TCCTAAcctagcctgacc")),
    {
        **IDS,
        "cds_start": per_strand(14, 22),
        "cds_end": per_strand(75, 83),
        "ref_cds_seq": "ATGGCCAAGCTGGGCTCCTAA",
        "ref_cds_length": 21,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [
            {"exon_number": 1, "length": 6},
            {"exon_number": 2, "length": 9},
            {"exon_number": 3, "length": 6},
        ],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 18,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [{"position": 18, "codon": "TAA"}],
        "ref_stop_codon_exons": [3],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 87,
        "transcript_seq": "GACCATGGCCAAGCTGGGCTCCTAACCTAGCCTGACC",
        "transcript_length": 37,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 25,
        "transcript_exons": [
            {"exon_number": 1, "length": 10},
            {"exon_number": 2, "length": 9},
            {"exon_number": 3, "length": 18},
        ],
        "utr3_length": 12,
        "utr5_length": 4,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# Three coding exons; exon 2 has 12 nt and ends with the G that also ends intron 1
SHORT_EXON_2 = Layout(
    Transcript(("gccacATGGCCAAGCTGCTG", "CAGCAGCTGCTG", "AAGCTGTAAgcccccccc")),
    {
        **IDS,
        "cds_start": per_strand(15, 19),
        "cds_end": per_strand(91, 95),
        "ref_cds_seq": "ATGGCCAAGCTGCTGCAGCAGCTGCTGAAGCTGTAA",
        "ref_cds_length": 36,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [
            {"exon_number": 1, "length": 15},
            {"exon_number": 2, "length": 12},
            {"exon_number": 3, "length": 9},
        ],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 33,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [{"position": 33, "codon": "TAA"}],
        "ref_stop_codon_exons": [3],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 100,
        "transcript_seq": "GCCACATGGCCAAGCTGCTGCAGCAGCTGCTGAAGCTGTAAGCCCCCCCC",
        "transcript_length": 50,
        "cds_start_in_transcript": 5,
        "cds_end_in_transcript": 41,
        "transcript_exons": [
            {"exon_number": 1, "length": 20},
            {"exon_number": 2, "length": 12},
            {"exon_number": 3, "length": 18},
        ],
        "utr3_length": 9,
        "utr5_length": 5,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# The coding region is exon 2: it begins with the start codon and ends with the stop codon. Exons 1 and 3 are UTR.
CODING_EXON_2 = Layout(
    Transcript(("gccaccgcag", "ATGGCCAAGCTGCTGAAGCTGTAA", "gccccccccccc")),
    {
        **IDS,
        "cds_start": per_strand(40, 42),
        "cds_end": per_strand(64, 66),
        "ref_cds_seq": "ATGGCCAAGCTGCTGAAGCTGTAA",
        "ref_cds_length": 24,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [{"exon_number": 2, "length": 24}],
        "cds_in_transcript": True,
        "start_codon_exon": 2,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 21,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [{"position": 21, "codon": "TAA"}],
        "ref_stop_codon_exons": [2],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 96,
        "transcript_seq": "GCCACCGCAGATGGCCAAGCTGCTGAAGCTGTAAGCCCCCCCCCCC",
        "transcript_length": 46,
        "cds_start_in_transcript": 10,
        "cds_end_in_transcript": 34,
        "transcript_exons": [
            {"exon_number": 1, "length": 10},
            {"exon_number": 2, "length": 24},
            {"exon_number": 3, "length": 12},
        ],
        "utr3_length": 12,
        "utr5_length": 10,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# Exon 1 is 5'UTR only. Intron 1 has 4 nt (GTAG), so its donor GT ends 2 nt before the coding row of exon 2.
UTR_EXON_1_SHORT_INTRON = Layout(
    Transcript(("gccaccgcag", "ATGGCCAAGCTGTAAgccgcc"), introns=("GTAG",)),
    {
        **IDS,
        "cds_start": per_strand(24, 16),
        "cds_end": per_strand(39, 31),
        "ref_cds_seq": "ATGGCCAAGCTGTAA",
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
        "transcript_start": 10,
        "transcript_end": 45,
        "transcript_seq": "GCCACCGCAGATGGCCAAGCTGTAAGCCGCC",
        "transcript_length": 31,
        "cds_start_in_transcript": 10,
        "cds_end_in_transcript": 25,
        "transcript_exons": [{"exon_number": 1, "length": 10}, {"exon_number": 2, "length": 21}],
        "utr3_length": 6,
        "utr5_length": 10,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)

# One exon whose 5'UTR ends in a, the first base of the start codon ATG
UTR5_ENDS_IN_A = Layout(
    Transcript(("ggaccaATGCTGCTGTAAggccggccgg",)),
    {
        **IDS,
        "cds_start": per_strand(16, 20),
        "cds_end": per_strand(28, 32),
        "ref_cds_seq": "ATGCTGCTGTAA",
        "ref_cds_length": 12,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [{"exon_number": 1, "length": 12}],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 9,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [{"position": 9, "codon": "TAA"}],
        "ref_stop_codon_exons": [1],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 38,
        "transcript_seq": "GGACCAATGCTGCTGTAAGGCCGGCCGG",
        "transcript_length": 28,
        "cds_start_in_transcript": 6,
        "cds_end_in_transcript": 18,
        "transcript_exons": [{"exon_number": 1, "length": 28}],
        "utr3_length": 10,
        "utr5_length": 6,
        "total_exon_count": 1,
        "likely_misannotated": False,
    },
)

# One exon whose 5'UTR ends in t, the second base of the start codon ATG
UTR5_ENDS_IN_T = Layout(
    Transcript(("ggacctATGCTGCTGTAAggccggccgg",)),
    {
        **IDS,
        "cds_start": per_strand(16, 20),
        "cds_end": per_strand(28, 32),
        "ref_cds_seq": "ATGCTGCTGTAA",
        "ref_cds_length": 12,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [{"exon_number": 1, "length": 12}],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 9,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [{"position": 9, "codon": "TAA"}],
        "ref_stop_codon_exons": [1],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 38,
        "transcript_seq": "GGACCTATGCTGCTGTAAGGCCGGCCGG",
        "transcript_length": 28,
        "cds_start_in_transcript": 6,
        "cds_end_in_transcript": 18,
        "transcript_exons": [{"exon_number": 1, "length": 28}],
        "utr3_length": 10,
        "utr5_length": 6,
        "total_exon_count": 1,
        "likely_misannotated": False,
    },
)

# One exon with the non-ATG start codon CTG, which the start_codon rows mark. The 5'UTR ends in a, so a C inserted
# before the CTG can shift 1 nt into the CDS, but not into the 5'UTR.
CTG_START_AFTER_A = Layout(
    Transcript(("gccaCTGGCCAAGCTGTAAgcc",)),
    {
        **IDS,
        "cds_start": per_strand(14, 13),
        "cds_end": per_strand(29, 28),
        "ref_cds_seq": "CTGGCCAAGCTGTAA",
        "ref_cds_length": 15,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [{"exon_number": 1, "length": 15}],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 12,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [{"position": 12, "codon": "TAA"}],
        "ref_stop_codon_exons": [1],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 32,
        "transcript_seq": "GCCACTGGCCAAGCTGTAAGCC",
        "transcript_length": 22,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 19,
        "transcript_exons": [{"exon_number": 1, "length": 22}],
        "utr3_length": 3,
        "utr5_length": 4,
        "total_exon_count": 1,
        "likely_misannotated": False,
    },
)

TOUCHES_NO_CODING_REGION = NoRow("touches no coding region")

CASES = [
    # EB-25
    Case(
        "deletion_of_a_whole_short_coding_exon_that_keeps_both_splice_sites_empties_the_exon",
        """
        ref 5' [gccac ATG GCC AAG CTG CTG]gtaagtcccccccctttcag[CAG CAG CTG CTG]|[AAG CTG TAA gcccccccc] 3'
        alt 5' [gccac ATG GCC AAG CTG CTG]gtaagtcccccccctttcag[--- --- --- ---]|[AAG CTG TAA gcccccccc] 3'
                                                               ^^^^^^^^^^^^^^^ 12 nt>-
                <-------- 20 nt -------->
                                                               <--- 12 nt --->
                                                                                 <------ 18 nt ------>
        The acceptor AG and the donor GT stay. The placement 1 nt to the left deletes the G of the AG instead of
        the last G of exon 2. It gives the same sequence but loses the acceptor, so it does not count.
        alt CDS: ATGGCCAAGCTGCTG|AAGCTGTAA, 4 codons shorter, in frame
        """,
        SHORT_EXON_2,
        Change("TTTCAG[CAGCAGCTGCTG>]GTAAGT"),
        {
            "variant_id": "var1",
            "ref": per_strand("GCAGCAGCTGCTG", "CCAGCAGCTGCTG"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(49, 47),
            "variant_end": per_strand(62, 60),
            "alt_cds_seq": "ATGGCCAAGCTGCTGAAGCTGTAA",
            "alt_cds_length": 24,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 15},
                {"exon_number": 2, "length": 0},
                {"exon_number": 3, "length": 9},
            ],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 21,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 21, "codon": "TAA"}],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GCCACATGGCCAAGCTGCTGAAGCTGTAAGCCCCCCCC",
            "alt_transcript_length": 38,
            "alt_cds_start_in_transcript": 5,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 20},
                {"exon_number": 2, "length": 0},
                {"exon_number": 3, "length": 18},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("TTTCA[GCAGCAGCTGCT>]GGTAAGT"),),
        marks=(Span("ref", 0, 20, "20 nt"), Span("ref", 20, 32, "12 nt"), Span("ref", 32, 50, "18 nt")),
    ),
    # EB-26
    Case(
        "delins_over_a_short_exon_that_only_the_left_matching_keeps_both_splice_sites",
        """
        ref 5' [gccac ATG GCC AAG CTG CTG]gtaagtcccccccctttcag[CAG CAG CTG CTG]g--taagtcccccccctttcag[AAG CTG TAA gcccccccc] 3'
        alt 5' [gccac ATG GCC AAG CTG CTG]gtaagtccccccccttttag[AGC AGC TGC TGC]gtctaagtcccccccctttcag[AAG CTG TAA gcccccccc] 3'
                                                           ^^^^^^^^^^^^^^^^^^^^^^^ 16 nt>18 nt
        CAG CAGCAGCTGCTG G > TAG AGCAGCTGCTGC GTC (16 nt > 18 nt): intron -3 to -1, exon 2 and donor +1
        Matched from the left: C>T at intron -3, AG stays, exon 2 becomes AGCAGCTGCTGC, and the donor reads GT
        (G, then the inserted T). Matched from the right: AG stays, but the donor reads CT, so that matching does
        not count.
        alt exon 2: AGCAGCTGCTGC, 12 nt, in frame
        """,
        SHORT_EXON_2,
        Change("CCTTT[CAGCAGCAGCTGCTGG>TAGAGCAGCTGCTGCGTC]TAAGTC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CAGCAGCAGCTGCTGG", "CCAGCAGCTGCTGCTG"),
            "alt": per_strand("TAGAGCAGCTGCTGCGTC", "GACGCAGCAGCTGCTCTA"),
            "variant_start": per_strand(47, 47),
            "variant_end": per_strand(63, 63),
            "alt_cds_seq": "ATGGCCAAGCTGCTGAGCAGCTGCTGCAAGCTGTAA",
            "alt_cds_length": 36,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 15},
                {"exon_number": 2, "length": 12},
                {"exon_number": 3, "length": 9},
            ],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 33,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 33, "codon": "TAA"}],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GCCACATGGCCAAGCTGCTGAGCAGCTGCTGCAAGCTGTAAGCCCCCCCC",
            "alt_transcript_length": 50,
            "alt_cds_start_in_transcript": 5,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 20},
                {"exon_number": 2, "length": 12},
                {"exon_number": 3, "length": 18},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("CCTTT[CAGCAGCAGCTGCTGGT>TAGAGCAGCTGCTGCGTCT]AAGTC"),),
    ),
    # EB-27
    Case(
        "deletion_of_a_whole_intron_between_two_coding_exons_destroys_the_splice_site",
        """
        ref 5' [gccac ATG GCC AAG CTG CTG]gtaagtcccccccctttcag[CAG CAG CTG CTG]|[AAG CTG TAA gcccccccc] 3'
        alt 5' [gccac ATG GCC AAG CTG CTG]--------------------[CAG CAG CTG CTG]|[AAG CTG TAA gcccccccc] 3'
                                          ^^^^^^^^^^^^^^^^^^^^ 20 nt>-
        Exon 1 now ends before CA instead of GT, and exon 2 starts after TG instead of AG. The placement 1 nt to
        the left deletes the last G of exon 1 instead. No placement keeps the splice sites.
        """,
        SHORT_EXON_2,
        Change("GCTGCTG[GTAAGTCCCCCCCCTTTCAG>]CAGCAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("GGTAAGTCCCCCCCCTTTCAG", "GCTGAAAGGGGGGGGACTTAC"),
            "alt": per_strand("G", "G"),
            "variant_start": per_strand(29, 59),
            "variant_end": per_strand(50, 80),
            **UNKNOWN_ALT,
            "unknown_reason": "splice_site_destroyed",
        },
        equivalent=(Change("GCTGCT[GGTAAGTCCCCCCCCTTTCA>]GCAGCAG"),),
    ),
    # EB-28
    Case(
        "snv_in_the_donor_of_a_utr_only_exon_touches_no_coding_region",
        """
        ref 5' [gccaccgcag]gtag[ATG GCC AAG CTG TAA gccgcc] 3'
        alt 5' [gccaccgcag]gcag[ATG GCC AAG CTG TAA gccgcc] 3'
                            ^ T>C
        T>C at donor +2 of exon 1
        Intron 1 has 4 nt: the donor GT of exon 1 and the acceptor AG of exon 2. The SNV lies 3 nt before the
        coding row, but changes only the donor of exon 1, which is an edge of no coding row.
        """,
        UTR_EXON_1_SHORT_INTRON,
        Change("GCAGG[T>C]AGATGG"),
        TOUCHES_NO_CODING_REGION,
    ),
    # EB-29
    Case(
        "snv_in_the_acceptor_before_an_exon_that_begins_with_the_start_codon_destroys_it",
        """
        ref 5' [gccaccgcag]gtaagtcccccccctttcag[ATG GCC ..15.. TAA]|[gccccccccccc] 3'
        alt 5' [gccaccgcag]gtaagtcccccccctttcac[ATG GCC ..15.. TAA]|[gccccccccccc] 3'
                                              ^ G>C
        G>C at acceptor -1
        The start codon begins exon 2, so its edge is also an exon edge with a splice site.
        """,
        CODING_EXON_2,
        Change("CTTTCA[G>C]ATGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(39, 66),
            "variant_end": per_strand(40, 67),
            **UNKNOWN_ALT,
            "unknown_reason": "splice_site_destroyed",
        },
    ),
    # EB-29
    Case(
        "snv_in_the_donor_after_an_exon_that_ends_with_the_stop_codon_destroys_it",
        """
        ref 5' [gccaccgcag]|[ATG ..15.. CTG TAA]gtaagtcccccccctttcag[gccccccccccc] 3'
        alt 5' [gccaccgcag]|[ATG ..15.. CTG TAA]ataagtcccccccctttcag[gccccccccccc] 3'
                                                ^ G>A
        G>A at donor +1
        The stop codon ends exon 2, so its edge is also an exon edge with a splice site.
        """,
        CODING_EXON_2,
        Change("GCTGTAA[G>A]TAAGTCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "variant_start": per_strand(64, 41),
            "variant_end": per_strand(65, 42),
            **UNKNOWN_ALT,
            "unknown_reason": "splice_site_destroyed",
        },
    ),
    # EB-30
    Case(
        "snv_at_intron_plus_5_touches_no_coding_region",
        """
        ref 5' [gccac ATG GCC AAG CTG CTG]gtaagtcccccccctttcag[CAG CAG CTG CTG]|[AAG CTG TAA gcccccccc] 3'
        alt 5' [gccac ATG GCC AAG CTG CTG]gtaaatcccccccctttcag[CAG CAG CTG CTG]|[AAG CTG TAA gcccccccc] 3'
                                              ^ G>A
        G>A at intron +5
        """,
        SHORT_EXON_2,
        Change("AAGCTGCTGGTAA[G>A]TCCCC"),
        TOUCHES_NO_CODING_REGION,
    ),
    # EB-30
    Case(
        "deletion_of_an_a_at_intron_plus_3_that_cannot_shift_into_the_donor_touches_no_coding_region",
        """
        ref 5' [gccac ATG GCC AAG CTG CTG]gtaagtcccccccctttcag[CAG CAG CTG CTG]|[AAG CTG TAA gcccccccc] 3'
        alt 5' [gccac ATG GCC AAG CTG CTG]gt-agtcccccccctttcag[CAG CAG CTG CTG]|[AAG CTG TAA gcccccccc] 3'
                                            ^ A>-
        A>- at intron +3 or +4
        The deletion shifts only within AA. It lies 3 nt after the coding row, but no placement changes the
        donor GT.
        """,
        SHORT_EXON_2,
        Change("AAGCTGCTGGT[A>]AGTCC"),
        TOUCHES_NO_CODING_REGION,
        equivalent=(Change("AAGCTGCTGGTA[A>]GTCC"),),
    ),
    # EB-31
    Case(
        "insertion_right_after_the_stop_codon_inside_an_exon_touches_no_coding_region",
        """
        ref 5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAA -gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAA cgccgcc] 3'
                                                     ^ ->C
        No placement of the C reaches the stop codon TAA, so the C goes into the 3'UTR.
        """,
        THREE_EXONS,
        Change("TCCTAA[>C]GCCGCC"),
        TOUCHES_NO_CODING_REGION,
        equivalent=(Change("TCCTAA[G>CG]CCGCC"),),
    ),
    # EB-32
    Case(
        "stop_codon_inserted_right_before_the_stop_codon_goes_into_the_3utr",
        """
        ref    5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC ---TAA gccgcc] 3'
        alt    5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAATAA gccgcc] 3'
                                                    ^^^ ->TAA
        alt tx                                      22 25
        alt TCC TAA TAA gccgcc. The 5'-most placement lies before the stop codon, the 3'-most after it. The
        stop codon edge takes the 3'-most, so the CDS is unchanged and the 3'UTR starts with taa.
        alt [gacc ATGGCC]|[AAGCTGGGC]|[TCC TAA taagccgcc]: first stop codon at tx 22 = annotated stop codon,
        distance 0
        """,
        THREE_EXONS,
        Change("TCC[>TAA]TAAGCCGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "A"),
            "alt": per_strand("CTAA", "ATTA"),
            "variant_start": per_strand(71, 18),
            "variant_end": per_strand(72, 19),
            "alt_cds_seq": "ATGGCCAAGCTGGGCTCCTAA",
            "alt_cds_length": 21,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 6},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 6},
            ],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 18,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 18, "codon": "TAA"}],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATGGCCAAGCTGGGCTCCTAATAAGCCGCC",
            "alt_transcript_length": 34,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 15},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("TCCTAA[>TAA]GCCGCC"),),
        marks=(Ruler((22, 25), "tx", "alt"),),
    ),
    # EB-33
    Case(
        "deletion_of_one_a_in_the_run_at_the_stop_codon_shortens_the_3utr",
        """
        ref 5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAA aaacctagcc] 3'
        alt 5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC T-A aaacctagcc] 3'
                                                  ^ A>-
        one A deleted anywhere in the run AAaaa
        The 3'-most placement deletes the last a of the run, so the CDS keeps TAA and the 3'UTR loses an a.
        """,
        UTR3_A_RUN,
        Change("CAGTCCT[A>]AAAACC"),
        {
            "variant_id": "var1",
            "ref": per_strand("TA", "TT"),
            "alt": per_strand("T", "T"),
            "variant_start": per_strand(72, 20),
            "variant_end": per_strand(74, 22),
            "alt_cds_seq": "ATGGCCAAGCTGGGCTCCTAA",
            "alt_cds_length": 21,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 6},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 6},
            ],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 18,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 18, "codon": "TAA"}],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATGGCCAAGCTGGGCTCCTAAAACCTAGCC",
            "alt_transcript_length": 34,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 15},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("TCCTAA[A>]AACCTAG"), Change("TCCTAAAA[A>]CCTAG")),
    ),
    # EB-34
    Case(
        "insertion_right_after_a_stop_codon_that_ends_an_exon_touches_no_coding_region",
        """
        ref 5' [gccaccgcag]|[ATG ..15.. CTG TAA--]gtaagtcccccccctttcag[gccccccccccc] 3'
        alt 5' [gccaccgcag]|[ATG ..15.. CTG TAACC]gtaagtcccccccctttcag[gccccccccccc] 3'
                                               ^^ ->CC
        The splice site puts the CC into exon 2, and the stop codon edge puts them into its 3'UTR part.
        """,
        CODING_EXON_2,
        Change("CTGTAA[>CC]GTAAGT"),
        TOUCHES_NO_CODING_REGION,
        equivalent=(Change("CTGTAA[G>CCG]TAAGT"),),
    ),
    # EB-35
    Case(
        "insertion_in_the_a_run_of_a_stop_codon_that_ends_an_exon_goes_into_the_3utr",
        """
        ref    5' [gccaccgcag]|[ATG GCC AAG CTG CTG AAG CTG TA--A]gtaagtcccccccctttcag[gccccccccccc] 3'
        alt    5' [gccaccgcag]|[ATG GCC AAG CTG CTG AAG CTG TAAAA]gtaagtcccccccctttcag[gccccccccccc] 3'
                                                              ^^ ->AA
        alt tx     0            10                          31 34
        AA inserted anywhere in the run AA
        The 3'-most placement inserts AA at the end of exon 2. The splice site keeps them in the exon, and the
        stop codon edge puts them into its 3'UTR part. The CDS is unchanged.
        alt [gccaccgcag]|[ATGGCCAAGCTGCTGAAGCTG TAA aa]|[gccccccccccc]
        """,
        CODING_EXON_2,
        Change("AAGCTGTA[>AA]AGTAAGTCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("AAA", "TTT"),
            "variant_start": per_strand(62, 42),
            "variant_end": per_strand(63, 43),
            "alt_cds_seq": "ATGGCCAAGCTGCTGAAGCTGTAA",
            "alt_cds_length": 24,
            "alt_cds_exons": [{"exon_number": 2, "length": 24}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 21,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 21, "codon": "TAA"}],
            "alt_stop_codon_exons": [2],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GCCACCGCAGATGGCCAAGCTGCTGAAGCTGTAAAAGCCCCCCCCCCC",
            "alt_transcript_length": 48,
            "alt_cds_start_in_transcript": 10,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 26},
                {"exon_number": 3, "length": 12},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("AAGCTGT[>AA]AAGTAAGTCC"), Change("AAGCTGTAA[>AA]GTAAGTCC")),
        marks=(Ruler((0, 10, 31, 34), "tx", "alt"),),
    ),
    # EB-36
    Case(
        "insertion_right_before_the_start_codon_whose_bases_do_not_start_with_atg_touches_no_coding_region",
        """
        ref 5' [gacc -ATG GCC]|[AAG CTG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc cATG GCC]|[AAG CTG GGC]|[TCC TAA gccgcc] 3'
                     ^ ->C
        No alt position of the start codon edge starts an ATG, except the one after the C. So the C goes into
        the 5'UTR.
        """,
        THREE_EXONS,
        Change("GACC[>C]ATGGCC"),
        TOUCHES_NO_CODING_REGION,
        equivalent=(Change("GA[>C]CCATGGCC"),),
    ),
    # EB-37
    Case(
        "insertion_right_before_a_start_codon_that_begins_an_exon_touches_no_coding_region",
        """
        ref 5' [gccaccgcag]gtaagtcccccccctttcag[--ATG GCC AAG ..12.. TAA]|[gccccccccccc] 3'
        alt 5' [gccaccgcag]gtaagtcccccccctttcag[CCATG GCC AAG ..12.. TAA]|[gccccccccccc] 3'
                                                ^^ ->CC
        The splice site puts the CC into exon 2, and the start codon edge puts them into its 5'UTR part.
        """,
        CODING_EXON_2,
        Change("TTTCAG[>CC]ATGGCC"),
        TOUCHES_NO_CODING_REGION,
        equivalent=(Change("TTTCAG[A>CCA]TGGCC"),),
    ),
    # EB-38
    Case(
        "atg_inserted_at_the_start_codon_starts_the_cds_at_the_5_most_atg",
        """
        ref 5' [gacc ---ATG GCC]|[AAG CTG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc atgATG GCC]|[AAG CTG GGC]|[TCC TAA gccgcc] 3'
                     ^^^ ->ATG
        alt gacc ATG ATG GCC. The 5'-most placement puts an ATG at the start codon edge, so the CDS gains a Met
        and keeps its start codon.
        """,
        THREE_EXONS,
        Change("GACC[>ATG]ATGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "T"),
            "alt": per_strand("CATG", "TCAT"),
            "variant_start": per_strand(13, 76),
            "variant_end": per_strand(14, 77),
            "alt_cds_seq": "ATGATGGCCAAGCTGGGCTCCTAA",
            "alt_cds_length": 24,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 9},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 6},
            ],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 21,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 21, "codon": "TAA"}],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATGATGGCCAAGCTGGGCTCCTAAGCCGCC",
            "alt_transcript_length": 34,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 13},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 12},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GACCATG[>ATG]GCCGTAAG"),),
    ),
    # EB-39
    Case(
        "atg_inserted_at_a_start_codon_that_begins_an_exon_starts_the_cds_at_the_5_most_atg",
        """
        ref 5' [gccaccgcag]gtaagtcccccccctttcag[---ATG GCC AAG CTG CTG AAG CTG TAA]|[gccccccccccc] 3'
        alt 5' [gccaccgcag]gtaagtcccccccctttcag[ATGATG GCC AAG CTG CTG AAG CTG TAA]|[gccccccccccc] 3'
                                                ^^^ ->ATG
        alt CAG|ATG ATG GCC. One placement inserts GAT into the acceptor AG and does not count. Of the others,
        the 5'-most puts an ATG at the start of exon 2, so the CDS gains a Met.
        """,
        CODING_EXON_2,
        Change("TTTCAG[>ATG]ATGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "T"),
            "alt": per_strand("GATG", "TCAT"),
            "variant_start": per_strand(39, 65),
            "variant_end": per_strand(40, 66),
            "alt_cds_seq": "ATGATGGCCAAGCTGCTGAAGCTGTAA",
            "alt_cds_length": 27,
            "alt_cds_exons": [{"exon_number": 2, "length": 27}],
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
            "alt_transcript_seq": "GCCACCGCAGATGATGGCCAAGCTGCTGAAGCTGTAAGCCCCCCCCCCC",
            "alt_transcript_length": 49,
            "alt_cds_start_in_transcript": 10,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 27},
                {"exon_number": 3, "length": 12},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("TTTCA[>GAT]GATGGCC"), Change("CAGATG[>ATG]GCCAAG")),
    ),
    # EB-40
    Case(
        "deletion_of_the_a_before_the_start_codon_shortens_the_5utr",
        """
        tx         0      6               18
        ref    5' [ggacca ATG CTG CTG TAA ggccggccgg] 3'
        alt    5' [ggacc- ATG CTG CTG TAA ggccggccgg] 3'
                        ^ A>-
        alt tx     0      5
        one A deleted in the run aA
        Only the placement in the 5'UTR keeps an ATG at the start codon edge. The CDS is unchanged.
        alt [ggacc ATGCTGCTGTAA ggccggccgg]: alt_cds_start_in_transcript = 5, cds_start_in_transcript = 6
        """,
        UTR5_ENDS_IN_A,
        Change("GGACC[A>]ATGCTG"),
        {
            "variant_id": "var1",
            "ref": per_strand("CA", "TT"),
            "alt": per_strand("C", "T"),
            "variant_start": per_strand(14, 31),
            "variant_end": per_strand(16, 33),
            "alt_cds_seq": "ATGCTGCTGTAA",
            "alt_cds_length": 12,
            "alt_cds_exons": [{"exon_number": 1, "length": 12}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 9,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 9, "codon": "TAA"}],
            "alt_stop_codon_exons": [1],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GGACCATGCTGCTGTAAGGCCGGCCGG",
            "alt_transcript_length": 27,
            "alt_cds_start_in_transcript": 5,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [{"exon_number": 1, "length": 27}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GGACCA[A>]TGCTG"),),
        marks=(Ruler((0, 5), "tx", "alt"),),
        ruler=Ruler((0, 6, 18)),
    ),
    # EB-41
    Case(
        "deletion_of_ta_over_the_start_codon_edge_without_atg_takes_the_placement_farthest_into_the_5utr",
        """
        ref 5' [ggacct ATG CTG CTG TAA ggccggccgg] 3'
        alt 5' [ggacc- -TG CTG CTG TAA ggccggccgg] 3'
                     ^^^ TA>-
        tA>- (or AT>-)
        Neither placement leaves an ATG at the start codon edge. The edge takes the placement farthest into the
        5'UTR, which deletes t and A: alt CDS TGCTGCTGTAA, a start loss.
        alt [ggacc TGCTGCTGTAA ggccggccgg] has no ATG, so the scan finds no ORF: alt_has_ptc and stop_loss
        False.
        """,
        UTR5_ENDS_IN_T,
        Change("GGACC[TA>]TGCTG"),
        {
            "variant_id": "var1",
            "ref": per_strand("CTA", "ATA"),
            "alt": per_strand("C", "A"),
            "variant_start": per_strand(14, 30),
            "variant_end": per_strand(17, 33),
            "alt_cds_seq": "TGCTGCTGTAA",
            "alt_cds_length": 11,
            "alt_cds_exons": [{"exon_number": 1, "length": 11}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GGACCTGCTGCTGTAAGGCCGGCCGG",
            "alt_transcript_length": 26,
            "alt_cds_start_in_transcript": 5,
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
            "alt_transcript_exons": [{"exon_number": 1, "length": 26}],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GGACCT[AT>]GCTG"),),
    ),
    # Pins the fallback of the start codon edge rule for an insertion that reaches into the CDS ("Variants at exon
    # boundaries"): "Without such a position, it takes the placement shifted farthest into the 5' UTR. So an
    # insertion right before the start codon changes the CDS only if its bases start with ATG." The start codon is
    # the non-ATG codon CTG ("It can be a non-ATG codon such as CTG"). With an ATG start codon, the placement of
    # such an insertion in the 5'UTR keeps the ATG at the edge, so the fallback is never reached. In EB-36, no
    # placement of the insertion reaches into the CDS, so it gives no row. EB-41 reaches the fallback with a
    # deletion.
    Case(
        "insertion_of_c_before_a_ctg_start_codon_that_can_shift_into_the_cds_goes_into_the_5utr",
        """
        tx         0     4                   19 22
        ref    5' [gcca -CTG GCC AAG CTG TAA gcc] 3'
        alt    5' [gcca cCTG GCC AAG CTG TAA gcc] 3'
                        ^ ->C
        alt tx     0     5
        ->C between the 5'UTR base a at tx 3 and the start codon CTG at tx 4. The C also fits after the C of CTG,
        inside the CDS, so the variant gives a row. It cannot shift farther 5', because tx 3 is a.
        No placement puts an ATG at the start codon edge: the edge reads CTG after the inserted C, and CCT at it.
        So the edge takes the placement farthest into the 5'UTR, and the C goes into the 5'UTR.
        alt [gccac CTGGCCAAGCTGTAA gcc]: the CDS is unchanged, alt_cds_start_in_transcript = 5, start_loss False,
        annotated_stop_distance = 0.
        """,
        CTG_START_AFTER_A,
        Change("gcca[>C]CTGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "G"),
            "alt": per_strand("AC", "GG"),
            "variant_start": per_strand(13, 27),
            "variant_end": per_strand(14, 28),
            "alt_cds_seq": "CTGGCCAAGCTGTAA",
            "alt_cds_length": 15,
            "alt_cds_exons": [{"exon_number": 1, "length": 15}],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 12,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 12, "codon": "TAA"}],
            "alt_stop_codon_exons": [1],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GCCACCTGGCCAAGCTGTAAGCC",
            "alt_transcript_length": 23,
            "alt_cds_start_in_transcript": 5,
            "alt_transcript_exons": [{"exon_number": 1, "length": 23}],
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("gccaC[>C]TGGCC"), Change("gcc[aC>aCC]TGGCC")),
        marks=(Ruler((0, 5), "tx", "alt"),),
        ruler=Ruler((0, 4, 19, 22)),
    ),
    # EB-42
    Case(
        "delins_ca_to_gatg_at_the_start_codon_edge_starts_the_cds_at_the_atg_of_the_left_matching",
        """
        tx      0    4
        ref 5' [gacc A--TG GCC]|[AAG CTG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacg ATGTG GCC]|[AAG CTG GGC]|[TCC TAA gccgcc] 3'
                   ^^^^^ CA>GATG
        cA>GATG: the last 5'UTR base and the A of the start codon
        Matched from the left, c>G stays in the 5'UTR and ATG starts at the start codon edge. Matched from the
        right, no ATG starts there.
        alt [gacg ATGTGGCC]|[AAGCTGGGC]|[TCCTAAgccgcc]: alt_cds_start_in_transcript = 4. The CDS gains 2 nt, and
        the shifted frame reads no stop codon up to the transcript end (nonstop).
        """,
        THREE_EXONS,
        Change("GAC[CA>GATG]TGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CA", "TG"),
            "alt": per_strand("GATG", "CATC"),
            "variant_start": per_strand(13, 76),
            "variant_end": per_strand(15, 78),
            "alt_cds_seq": "ATGTGGCCAAGCTGGGCTCCTAA",
            "alt_cds_length": 23,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 8},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 6},
            ],
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
            "alt_transcript_seq": "GACGATGTGGCCAAGCTGGGCTCCTAAGCCGCC",
            "alt_transcript_length": 33,
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
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 12},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 12},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GA[CCA>CGATG]TGGCC"), Change("GAC[CAT>GATGT]GGCC")),
        ruler=Ruler((0, 4)),
    ),
    # Pins the coding edge at the start of REF ("Variants at exon boundaries"): "An edge at an end of REF maps to the
    # same end of ALT, also for a delins whose REF and ALT differ in length." and "A>GGC at the first base of the start
    # codon ATG gives one that starts with GGCTG, and `start_loss` is True." Neither placement puts an ATG at the
    # edge, and the alt transcript has no ATG: "Both flags are then False, and `annotated_stop_distance` is null." EB-42
    # has the start codon edge strictly inside REF, where the matchings map it base for base.
    Case(
        "delins_a_to_ggc_at_the_first_base_of_the_start_codon_keeps_all_alt_bases_in_the_cds",
        """
        tx      0    4
        ref 5' [gacc A--TG GCC]|[AAG CTG GGC]|[TCC TAA gccgcc] 3'
        alt 5' [gacc GGCTG GCC]|[AAG CTG GGC]|[TCC TAA gccgcc] 3'
                     ^^^ A>GGC
        A>GGC at the first base of the start codon. REF and ALT share no first or last base, so the delins has two
        matchings. The start codon edge lies at the start of REF, so both matchings map it to the start of ALT, and
        GGC is coding. No ATG starts at the edge, and the alt transcript has no ATG.
        alt [gacc GGCTGGCC]|[AAGCTGGGC]|[TCCTAAgccgcc]: a start loss without an ATG, so both flags are False and
        annotated_stop_distance is null.
        """,
        THREE_EXONS,
        Change("GACC[A>GGC]TGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("GGC", "GCC"),
            "variant_start": per_strand(14, 76),
            "variant_end": per_strand(15, 77),
            "alt_cds_seq": "GGCTGGCCAAGCTGGGCTCCTAA",
            "alt_cds_length": 23,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 8},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 6},
            ],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "GACCGGCTGGCCAAGCTGGGCTCCTAAGCCGCC",
            "alt_transcript_length": 33,
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
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 12},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 12},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GACC[AT>GGCT]GGCC"),),
        ruler=Ruler((0, 4)),
    ),
    # EB-43
    Case(
        "delins_aa_to_ccc_over_the_stop_codon_end_puts_the_length_change_into_the_3utr",
        """
        tx         0    4                           22  25
        ref    5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAA a-aacctagcc] 3'
        alt    5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAC ccaacctagcc] 3'
                                                      ^^^^ AA>CCC
        alt tx                                      22  25 28 31
        Aa>CCC: the last stop codon base and the first 3'UTR base
        The matching with its length change at the 3'UTR end counts: C replaces the last stop codon base, and
        CC goes into the 3'UTR. The stop codon becomes TAC and is lost.
        alt [...TCC TAC cca acc tag cc]: the next in-frame stop codon is the tag at tx 31, distance 22 - 31 = -9
        """,
        UTR3_A_RUN,
        Change("TCCTA[AA>CCC]AACCTAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("AA", "TT"),
            "alt": per_strand("CCC", "GGG"),
            "variant_start": per_strand(74, 19),
            "variant_end": per_strand(76, 21),
            "alt_cds_seq": "ATGGCCAAGCTGGGCTCCTAC",
            "alt_cds_length": 21,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 6},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 6},
            ],
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
            "alt_transcript_seq": "GACCATGGCCAAGCTGGGCTCCTACCCAACCTAGCC",
            "alt_transcript_length": 36,
            "alt_cds_start_in_transcript": 4,
            "alt_scan_start_codon_pos": 4,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAG",
            "alt_scan_first_stop_pos": 31,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 31, "codon": "TAG"}],
            "alt_scan_stop_codon_exons": [3],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -9,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 17},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("TCCT[AAA>ACCC]AACCTAG"),),
        marks=(Ruler((22, 25, 28, 31), "tx", "alt"),),
        ruler=Ruler((0, 4, 22, 25)),
    ),
    # Pins the coding edge at the end of REF ("Variants at exon boundaries"): "An edge at an end of REF maps to the
    # same end of ALT, also for a delins whose REF and ALT differ in length." and "Both placements map it to the same
    # end of ALT, so the extra ALT bases stay coding. E.g. A>CG at the last base of the stop codon TAA gives an alt CDS
    # that ends with TACG." The stop codon is lost, and the alt transcript has no in-frame stop codon: a nonstop. EB-43
    # has the stop codon edge strictly inside REF, where the matchings map it base for base.
    Case(
        "delins_a_to_cg_at_the_last_base_of_the_stop_codon_keeps_all_alt_bases_in_the_cds",
        """
        tx      0    4                           22
        ref 5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAA- gccgcc] 3'
        alt 5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TACg gccgcc] 3'
                                                   ^^ A>CG
        A>CG at the last base of the stop codon. REF and ALT share no first or last base, so the delins has two
        matchings. The stop codon edge lies at the end of REF, so both matchings map it to the end of ALT, and C and
        G are coding. The alt line aligns the G with the 3'UTR, but the alt CDS ends with TACG.
        alt [gaccATGGCC]|[AAGCTGGGC]|[TCCTACG gccgcc]: the stop codon is lost, and the alt transcript reads no
        in-frame stop codon from the start codon to its end (nonstop).
        """,
        THREE_EXONS,
        Change("TCCTA[A>CG]GCCGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("CG", "CG"),
            "variant_start": per_strand(74, 16),
            "variant_end": per_strand(75, 17),
            "alt_cds_seq": "ATGGCCAAGCTGGGCTCCTACG",
            "alt_cds_length": 22,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 6},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 7},
            ],
            "alt_last_codon": "ACG",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": "GACCATGGCCAAGCTGGGCTCCTACGGCCGCC",
            "alt_transcript_length": 32,
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
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 13},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("TCCT[AA>ACG]GCCGCC"),),
        ruler=Ruler((0, 4, 22)),
    ),
    # EB-44
    Case(
        "delins_over_the_5utr_the_start_codon_and_a_donor_takes_the_start_codon_edge_from_the_matching_that_keeps_the_donor",
        """
        ref 5' cccccccccc[gacc ATG GCC]gtaagt----------------cccccccctttcag[AAG CTG GGC]|[TCC TAA gccgcc] 3'
        alt 5' cccccccccc[cccc CCC CCC]gtcccccccccccccccccccccccccccctttcag[AAG CTG GGC]|[TCC TAA gccgcc] 3'
                          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 16 nt>32 nt
        exon 1 and intron +1 to +6 > 10 C, GT, 20 C
        Matched from the left, GT stays after exon 1, which becomes 10 C. Matched from the right, the donor reads
        CC, so that matching does not count, although it lies farther into the 5'UTR. No placement puts an ATG
        at the start codon edge, so the left matching decides: the CDS starts at tx 4 with CCCCCC, a start loss.
        alt [cccc CCCCCC]|[AAGCTGGGC]|[TCCTAAgccgcc] has no ATG, so the scan finds no ORF: alt_has_ptc and
        stop_loss False.
        """,
        THREE_EXONS,
        Change("CCCCC[GACCATGGCCGTAAGT>" + "C" * 10 + "GT" + "C" * 20 + "]CCCCCCCCTTTCAGAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("GACCATGGCCGTAAGT", "ACTTACGGCCATGGTC"),
            "alt": per_strand("C" * 10 + "GT" + "C" * 20, "G" * 20 + "AC" + "G" * 10),
            "variant_start": per_strand(10, 65),
            "variant_end": per_strand(26, 81),
            "alt_cds_seq": "CCCCCCAAGCTGGGCTCCTAA",
            "alt_cds_length": 21,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 6},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 6},
            ],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 18,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 18, "codon": "TAA"}],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "CCCCCCCCCCAAGCTGGGCTCCTAAGCCGCC",
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
        equivalent=(Change("CCCC[CGACCATGGCCGTAAGT>" + "C" * 11 + "GT" + "C" * 20 + "]CCCCCCCCTTTCAGAAG"),),
    ),
    # EB-45
    Case(
        "delins_over_the_stop_codon_and_the_donor_keeps_the_donor_only_in_the_right_matching",
        """
        ref 5' [gccaccgcag]|[ATG ..12.. AAG CTG TAA]gta--agtcccccccctttcag[gccccccccccc] 3'
        alt 5' [gccaccgcag]|[ATG ..12.. AAG CTG TGT]aagttagtcccccccctttcag[gccccccccccc] 3'
                                                 ^^^^^^^^ AAGTA>GTAAGTT
        Matched from the left, the donor reads AA. Matched from the right, GT follows the alt GTAA, so exon 2
        ends with CTG TGT AA, 2 nt longer. The stop codon ends the exon, so the CDS holds every alt exon base.
        alt CDS ATG...CTG TGT AA: the shifted frame reads no stop codon up to the transcript end (nonstop)
        """,
        CODING_EXON_2,
        Change("GCTGT[AAGTA>GTAAGTT]AGTCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("AAGTA", "TACTT"),
            "alt": per_strand("GTAAGTT", "AACTTAC"),
            "variant_start": per_strand(62, 39),
            "variant_end": per_strand(67, 44),
            "alt_cds_seq": "ATGGCCAAGCTGCTGAAGCTGTGTAA",
            "alt_cds_length": 26,
            "alt_cds_exons": [{"exon_number": 2, "length": 26}],
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
            "alt_transcript_seq": "GCCACCGCAGATGGCCAAGCTGCTGAAGCTGTGTAAGCCCCCCCCCCC",
            "alt_transcript_length": 48,
            "alt_cds_start_in_transcript": 10,
            "alt_scan_start_codon_pos": 10,
            "alt_scan_start_codon_exon": 2,
            "alt_scan_first_stop_codon": None,
            "alt_scan_first_stop_pos": None,
            "alt_scan_stop_codon_count": 0,
            "alt_scan_stop_codons": [],
            "alt_scan_stop_codon_exons": [],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": None,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 26},
                {"exon_number": 3, "length": 12},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("AAGCTG[TAAGTA>TGTAAGTT]AGTCC"),),
    ),
    # EB-46
    Case(
        "deletion_from_the_stop_codon_into_the_3utr_reads_on_into_the_3utr",
        """
        tx         0    4                           22  25
        ref    5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAA cctagcctgacc] 3'
        alt    5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC T-- --tagcctgacc] 3'
                                                     ^^^^^ AACC>-
        alt tx                                      22      25 28
        The CDS loses AA and ends in TCC T, and the 3'UTR loses cc.
        alt [...GGC TCC Tta gcc tga cc]: the next in-frame stop codon is the tga at tx 28, distance 22 - 28 = -6
        """,
        UTR3_CCTAG,
        Change("CAGTCCT[AACC>]TAGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("TAACC", "AGGTT"),
            "alt": per_strand("T", "A"),
            "variant_start": per_strand(72, 19),
            "variant_end": per_strand(77, 24),
            "alt_cds_seq": "ATGGCCAAGCTGGGCTCCT",
            "alt_cds_length": 19,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 6},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 4},
            ],
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
            "alt_transcript_seq": "GACCATGGCCAAGCTGGGCTCCTTAGCCTGACC",
            "alt_transcript_length": 33,
            "alt_cds_start_in_transcript": 4,
            "alt_scan_start_codon_pos": 4,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TGA",
            "alt_scan_first_stop_pos": 28,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 28, "codon": "TGA"}],
            "alt_scan_stop_codon_exons": [3],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -6,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 14},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("CAGTCCT[AACCT>T]AGCC"),),
        marks=(Ruler((22, 25, 28), "tx", "alt"),),
        ruler=Ruler((0, 4, 22, 25)),
    ),
    # EB-47
    Case(
        "deletion_across_the_stop_codon_lands_a_tag_of_the_3utr_at_the_stop_codon_position",
        """
        tx         0    4                           22  25
        ref    5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAA cctagcctgacc] 3'
        alt    5' [gacc ATG GCC]|[AAG CTG GGC]|[T-- --- cctagcctgacc] 3'
                                                 ^^^^^^ CCTAA>-
        alt tx                                            22
        5 nt deleted: CCTAA, CTAAc, TAAcc, AActa or Accta
        Each placement gives TCC TAG cctgacc. The 3'-most one (Accta) keeps TA of the stop codon, so the stop
        codon edge takes it and the CDS ends in TCC TA. The TAG at tx 22 is the annotated stop codon: neither a
        PTC nor a stop loss, distance 0.
        """,
        UTR3_CCTAG,
        Change("CAGT[CCTAA>]CCTAGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("TCCTAA", "GTTAGG"),
            "alt": per_strand("T", "G"),
            "variant_start": per_strand(69, 21),
            "variant_end": per_strand(75, 27),
            "alt_cds_seq": "ATGGCCAAGCTGGGCTCCTA",
            "alt_cds_length": 20,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 6},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 5},
            ],
            "alt_last_codon": "CTA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "GACCATGGCCAAGCTGGGCTCCTAGCCTGACC",
            "alt_transcript_length": 32,
            "alt_cds_start_in_transcript": 4,
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 13},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("CAGTCC[TAACC>]TAGCCTG"), Change("CAGTCCTA[ACCTA>]GCCTGACC")),
        marks=(Ruler((22,), "tx", "alt"),),
        ruler=Ruler((0, 4, 22, 25)),
    ),
    # EB-48
    Case(
        "delins_from_the_stop_codon_into_an_a_run_of_the_3utr_puts_tga_right_after_the_cds",
        """
        tx         0    4                           22  25
        ref    5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAA aa--acctagcc] 3'
        alt    5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC CCC tgacacctagcc] 3'
                                                    ^^^^^^^^ TAAAA>CCCTGAC
        alt tx                                      22  25
        TAAaa>CCCTGAC: the stop codon and the first 2 3'UTR bases
        Matched from the left, the length change lies at the 3'UTR end: CCC replaces the stop codon, and the
        3'UTR starts with TGAC.
        alt [...TCC CCC tga cac ctagcc]: the next in-frame stop codon is the TGA at tx 25, distance 22 - 25 = -3
        """,
        UTR3_A_RUN,
        Change("CAGTCC[TAAAA>CCCTGAC]ACCTAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("TAAAA", "TTTTA"),
            "alt": per_strand("CCCTGAC", "GTCAGGG"),
            "variant_start": per_strand(72, 18),
            "variant_end": per_strand(77, 23),
            "alt_cds_seq": "ATGGCCAAGCTGGGCTCCCCC",
            "alt_cds_length": 21,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 6},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 6},
            ],
            "alt_last_codon": "CCC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": "GACCATGGCCAAGCTGGGCTCCCCCTGACACCTAGCC",
            "alt_transcript_length": 37,
            "alt_cds_start_in_transcript": 4,
            "alt_scan_start_codon_pos": 4,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TGA",
            "alt_scan_first_stop_pos": 25,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 25, "codon": "TGA"}],
            "alt_scan_stop_codon_exons": [3],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -3,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 18},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("CAGTCC[TAAAAA>CCCTGACA]CCTAG"),),
        marks=(Ruler((22, 25), "tx", "alt"),),
        ruler=Ruler((0, 4, 22, 25)),
    ),
    # EB-48
    Case(
        "delins_from_the_stop_codon_into_the_3utr_puts_tga_right_after_the_cds",
        """
        tx         0    4                           22  25
        ref    5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAA g--ccgcc] 3'
        alt    5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC CCC tgaccgcc] 3'
                                                    ^^^^^^^ TAAG>CCCTGA
        alt tx                                      22  25
        TAAg>CCCTGA: the stop codon and the first 3'UTR base
        Matched from the left, the length change lies at the 3'UTR end: CCC replaces the stop codon, and the
        3'UTR starts with TGA.
        alt [...TCC CCC tga ccg cc]: the next in-frame stop codon is the TGA at tx 25, distance 22 - 25 = -3
        """,
        THREE_EXONS,
        Change("TCC[TAAG>CCCTGA]CCGCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("TAAG", "CTTA"),
            "alt": per_strand("CCCTGA", "TCAGGG"),
            "variant_start": per_strand(72, 15),
            "variant_end": per_strand(76, 19),
            "alt_cds_seq": "ATGGCCAAGCTGGGCTCCCCC",
            "alt_cds_length": 21,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 6},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 6},
            ],
            "alt_last_codon": "CCC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": True,
            "alt_transcript_seq": "GACCATGGCCAAGCTGGGCTCCCCCTGACCGCC",
            "alt_transcript_length": 33,
            "alt_cds_start_in_transcript": 4,
            "alt_scan_start_codon_pos": 4,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TGA",
            "alt_scan_first_stop_pos": 25,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [{"position": 25, "codon": "TGA"}],
            "alt_scan_stop_codon_exons": [3],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -3,
            **NO_RULE,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 14},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("TC[CTAAG>CCCCTGA]CCGCCC"),),
        marks=(Ruler((22, 25), "tx", "alt"),),
        ruler=Ruler((0, 4, 22, 25)),
    ),
    # EB-49
    Case(
        "deletion_from_the_stop_codon_past_the_transcript_end_shortens_the_last_exon",
        """
        tx      0    4                           22  25    31
        ref 5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC TAA gccgcc]cccccccccc 3'
        alt 5' [gacc ATG GCC]|[AAG CTG GGC]|[TCC T-- ------]-ccccccccc 3'
                                                  ^^^^^^^^^^^ AAGCCGCCC>-
        AAgccgcc and 1 base after the transcript end deleted
        The transcript end has no splice site, so the deletion is known. The last exon ends after TCC T.
        alt [gaccATGGCC]|[AAGCTGGGC]|[TCCT]: no stop codon up to the transcript end (nonstop)
        """,
        THREE_EXONS,
        Change("CAGTCCT[AAGCCGCCC>]CCCCCCCCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("TAAGCCGCCC", "GGGGCGGCTT"),
            "alt": per_strand("T", "G"),
            "variant_start": per_strand(72, 8),
            "variant_end": per_strand(82, 18),
            "alt_cds_seq": "ATGGCCAAGCTGGGCTCCT",
            "alt_cds_length": 19,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 6},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 4},
            ],
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
            "alt_transcript_seq": "GACCATGGCCAAGCTGGGCTCCT",
            "alt_transcript_length": 23,
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
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 9},
                {"exon_number": 3, "length": 4},
            ],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("CAGTCCT[AAGCCGCCCC>C]CCCCCCCC"),),
        ruler=Ruler((0, 4, 22, 25, 31)),
    ),
]
