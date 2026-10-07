"""
Conformance cases of variant placement and exon boundaries. They cover SNVs in and next to the splice
dinucleotides, the reach of a variant into a splice site, indels and delins at an exon edge, the transcript start,
and the rows with unknown_reason ("Technical Notes.md", section "Variants at exon boundaries").

The drawings show the transcript 5' to 3', also for the minus strand. The runner renders their layout block from the
case. The line ref is the transcript and the line alt is the transcript with the change. The two lines are aligned,
and `-` fills a gap. `[...]` is an exon and `|` is an exon junction. Upper case is the coding region (the CDS with
the stop codon), in codons of the annotated frame; lower case is UTR. An intron or a flank that the change touches is
drawn in lower case outside the brackets. `..N..` leaves out N bases. `^` marks the change as ref>alt. A ruler gives
tx or layout positions, as labelled. The layout position counts from the first base of the 5' flank (10 nt).
"""

from .runner import (
    IDS,
    INTRON,
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
    Transcript,
    per_strand,
)

INTRON_1 = "GTGTAAGCCCCCCTTTTCAG"
EXON_2 = "AGTGAACGTTGGAAGC"

# layout  10         18                  38                    54                74      83           92
# tx      0                              8                                       24      33           42
# ref 5' [ATG GCT CT]gtgtaag..6..ttttcag[A GTG AAC GTT GGA AGC]gtaagt..8..tttcag[CTG CGT TAA aaagctgcc] 3'
#         exon 1     CT|GTGT donor       exon 2, CAG|AG acceptor                 exon 3
MAIN_REF = {
    **IDS,
    "cds_start": per_strand(10, 19),
    "cds_end": per_strand(83, 92),
    "ref_cds_seq": "ATGGCTCTAGTGAACGTTGGAAGCCTGCGTTAA",
    "ref_cds_length": 33,
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "ref_cds_exons": [
        {"exon_number": 1, "length": 8},
        {"exon_number": 2, "length": 16},
        {"exon_number": 3, "length": 9},
    ],
    "cds_in_transcript": True,
    "start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_first_stop_pos": 30,
    "ref_stop_codon_count": 1,
    "ref_stop_codons": [{"position": 30, "codon": "TAA"}],
    "ref_stop_codon_exons": [3],
    "ref_has_ptc": False,
    "transcript_start": 10,
    "transcript_end": 92,
    "transcript_seq": "ATGGCTCTAGTGAACGTTGGAAGCCTGCGTTAAAAAGCTGCC",
    "transcript_length": 42,
    "cds_start_in_transcript": 0,
    "cds_end_in_transcript": 33,
    "transcript_exons": [
        {"exon_number": 1, "length": 8},
        {"exon_number": 2, "length": 16},
        {"exon_number": 3, "length": 18},
    ],
    "utr3_length": 9,
    "utr5_length": 0,
    "total_exon_count": 3,
    "likely_misannotated": False,
}
MAIN = Layout(Transcript(("ATGGCTCT", EXON_2, "CTGCGTTAAaaagctgcc"), introns=(INTRON_1, INTRON)), MAIN_REF)
# The same transcript, with an A as the last base of the 5' flank: CCCCCCCCCA|[ATGG...
MAIN_AFTER_A = Layout(
    Transcript(("ATGGCTCT", EXON_2, "CTGCGTTAAaaagctgcc"), introns=(INTRON_1, INTRON), flanks=("CCCCCCCCCA", "C" * 10)),
    MAIN_REF,
)

# The alt columns of a known row whose alt CDS keeps 33 nt and its only in-frame stop codon, the annotated one
KEEPS_THE_STOP_CODON = {
    "alt_cds_length": 33,
    "alt_cds_exons": [
        {"exon_number": 1, "length": 8},
        {"exon_number": 2, "length": 16},
        {"exon_number": 3, "length": 9},
    ],
    "alt_last_codon": "TAA",
    "alt_valid_stop": True,
    "alt_first_stop_codon": "TAA",
    "alt_first_stop_pos": 30,
    "alt_stop_codon_count": 1,
    "alt_stop_codons": [{"position": 30, "codon": "TAA"}],
    "alt_stop_codon_exons": [3],
    "alt_has_ptc": False,
    "start_loss": False,
    "stop_loss": False,
    "alt_transcript_length": 42,
    "alt_cds_start_in_transcript": 0,
    **NOT_SCANNED,
    "unknown_reason": None,
    **NO_PTC_FEATURES,
    "annotated_stop_distance": 0,
    **NO_RULE,
}
# The alt columns of a known row after a frameshift in the CDS of MAIN: the alt CDS ends in ...GCG TTA A, out of
# frame, and the alt transcript reads on through ...TTA AAA AGC TGC C without a stop codon (nonstop)
FRAMESHIFT_NONSTOP = {
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
    "alt_cds_start_in_transcript": 0,
    "alt_scan_start_codon_pos": 0,
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
}
DESTROYED = {**UNKNOWN_ALT, "unknown_reason": "splice_site_destroyed"}
AMBIGUOUS = {**UNKNOWN_ALT, "unknown_reason": "exon_boundary_ambiguous"}

# layout  10  13      19 22
# ref 5' [gcc ATG TAA gcc] 3'   one exon, 6 nt coding region
SHORT_CDS = Layout(
    Transcript(("gccATGTAAgcc",)),
    {
        **IDS,
        "cds_start": 13,
        "cds_end": 19,
        "ref_cds_seq": "ATGTAA",
        "ref_cds_length": 6,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [{"exon_number": 1, "length": 6}],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 3,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [{"position": 3, "codon": "TAA"}],
        "ref_stop_codon_exons": [1],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 22,
        "transcript_seq": "GCCATGTAAGCC",
        "transcript_length": 12,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 9,
        "transcript_exons": [{"exon_number": 1, "length": 12}],
        "utr3_length": 3,
        "utr5_length": 3,
        "total_exon_count": 1,
        "likely_misannotated": False,
    },
)


# One exon that ends with the stop codon TAA, 1 nt before the chromosome end: the 3' flank is the single base A. On
# the minus strand, the transcript ends 1 nt after the chromosome start.
STOP_CODON_1_NT_FROM_THE_CHROMOSOME_EDGE = Layout(
    Transcript(("gccATGGCCAAGCTGTAA",), flanks=("C" * 10, "A")),
    {
        **IDS,
        "cds_start": per_strand(13, 1),
        "cds_end": per_strand(28, 16),
        "ref_cds_seq": "ATGGCCAAGCTGTAA",
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
        "transcript_start": per_strand(10, 1),
        "transcript_end": per_strand(28, 19),
        "transcript_seq": "GCCATGGCCAAGCTGTAA",
        "transcript_length": 18,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 18,
        "transcript_exons": [{"exon_number": 1, "length": 18}],
        "utr3_length": 0,
        "utr5_length": 3,
        "total_exon_count": 1,
        "likely_misannotated": False,
    },
)
# Exon 1 is the start codon ATG, 1 nt after the chromosome start: the 5' flank is the single base A. On the minus
# strand, the transcript starts 1 nt before the chromosome end.
START_CODON_1_NT_FROM_THE_CHROMOSOME_EDGE = Layout(
    Transcript(("ATG", "GCCAAGCTGTAAgcc"), flanks=("A", "C" * 10)),
    {
        **IDS,
        "cds_start": per_strand(1, 13),
        "cds_end": per_strand(36, 48),
        "ref_cds_seq": "ATGGCCAAGCTGTAA",
        "ref_cds_length": 15,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [{"exon_number": 1, "length": 3}, {"exon_number": 2, "length": 12}],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 12,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [{"position": 12, "codon": "TAA"}],
        "ref_stop_codon_exons": [2],
        "ref_has_ptc": False,
        "transcript_start": per_strand(1, 10),
        "transcript_end": per_strand(39, 48),
        "transcript_seq": "ATGGCCAAGCTGTAAGCC",
        "transcript_length": 18,
        "cds_start_in_transcript": 0,
        "cds_end_in_transcript": 15,
        "transcript_exons": [{"exon_number": 1, "length": 3}, {"exon_number": 2, "length": 15}],
        "utr3_length": 3,
        "utr5_length": 0,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)


def record(ref, alt, start, end):
    """The columns that echo the VCF record. Each argument is a value or a per_strand pair."""
    return {"variant_id": "var1", "ref": ref, "alt": alt, "start": start, "end": end}


CASES = [
    Case(
        "missense_snv_inside_an_internal_coding_exon_changes_one_cds_base",
        """
        ref 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]|[A GTG AAC GCT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                                        ^ T>C
        T>C at CDS 16: GTT>GCT, 8 nt from the exon start, 7 nt from the exon end
        one placement; the alt CDS is the ref CDS with the SNV, unknown_reason null
        """,
        MAIN,
        Change("AACG[T>C]TGGA"),
        {
            **record(per_strand("T", "A"), per_strand("C", "G"), per_strand(46, 55), per_strand(47, 56)),
            **KEEPS_THE_STOP_CODON,
            "alt_cds_seq": "ATGGCTCTAGTGAACGCTGGAAGCCTGCGTTAA",
            "alt_transcript_seq": "ATGGCTCTAGTGAACGCTGGAAGCCTGCGTTAAAAAGCTGCC",
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "alt_transcript",
        },
        equivalent=(Change("AAC[GT>GC]TGGA"),),
    ),
    Case(
        "snv_at_donor_plus_1_destroys_the_splice_site",
        """
        ref 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtaagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]ataagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
                                                   ^ G>A
        G>A at donor +1 after coding exon 2: GT>AT
        no valid placement: splice_site_destroyed, every alt column null
        """,
        MAIN,
        Change("GGAAGC[G>A]TAAGTCC"),
        {
            **record(per_strand("G", "C"), per_strand("A", "T"), per_strand(54, 47), per_strand(55, 48)),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
    ),
    Case(
        "snv_at_donor_plus_2_destroys_the_splice_site",
        """
        ref 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtaagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gaaagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
                                                    ^ T>A
        T>A at donor +2: GT>GA
        """,
        MAIN,
        Change("GGAAGCG[T>A]AAGTCC"),
        {
            **record(per_strand("T", "A"), per_strand("A", "T"), per_strand(55, 46), per_strand(56, 47)),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
    ),
    Case(
        "snv_at_acceptor_minus_1_destroys_the_splice_site",
        """
        ref 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtaagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtaagtcccccccctttcac[CTG CGT TAA aaagctgcc] 3'
                                                                      ^ G>C
        G>C at acceptor -1 of coding exon 3: AG>AC
        """,
        MAIN,
        Change("CCTTTCA[G>C]CTGCGT"),
        {
            **record(per_strand("G", "C"), per_strand("C", "G"), per_strand(73, 28), per_strand(74, 29)),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
    ),
    Case(
        "snv_at_acceptor_minus_2_destroys_the_splice_site",
        """
        ref 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtaagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtaagtcccccccctttcgg[CTG CGT TAA aaagctgcc] 3'
                                                                     ^ A>G
        A>G at acceptor -2: AG>GG
        """,
        MAIN,
        Change("CCTTTC[A>G]GCTGCGT"),
        {
            **record(per_strand("A", "T"), per_strand("G", "C"), per_strand(72, 29), per_strand(73, 30)),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
    ),
    # The +2 side is snv_at_donor_plus_2_destroys_the_splice_site.
    Case(
        "snv_at_donor_plus_3_touches_no_coding_region",
        """
        ref 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtaagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtcagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
                                                     ^ A>C
        A>C at donor +3, outside the splice dinucleotide: no row
        """,
        MAIN,
        Change("GAAGCGT[A>C]AGTCC"),
        NoRow("touches no coding region"),
    ),
    # An insertion inside the donor dinucleotide reaches the coding row
    Case(
        "insertion_between_donor_plus_1_and_plus_2_destroys_the_splice_site",
        """
        ref 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]g-taagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gataagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
                                                    ^ ->A
        A inserted between donor +1 and +2: GT>GAT
        """,
        MAIN,
        Change("GGAAGCG[>A]TAAGTCC"),
        {
            **record(per_strand("G", "A"), per_strand("GA", "AT"), per_strand(54, 46), per_strand(55, 47)),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("GGAAGCG[T>AT]AAGTCC"),),
    ),
    # An insertion right after the donor dinucleotide changes neither it nor the coding row
    Case(
        "insertion_between_donor_plus_2_and_plus_3_touches_no_coding_region",
        """
        ref 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gt-aagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtcaagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
                                                     ^ ->C
        C inserted between donor +2 and +3: GT stays, no row
        """,
        MAIN,
        Change("GGAAGCGT[>C]AAGTCC"),
        NoRow("touches no coding region"),
    ),
    Case(
        "deletion_from_the_last_exon_base_to_donor_plus_5_destroys_the_splice_site",
        """
        ref 5' [ATG GCT CT]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT C-]-----agccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                         ^^^^^^^ TGTGTA>-
        TGTGTA deleted from exon -1 to donor +5: ...GGCTC|AGCC, no GT after the exon end
        one placement, not valid: splice_site_destroyed
        """,
        MAIN,
        Change("GGCTC[TGTGTA>]AGCCCC"),
        {
            **record(per_strand("CTGTGTA", "TTACACA"), per_strand("C", "T"), per_strand(16, 78), per_strand(23, 85)),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("GGCTC[TGTGTAA>A]GCCCC"),),
    ),
    Case(
        "deletion_of_tg_over_a_donor_with_an_intronic_placement_keeps_the_exon",
        """
        ref 5' [ATG GCT CT]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT C-]-tgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                         ^^^ TG>-
        TG deleted over the exon end: CT|GTGT > C|TGTAAG
        the same deletion as GT at donor +1+2 (or TG, GT further 3'): those placements keep GT after the exon end,
        so exon 1 stays as it is; the alt CDS equals the ref CDS
        """,
        MAIN,
        Change("GGCTC[TG>]TGTAAGC"),
        {
            **record(per_strand("CTG", "ACA"), per_strand("C", "A"), per_strand(16, 82), per_strand(19, 85)),
            **KEEPS_THE_STOP_CODON,
            "alt_cds_seq": "ATGGCTCTAGTGAACGTTGGAAGCCTGCGTTAA",
            "alt_transcript_seq": "ATGGCTCTAGTGAACGTTGGAAGCCTGCGTTAAAAAGCTGCC",
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "alt_transcript",
        },
        equivalent=(Change("GGCTCT[GT>]GTAAGC"), Change("GGCTCTGT[GT>]AAGC")),
    ),
    Case(
        "deletion_of_ag_at_a_cag_ag_acceptor_keeps_the_acceptor_only_if_the_exon_loses_ag",
        """
        ref 5' [ATG GCT CT]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]gtgtaagccccccttttc--[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                                             ^^ AG>-
        AG deleted before exon 2, the same as GA or the exon's AG
        only the placement that deletes the exon's AG keeps AG before the exon: exon 2 loses 2 nt (tx 8, 9)
        frameshift: ATG GCT CTT GAA CGT TGG AAG CCT GCG TTA AAA AGC TGC C, no stop codon (nonstop)
        """,
        MAIN,
        Change("CCTTTTC[AG>]AGTGAAC"),
        {
            **record(per_strand("CAG", "TCT"), per_strand("C", "T"), per_strand(35, 63), per_strand(38, 66)),
            **FRAMESHIFT_NONSTOP,
            "alt_cds_seq": "ATGGCTCTTGAACGTTGGAAGCCTGCGTTAA",
            "alt_cds_length": 31,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 8},
                {"exon_number": 2, "length": 14},
                {"exon_number": 3, "length": 9},
            ],
            "alt_transcript_seq": "ATGGCTCTTGAACGTTGGAAGCCTGCGTTAAAAAGCTGCC",
            "alt_transcript_length": 40,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 8},
                {"exon_number": 2, "length": 14},
                {"exon_number": 3, "length": 18},
            ],
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "alt_transcript",
        },
        equivalent=(Change("CCTTTTCA[GA>]GTGAAC"), Change("CCTTTTCAG[AG>]TGAAC")),
    ),
    Case(
        "insertion_between_an_acceptor_and_a_coding_exon_goes_into_the_exon",
        """
        ref 5' [ATG GCT CT]gtgtaagccccccttttcag[-A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]gtgtaagccccccttttcag[TA GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                                                ^ ->T
        T inserted at the exon 2 start: exon 2 is TAGTGAACGTTGGAAGC (17 nt)
        frameshift: ATG GCT CTT AGT GAA CGT TGG AAG CCT GCG TTA AAA AGC TGC C, nonstop
        """,
        MAIN,
        Change("CCTTTTCAG[>T]AGTGAAC"),
        {
            **record(per_strand("G", "T"), per_strand("GT", "TA"), per_strand(37, 63), per_strand(38, 64)),
            **FRAMESHIFT_NONSTOP,
            "alt_cds_seq": "ATGGCTCTTAGTGAACGTTGGAAGCCTGCGTTAA",
            "alt_cds_length": 34,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 8},
                {"exon_number": 2, "length": 17},
                {"exon_number": 3, "length": 9},
            ],
            "alt_transcript_seq": "ATGGCTCTTAGTGAACGTTGGAAGCCTGCGTTAAAAAGCTGCC",
            "alt_transcript_length": 43,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 8},
                {"exon_number": 2, "length": 17},
                {"exon_number": 3, "length": 18},
            ],
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "alt_transcript",
        },
        equivalent=(Change("CCTTTTCAG[A>TA]GTGAAC"),),
    ),
    Case(
        "insertion_between_a_coding_exon_and_a_donor_goes_into_the_exon",
        """
        ref 5' [ATG GCT CT-]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CTA]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                          ^ ->A
        A inserted at the exon 1 end: exon 1 is ATGGCTCTA (9 nt)
        frameshift: ATG GCT CTA AGT GAA CGT TGG AAG CCT GCG TTA AAA AGC TGC C, nonstop
        """,
        MAIN,
        Change("GGCTCT[>A]GTGTAAG"),
        {
            **record(per_strand("T", "C"), per_strand("TA", "CT"), per_strand(17, 83), per_strand(18, 84)),
            **FRAMESHIFT_NONSTOP,
            "alt_cds_seq": "ATGGCTCTAAGTGAACGTTGGAAGCCTGCGTTAA",
            "alt_cds_length": 34,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 9},
                {"exon_number": 2, "length": 16},
                {"exon_number": 3, "length": 9},
            ],
            "alt_transcript_seq": "ATGGCTCTAAGTGAACGTTGGAAGCCTGCGTTAAAAAGCTGCC",
            "alt_transcript_length": 43,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 9},
                {"exon_number": 2, "length": 16},
                {"exon_number": 3, "length": 18},
            ],
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "alt_transcript",
        },
        equivalent=(Change("GGCTCT[G>AG]TGTAAG"),),
    ),
    Case(
        "insertion_of_ag_at_a_cag_ag_acceptor_is_exon_boundary_ambiguous",
        """
        ref 5' [ATG GCT CT]gtgtaagccccccttttcag[--A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]gtgtaagccccccttttcag[AGA GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                                                ^^ ->AG
        AG inserted: TTTTCAGAG|AG... or TTTTCAG|AGAG...
        every placement keeps an AG before the exon, but exon 2 starts at either AG: exon_boundary_ambiguous
        """,
        MAIN,
        Change("CCTTTTCAG[>AG]AGTGAAC"),
        {
            **record(per_strand("G", "T"), per_strand("GAG", "TCT"), per_strand(37, 63), per_strand(38, 64)),
            **AMBIGUOUS,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("CCTTTTC[>AG]AGAGTGAAC"), Change("CCTTTTCAGAG[>AG]TGAAC")),
    ),
    Case(
        "insertion_of_gt_at_a_ct_gtgt_donor_is_exon_boundary_ambiguous",
        """
        ref 5' [ATG GCT CT--]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CTGT]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                          ^^ ->GT
        GT inserted: exon 1 ends ATGGCTCT|GT... or ATGGCTCTGT|GT...
        every placement keeps a GT after the exon, but at different positions: exon_boundary_ambiguous
        """,
        MAIN,
        Change("GGCTCT[>GT]GTGTAAG"),
        {
            **record(per_strand("T", "C"), per_strand("TGT", "CAC"), per_strand(17, 83), per_strand(18, 84)),
            **AMBIGUOUS,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("GGCTC[>TG]TGTGTAAG"), Change("GGCTCTGTGT[>GT]AAGCC")),
    ),
    Case(
        "deletion_over_the_transcript_start_shortens_exon_1_and_loses_the_start_codon",
        """
        tx                0
        ref 5' cccccccccc[ATG GCT CT]|[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' ccccccccc-[-TG GCT CT]|[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                        ^^^ CA>-
        CA deleted: the last flank base and the A of the start codon
        a transcript start has no splice dinucleotide, so the placement is valid: exon 1 is TGGCTCT (7 nt)
        start loss; the alt transcript TGGCTCTAG... has no ATG: alt_has_ptc and stop_loss False,
        annotated_stop_distance null
        """,
        MAIN,
        Change("CCC[CA>]TGGCTCT"),
        {
            **record(per_strand("CCA", "ATG"), per_strand("C", "A"), per_strand(8, 90), per_strand(11, 93)),
            "alt_cds_seq": "TGGCTCTAGTGAACGTTGGAAGCCTGCGTTAA",
            "alt_cds_length": 32,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 7},
                {"exon_number": 2, "length": 16},
                {"exon_number": 3, "length": 9},
            ],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 6,
            "alt_stop_codon_count": 2,
            "alt_stop_codons": [{"position": 6, "codon": "TAG"}, {"position": 9, "codon": "TGA"}],
            "alt_stop_codon_exons": [1, 2],
            "alt_has_ptc": False,
            "start_loss": True,
            "stop_loss": False,
            "alt_transcript_seq": "TGGCTCTAGTGAACGTTGGAAGCCTGCGTTAAAAAGCTGCC",
            "alt_transcript_length": 41,
            "alt_cds_start_in_transcript": 0,
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
                {"exon_number": 1, "length": 7},
                {"exon_number": 2, "length": 16},
                {"exon_number": 3, "length": 18},
            ],
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "start_loss_scan",
        },
        equivalent=(Change("CC[CCA>C]TGGCTCT"),),
        ruler=Ruler((0,)),
    ),
    Case(
        "deletion_in_an_a_run_across_the_transcript_start_is_exon_boundary_ambiguous",
        """
        tx                0
        ref 5' ccccccccca[ATG GCT CT]|[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' ccccccccca[-TG GCT CT]|[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                          ^ A>-
        one A deleted: the flank A or the A of the start codon
        A transcript start has no splice dinucleotide, so both placements are valid. They put the transcript start
        at different positions: exon_boundary_ambiguous
        """,
        MAIN_AFTER_A,
        Change("CCCCA[A>]TGGCTCT"),
        {
            **record(per_strand("AA", "AT"), per_strand("A", "A"), per_strand(9, 90), per_strand(11, 92)),
            **AMBIGUOUS,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("CCCC[A>]ATGGCTCT"),),
        ruler=Ruler((0,)),
    ),
    # Pins "The shift stops at the chromosome ends. So ... an insertion can go before the first base or after the
    # last base." with "A transcript end has no splice dinucleotide, so placements that put it at different positions
    # make it ambiguous too." The A run is the last 2 bases of the stop codon and the flank A, the last chromosome
    # base. The placement after that base leaves the inserted A outside the transcript, and the placements inside the
    # stop codon put it into the transcript. On the minus strand, the placement before the first chromosome base takes
    # that part. The case deletion_in_a_run_to_the_chromosome_end_makes_the_transcript_end_ambiguous in cases_vcf has
    # its run 10 nt past the transcript end, so its ambiguity does not need the chromosome edge.
    Case(
        "insertion_in_an_a_run_that_ends_at_the_chromosome_edge_1_nt_past_the_transcript_end_is_exon_boundary_ambiguous",
        """
        tx      0   3               15
        ref 5' [gcc ATG GCC AAG CTG TAA-]a 3'
        alt 5' [gcc ATG GCC AAG CTG TAAA]a 3'
                                       ^ ->A
        ->A in the A run of the last 2 stop codon bases and the flank A, the last chromosome base. The A has 4
        placements: before tx 16, before tx 17, before the flank A, and after it, after the last chromosome base.
        The last one leaves the transcript as it is, and the first two put the A into the stop codon. So they put
        the transcript end at different positions: exon_boundary_ambiguous. On the minus strand, the run starts at
        the first chromosome base.
        """,
        STOP_CODON_1_NT_FROM_THE_CHROMOSOME_EDGE,
        Change("CTGTAA[>A]A"),
        {
            **record(per_strand("A", "T"), per_strand("AA", "TT"), per_strand(27, 0), per_strand(28, 1)),
            **AMBIGUOUS,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("CTGT[>A]AA"), Change("CTGTA[>A]A"), Change("CTGTAA[A>AA]")),
        ruler=Ruler((0, 3, 15)),
    ),
    # Pins "The shift stops at the chromosome ends. So a deletion can reach the first or the last base of the
    # chromosome" with "A transcript end has no splice dinucleotide, so placements that put it at different positions
    # make it ambiguous too." Only the placement that deletes the flank A, the last chromosome base, keeps the
    # transcript end; on the minus strand, the flank A is the first chromosome base. The insertion case before pins
    # the other kind of indel.
    Case(
        "deletion_in_an_a_run_that_ends_at_the_chromosome_edge_1_nt_past_the_transcript_end_is_exon_boundary_ambiguous",
        """
        tx      0   3               15
        ref 5' [gcc ATG GCC AAG CTG TAA]a 3'
        alt 5' [gcc ATG GCC AAG CTG TA-]a 3'
                                      ^ A>-
        One A deleted in the A run of the last 2 stop codon bases and the flank A, the last chromosome base. The
        deletion has 3 placements: tx 16, tx 17 and the flank A. Deleting the flank A keeps the transcript, and the
        other two shorten it by 1 nt. So they put the transcript end at different positions:
        exon_boundary_ambiguous. On the minus strand, the run starts at the first chromosome base.
        """,
        STOP_CODON_1_NT_FROM_THE_CHROMOSOME_EDGE,
        Change("CTGTA[A>]A"),
        {
            **record(per_strand("AA", "TT"), per_strand("A", "T"), per_strand(26, 0), per_strand(28, 2)),
            **AMBIGUOUS,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("CTGT[A>]AA"), Change("CTGTA[AA>A]")),
        ruler=Ruler((0, 3, 15)),
    ),
    # Pins "The shift stops at the chromosome ends. So ... an insertion can go before the first base or after the
    # last base. The rules of this section hold for a transcript that starts or ends at a chromosome end too." The A
    # run is the flank A, the first chromosome base, and the A of the start codon. The placement before the first
    # chromosome base leaves the inserted A outside the transcript, and the placement inside the start codon puts it
    # into the transcript. So the two put the transcript start at different positions: "A transcript end has no splice
    # dinucleotide, so placements that put it at different positions make it ambiguous too." On the minus strand, the
    # placement after the last chromosome base takes that part. The case
    # deletion_in_an_a_run_across_the_transcript_start_is_exon_boundary_ambiguous deletes an A in a run across the
    # transcript start, 10 nt from the chromosome start.
    Case(
        "insertion_in_an_a_run_that_starts_at_the_chromosome_edge_1_nt_before_the_transcript_start_is_exon_boundary_ambiguous",
        """
        tx        0     3
        ref 5' a[-ATG]|[GCC AAG CTG TAA gcc] 3'
        alt 5' a[AATG]|[GCC AAG CTG TAA gcc] 3'
                 ^ ->A
        ->A in the A run of the flank A, the first chromosome base, and the A of the start codon. Exon 1 is the
        start codon ATG. The A has 3 placements: before the first chromosome base, before tx 0, and after tx 0,
        inside the start codon. The first leaves the transcript as it is, and the last puts the A into the
        transcript. So they put the transcript start at different positions: exon_boundary_ambiguous. On the minus
        strand, the run ends at the last chromosome base.
        """,
        START_CODON_1_NT_FROM_THE_CHROMOSOME_EDGE,
        Change("[>A]ATGGTAAG"),
        {
            **record(per_strand("A", "T"), per_strand("AA", "TT"), per_strand(0, 47), per_strand(1, 48)),
            **AMBIGUOUS,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("A[>A]TGGTAAG"), Change("[A>AA]ATGGTAAG")),
        ruler=Ruler((0, 3)),
    ),
    Case(
        "mnv_over_a_donor_that_changes_only_the_exon_base_keeps_the_donor",
        """
        ref 5' [ATG GCT CT]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CC]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                         ^^^ TG>CG
        TG>CG: the last exon base T>C, donor +1 G stays; CTA>CCA
        """,
        MAIN,
        Change("GGCTC[TG>CG]TGTAAG"),
        {
            **record(per_strand("TG", "CA"), per_strand("CG", "CG"), per_strand(17, 83), per_strand(19, 85)),
            **KEEPS_THE_STOP_CODON,
            "alt_cds_seq": "ATGGCTCCAGTGAACGTTGGAAGCCTGCGTTAA",
            "alt_transcript_seq": "ATGGCTCCAGTGAACGTTGGAAGCCTGCGTTAAAAAGCTGCC",
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "alt_transcript",
        },
        equivalent=(Change("GGCTC[T>C]GTGTAAG"),),
    ),
    # An equal-length delins maps base for base and keeps GT
    Case(
        "equal_length_delins_over_a_donor_that_keeps_gt_base_for_base",
        """
        ref 5' [ATG GCT CT]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CA]gtctaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                         ^^^^^ TGTG>AGTC
        TGTG>AGTC: CT|GTGT > CA|GTCT; CTA>CAA
        """,
        MAIN,
        Change("GGCTC[TGTG>AGTC]TAAGCC"),
        {
            **record(per_strand("TGTG", "CACA"), per_strand("AGTC", "GACT"), per_strand(17, 81), per_strand(21, 85)),
            **KEEPS_THE_STOP_CODON,
            "alt_cds_seq": "ATGGCTCAAGTGAACGTTGGAAGCCTGCGTTAA",
            "alt_transcript_seq": "ATGGCTCAAGTGAACGTTGGAAGCCTGCGTTAAAAAGCTGCC",
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "alt_transcript",
        },
        equivalent=(Change("GGCT[CTGTG>CAGTC]TAAGCC"),),
    ),
    # An MNV that changes the GT
    Case(
        "mnv_over_a_donor_that_changes_gt_destroys_the_splice_site",
        """
        ref 5' [ATG GCT CT]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CA]ttgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                         ^^^ TG>AT
        TG>AT: CT|GT > CA|TT
        """,
        MAIN,
        Change("GGCTC[TG>AT]TGTAAG"),
        {
            **record(per_strand("TG", "CA"), per_strand("AT", "AT"), per_strand(17, 83), per_strand(19, 85)),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
    ),
    # A GT that the ALT forms elsewhere does not count
    Case(
        "mnv_over_a_donor_whose_alt_forms_gt_one_base_upstream_destroys_the_splice_site",
        """
        ref 5' [ATG GCT CT]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CG]ttgtaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                         ^^^ TG>GT
        TG>GT: CT|GT > CG|TT; the new GT starts 1 nt before the exon end, which does not count
        """,
        MAIN,
        Change("GGCTC[TG>GT]TGTAAG"),
        {
            **record(per_strand("TG", "CA"), per_strand("GT", "AC"), per_strand(17, 83), per_strand(19, 85)),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
    ),
    Case(
        "delins_over_a_donor_that_neither_matching_keeps_destroys_the_splice_site",
        """
        ref 5' [ATG GCT CT]gtg-taagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CG]tccctaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                         ^^^^^^ TGTG>GTCCC
        matched from the left: CG|TCCC TAAG; from the right: CGT|CCC TAAG; neither keeps GT after the exon end.
        The GT that the ALT forms 1 nt before the exon end does not count.
        """,
        MAIN,
        Change("GGCTC[TGTG>GTCCC]TAAGCC"),
        {
            **record(per_strand("TGTG", "CACA"), per_strand("GTCCC", "GGGAC"), per_strand(17, 81), per_strand(21, 85)),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("GGCT[CTGTG>CGTCCC]TAAGCC"),),
    ),
    Case(
        "delins_over_a_donor_that_only_the_right_matching_keeps",
        """
        ref 5' [ATG GCT CT]gtg-taagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CC]cgtataagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                         ^^^^^^ TGTG>CCGTA
        matched from the left: CC|CGTA TAAG, no GT; from the right: CCC|GTA TAAG keeps GT: exon 1 is ATGGCTCCC
        frameshift: ATG GCT CCC AGT GAA CGT TGG AAG CCT GCG TTA AAA AGC TGC C, nonstop
        """,
        MAIN,
        Change("GGCTC[TGTG>CCGTA]TAAGCC"),
        {
            **record(per_strand("TGTG", "CACA"), per_strand("CCGTA", "TACGG"), per_strand(17, 81), per_strand(21, 85)),
            **FRAMESHIFT_NONSTOP,
            "alt_cds_seq": "ATGGCTCCCAGTGAACGTTGGAAGCCTGCGTTAA",
            "alt_cds_length": 34,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 9},
                {"exon_number": 2, "length": 16},
                {"exon_number": 3, "length": 9},
            ],
            "alt_transcript_seq": "ATGGCTCCCAGTGAACGTTGGAAGCCTGCGTTAAAAAGCTGCC",
            "alt_transcript_length": 43,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 9},
                {"exon_number": 2, "length": 16},
                {"exon_number": 3, "length": 18},
            ],
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "alt_transcript",
        },
        equivalent=(Change("GGCT[CTGTG>CCCGTA]TAAGCC"),),
    ),
    Case(
        "delins_over_a_donor_that_only_the_left_matching_keeps",
        """
        ref 5' [ATG GCT CT]gtg-taagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CA]gtcctaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                         ^^^^^^ TGTG>AGTCC
        matched from the left: CA|GTCC TAAG keeps GT; from the right: CAG|TCC TAAG, no GT: exon 1 is ATGGCTCA
        the other 4 ALT bases go into the intron; CTA>CAA
        """,
        MAIN,
        Change("GGCTC[TGTG>AGTCC]TAAGCC"),
        {
            **record(per_strand("TGTG", "CACA"), per_strand("AGTCC", "GGACT"), per_strand(17, 81), per_strand(21, 85)),
            **KEEPS_THE_STOP_CODON,
            "alt_cds_seq": "ATGGCTCAAGTGAACGTTGGAAGCCTGCGTTAA",
            "alt_transcript_seq": "ATGGCTCAAGTGAACGTTGGAAGCCTGCGTTAAAAAGCTGCC",
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "alt_transcript",
        },
        equivalent=(Change("GGCT[CTGTG>CAGTCC]TAAGCC"),),
    ),
    Case(
        "delins_over_a_donor_that_both_matchings_keep_at_different_positions_is_exon_boundary_ambiguous",
        """
        ref 5' [ATG GCT CT]gtg--taagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CA]gtgtctaagccccccttttcag[A GTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                         ^^^^^^^ TGTG>AGTGTC
        matched from the left: CA|GTGTC TAAG; from the right: CAGT|GTC TAAG; both keep GT, 2 nt apart
        """,
        MAIN,
        Change("GGCTC[TGTG>AGTGTC]TAAGCC"),
        {
            **record(
                per_strand("TGTG", "CACA"), per_strand("AGTGTC", "GACACT"), per_strand(17, 81), per_strand(21, 85)
            ),
            **AMBIGUOUS,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("GGCT[CTGTG>CAGTGTC]TAAGCC"),),
    ),
    # Mixing the two matchings would give exon 2 a negative length
    Case(
        "delins_over_an_exon_whose_matchings_keep_one_splice_site_each_destroys_the_splice_site",
        """
        ref 5' [ATG GCT CT]gtgtaagccccccttttcag[A GTG AAC GTT GGA AGC]gtaagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]gtgtaagcccccctttccag[T C-- --- --- --- ---]---agtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
                                           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 23 nt>CCAGTC
        TCAG + exon 2 + GTA > CCAGTC (23 nt > 6 nt)
        matched from the left: TTT CCAG|TC AGTCC keeps the acceptor AG, not the donor GT
        matched from the right: TTT CCA|GTC AGTCC keeps the donor GT, not the acceptor AG
        no placement keeps both; taking each site from another placement gives exon 2 a length of -1
        """,
        MAIN,
        Change("TTT[TCAG" + EXON_2 + "GTA>CCAGTC]AGTCC"),
        {
            **record(
                per_strand("TCAG" + EXON_2 + "GTA", "TACGCTTCCAACGTTCACTCTGA"),
                per_strand("CCAGTC", "GACTGG"),
                per_strand(34, 45),
                per_strand(57, 68),
            ),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("TT[TTCAG" + EXON_2 + "GTA>TCCAGTC]AGTCC"),),
    ),
    # Mixing the two matchings would make exon 1 and exon 2 overlap
    Case(
        "delins_over_an_intron_whose_matchings_keep_one_splice_site_each_destroys_the_splice_site",
        """
        ref 5' [ATG GCT CT]gtgtaagccccccttttcag[A G-----------------------TG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CC]ccccccccccccccccccag[C CCGTCCCCCCCCCCCCCCCCCCCCTG AAC GTT GGA AGC]|[CTG CGT TAA aaagctgcc] 3'
                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 23 nt>46 nt
        T + intron 1 + AG > 19 C, AGCCCGT, 20 C (23 nt > 46 nt); ALT positions count from 0
        matched from the left: the exon 1 end maps to ALT 1 (CC follows, no GT); the exon 2 start maps to ALT 21,
            after the AG at ALT 19 and 20
        matched from the right: the exon 1 end maps to ALT 24, before the GT at ALT 24 and 25; the exon 2 start
            maps to ALT 44 (CC before it, no AG)
        no placement keeps both; taking each site from another placement makes the exons share CCC (ALT 21 to 23)
        """,
        MAIN,
        Change("GGCTC[T" + INTRON_1 + "AG>" + "C" * 19 + "AGCCCGT" + "C" * 20 + "]TGAAC"),
        {
            **record(
                per_strand("T" + INTRON_1 + "AG", "CTCTGAAAAGGGGGGCTTACACA"),
                per_strand("C" * 19 + "AGCCCGT" + "C" * 20, "G" * 20 + "ACGGGCT" + "G" * 19),
                per_strand(17, 62),
                per_strand(40, 85),
            ),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
    ),
    Case(
        "delins_over_a_whole_short_coding_region_whose_edge_rules_take_different_matchings",
        """
        layout         13      19
        ref    5' [gcc ATG TAA g-------cc] 3'
        alt    5' [gct TTT TTT ttttttttcc] 3'
                     ^^^^^^^^^^^^^^^^^^ CATGTAAG>15 nt
        cATGTAAg > 15 T (8 nt > 15 nt, 7 nt longer; the coding region ATGTAA has 6 nt)
        no splice site is near, so both matchings are valid. Neither puts an ATG at the start codon edge, so that
        edge takes the matching from the right: layout 13 > 20. The stop codon edge takes the matching from the
        left: layout 19 > 19. The coding row would end before it starts: splice_site_destroyed
        """,
        SHORT_CDS,
        Change("GC[CATGTAAG>TTTTTTTTTTTTTTT]CC"),
        {
            **record(per_strand("CATGTAAG", "CTTACATG"), per_strand("T" * 15, "A" * 15), 12, 20),
            **DESTROYED,
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": None,
        },
        equivalent=(Change("G[CCATGTAAG>CTTTTTTTTTTTTTTT]CC"),),
        ruler=Ruler((13, 19), "layout"),
    ),
    # Pins "Each placement maps each exon edge base for base to a position in the alt sequence. So a splice
    # dinucleotide that the ALT bases form elsewhere does not count." for a delins with a longer ALT whose REF is the
    # first base of an exon. Matched from the right, the ALT bases GAG put an AG right before the exon base. The exon
    # start still maps to the REF start in both matchings, after the acceptor AG of the reference, so the delins goes
    # into the exon. The closest cases change bases at an exon edge with one placement (the MNV and delins cases over
    # a donor in this file, and in cases_misc
    # mnv_at_cag_ag_whose_alt_forms_another_ag_is_not_ambiguous), or have the exon edge strictly inside REF (the delins
    # cases over a donor in this file).
    Case(
        "delins_at_the_first_exon_base_whose_longer_alt_forms_an_ag_keeps_the_exon_start",
        """
        ref 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtaagtcccccccctttcag[C---TG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtaagtcccccccctttcag[GAGTTG CGT TAA aaagctgcc] 3'
                                                                        ^^^^ C>GAGT
        C>GAGT at the first base of exon 3 (CDS 24), in frame. REF and ALT share no first or last base, so the
        delins has two matchings. Matched from the right, the ALT bases GAG form an AG right before the exon base.
        That AG does not count: in both matchings, the exon start maps to the REF start, after the acceptor AG of
        the reference. Exon 3 is GAGTTGCGTTAAaaagctgcc (21 nt).
        in frame: ATG GCT CTA GTG AAC GTT GGA AGC GAG TTG CGT TAA, the annotated stop codon at CDS 33
        """,
        MAIN,
        Change("TTTCAG[C>GAGT]TGCGT"),
        {
            **record(per_strand("C", "G"), per_strand("GAGT", "ACTC"), per_strand(74, 27), per_strand(75, 28)),
            "alt_cds_seq": "ATGGCTCTAGTGAACGTTGGAAGCGAGTTGCGTTAA",
            "alt_cds_length": 36,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 8},
                {"exon_number": 2, "length": 16},
                {"exon_number": 3, "length": 12},
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
            "alt_transcript_seq": "ATGGCTCTAGTGAACGTTGGAAGCGAGTTGCGTTAAAAAGCTGCC",
            "alt_transcript_length": 45,
            "alt_cds_start_in_transcript": 0,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 8},
                {"exon_number": 2, "length": 16},
                {"exon_number": 3, "length": 21},
            ],
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "alt_transcript",
        },
        equivalent=(Change("TTTCA[GC>GGAGT]TGCGT"),),
    ),
    # Pins "A placement is valid if the reference splice dinucleotide stays next to each of its positions. If all
    # valid placements agree, the exon edges lie there." for a deletion from coding exon 2 across intron 2 into
    # coding exon 3. Both exon edges map to its start. The AG of exon 2 stays before that position and the GT of
    # exon 3 after it, so the alt transcript is known, and the alt bases come from two coding rows. The closest
    # cases delete a whole intron and lose both splice sites
    # (deletion_of_a_whole_intron_between_two_coding_exons_destroys_the_splice_site in cases_coding_region_edges),
    # or delete a whole exon, which joins no two coding rows
    # (deletion_of_a_whole_short_coding_exon_that_keeps_both_splice_sites_empties_the_exon there).
    Case(
        "deletion_across_an_intron_that_keeps_ag_and_gt_joins_two_coding_exons",
        """
        ref 5' [ATG GCT CT]|[A GTG AAC GTT GGA AGC]gtaagtcccccccctttcag[CTG CGT TAA aaagctgcc] 3'
        alt 5' [ATG GCT CT]|[A G-- --- --- --- ---]--------------------[--- -GT TAA aaagctgcc] 3'
                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 38 nt>-
        The deletion runs from exon 2 after its AG, across intron 2, into exon 3 up to its GT: 38 nt, of them 18
        coding nt. It has one placement. Both exon edges map to its start, between the AG of exon 2 and the GT of
        exon 3, so both splice dinucleotides stay. Exon 2 keeps AG (2 nt), and exon 3 keeps GTTAAaaagctgcc (14 nt).
        in frame: ATG GCT CTA GGT TAA, the annotated stop codon at CDS 12
        """,
        MAIN,
        Change("TCAGAG[TGAACGTTGGAAGC" + INTRON + "CTGC>]GTTAA"),
        {
            **record(
                per_strand("GTGAACGTTGGAAGC" + INTRON + "CTGC", "CGCAGCTGAAAGGGGGGGGACTTACGCTTCCAACGTTCA"),
                per_strand("G", "C"),
                per_strand(39, 23),
                per_strand(78, 62),
            ),
            "alt_cds_seq": "ATGGCTCTAGGTTAA",
            "alt_cds_length": 15,
            "alt_cds_exons": [
                {"exon_number": 1, "length": 8},
                {"exon_number": 2, "length": 2},
                {"exon_number": 3, "length": 5},
            ],
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 12,
            "alt_stop_codon_count": 1,
            "alt_stop_codons": [{"position": 12, "codon": "TAA"}],
            "alt_stop_codon_exons": [3],
            "alt_has_ptc": False,
            "start_loss": False,
            "stop_loss": False,
            "alt_transcript_seq": "ATGGCTCTAGGTTAAAAAGCTGCC",
            "alt_transcript_length": 24,
            "alt_cds_start_in_transcript": 0,
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 8},
                {"exon_number": 2, "length": 2},
                {"exon_number": 3, "length": 14},
            ],
            **NOT_SCANNED,
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": 0,
            **NO_RULE,
            "nmd_model_status": "no_ptc",
            "ptc_pos_in_alt_transcript": None,
            "ptc_exon_number": None,
            "stop_classification": "alt_transcript",
        },
        equivalent=(Change("TCAGAG[TGAACGTTGGAAGC" + INTRON + "CTGCG>G]TTAA"),),
    ),
]
