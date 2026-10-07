"""
Conformance cases of how the VCF describes a variant. A symbolic allele, a breakend or the ALT "." or "*" gives no row (SY). The same variant gives
the same row in each of its VCF descriptions: left- or right-aligned, padded, as a delins, in lower case (VD). A REF
equal to the ALT, or a REF that does not match the genome, gives no row.

The drawings show the transcript 5' to 3', also for the minus strand. The runner renders their layout block from the
case. The line ref is the transcript and the line alt is the transcript with the change. The two lines are aligned,
and `-` fills a gap. `[...]` is an exon and `|` is an exon junction. Upper case is the coding region, in codons of the
annotated frame; lower case is UTR. A flank that a description of the change touches is drawn in lower case outside
the brackets. `..N..` leaves out N bases. `^` marks the change as ref>alt. A mark line puts a character under bases,
e.g. `***` under a stop codon, and `<-- label -->` spans a length. A ruler gives tx positions. A case lists the other
descriptions of its variant in equivalent=; the runner checks that each gives the same row. The VCF record is on the
plus strand, so on the minus strand a left-aligned record lies at the 3' end of the change in the transcript. In the
VCF record, vcf_ref replaces the REF and vcf_alt replaces the ALT of the change text, e.g. by a symbolic allele. The
layout block shows the change text.
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
    NoRow,
    Ruler,
    Span,
    Transcript,
    per_strand,
)

# 5' [gacc ATG CAA]|[CTG GCC GCC GCC TTC AAG]|[TGG TAA acg] 3'
#     0    4         10                        28  31  34     tx
#          0         6                         24  27         CDS
THREE_EXONS = Layout(
    Transcript(("gaccATGCAA", "CTGGCCGCCGCCTTCAAG", "TGGTAAacg")),
    {
        **IDS,
        "cds_start": per_strand(14, 13),
        "cds_end": per_strand(84, 83),
        "ref_cds_seq": "ATGCAACTGGCCGCCGCCTTCAAGTGGTAA",
        "ref_cds_length": 30,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [(1, 6), (2, 18), (3, 6)],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 27,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [(27, "TAA")],
        "ref_stop_codon_exons": [3],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 87,
        "transcript_seq": "GACCATGCAACTGGCCGCCGCCTTCAAGTGGTAAACG",
        "transcript_length": 37,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 34,
        "transcript_exons": [(1, 10), (2, 18), (3, 9)],
        "utr3_length": 3,
        "utr5_length": 4,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)
# The alt columns of a variant in THREE_EXONS that changes no CDS length and no stop codon
THREE_EXONS_UNCHANGED_STOP = {
    "alt_cds_length": 30,
    "alt_cds_exons": [(1, 6), (2, 18), (3, 6)],
    "alt_last_codon": "TAA",
    "alt_valid_stop": True,
    "alt_first_stop_codon": "TAA",
    "alt_first_stop_pos": 27,
    "alt_stop_codon_count": 1,
    "alt_stop_codons": [(27, "TAA")],
    "alt_stop_codon_exons": [3],
    "alt_has_ptc": False,
    "start_loss": False,
    "stop_loss": False,
    "alt_transcript_length": 37,
    "alt_cds_start_in_transcript": 4,
    **NOT_SCANNED,
    "unknown_reason": None,
    **NO_PTC_FEATURES,
    "annotated_stop_distance": 0,
    **NO_RULE,
}
# The layout block of the SNV T>A at the second T of TTC in THREE_EXONS, tx 23, CDS 19. A case adds its prose.
THREE_EXONS_SNV_DRAWING = """
ref 5' [gacc ATG CAA]|[CTG GCC GCC GCC TTC AAG]|[TGG TAA acg] 3'
alt 5' [gacc ATG CAA]|[CTG GCC GCC GCC TAC AAG]|[TGG TAA acg] 3'
                                        ^ T>A
"""
THREE_EXONS_SNV_RECORD = {
    "variant_id": "var1",
    "ref": per_strand("T", "A"),
    "variant_start": per_strand(53, 43),
    "variant_end": per_strand(54, 44),
}

# The transcripts of the stop codon cases: the CDS ATG GCC GCC GCC | CTG ACC GCC GCC <last codon>, the stop codon and
# a 3'UTR. Only the last codon, the stop codon and the 3'UTR differ.
#
# 5' [ccgccgccaccgc ATG GCC GCC GCC]|[CTG ACC GCC GCC <last> <stop> <3'UTR>] 3'
#     0             13                25              37     40     43          tx
#                   0                 12              24     27     30          CDS
STOP_UTR5 = "ccgccgccaccgc"
STOP_CDS_HEAD = "ATGGCCGCCGCC"
STOP_CDS_BODY = "CTGACCGCCGCC"
STOP_SHARED = {
    **IDS,
    "ref_cds_length": 30,
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "ref_cds_exons": [(1, 12), (2, 18)],
    "cds_in_transcript": True,
    "start_codon_exon": 1,
    "ref_valid_stop": True,
    "ref_first_stop_pos": 27,
    "ref_stop_codon_count": 1,
    "ref_stop_codon_exons": [2],
    "ref_has_ptc": False,
    "transcript_start": 10,
    "cds_start_in_transcript": 13,
    "cds_end_in_transcript": 43,
    "utr5_length": 13,
    "total_exon_count": 2,
    "likely_misannotated": False,
}
# The alt columns of a variant in a stop codon layout that keeps the alt CDS as the ref CDS (TAA at CDS 27)
STOP_ALT_CDS_UNCHANGED = {
    "alt_cds_length": 30,
    "alt_cds_exons": [(1, 12), (2, 18)],
    "alt_last_codon": "TAA",
    "alt_valid_stop": True,
    "alt_first_stop_codon": "TAA",
    "alt_first_stop_pos": 27,
    "alt_stop_codon_count": 1,
    "alt_stop_codons": [(27, "TAA")],
    "alt_stop_codon_exons": [2],
    "alt_has_ptc": False,
    "start_loss": False,
    "stop_loss": False,
    "alt_cds_start_in_transcript": 13,
    **NOT_SCANNED,
    "unknown_reason": None,
    **NO_PTC_FEATURES,
    "annotated_stop_distance": 0,
    **NO_RULE,
}


def stop_layout(last_codon, stop_codon, utr3, ref_values):
    """A stop codon layout with the given last codon, stop codon and 3'UTR, and its expected REF_COLUMNS values."""
    exons = (STOP_UTR5 + STOP_CDS_HEAD, STOP_CDS_BODY + last_codon + stop_codon + utr3)
    return Layout(Transcript(exons), {**STOP_SHARED, **ref_values})


# Last codon GCC, stop codon TAA, 3'UTR cccc tga cc tag ccc
STOP_GCC_TAA = stop_layout(
    "GCC",
    "TAA",
    "cccctgacctagccc",
    {
        "cds_start": per_strand(23, 25),
        "cds_end": per_strand(73, 75),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTAA",
        "ref_last_codon": "TAA",
        "ref_first_stop_codon": "TAA",
        "ref_stop_codons": [(27, "TAA")],
        "transcript_end": 88,
        "transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCGCCTAA" + "CCCCTGACCTAGCCC",
        "transcript_length": 58,
        "transcript_exons": [(1, 25), (2, 33)],
        "utr3_length": 15,
    },
)
# Last codon TCC, else as STOP_GCC_TAA
STOP_TCC_TAA = stop_layout(
    "TCC",
    "TAA",
    "cccctgacctagccc",
    {
        "cds_start": per_strand(23, 25),
        "cds_end": per_strand(73, 75),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTCCTAA",
        "ref_last_codon": "TAA",
        "ref_first_stop_codon": "TAA",
        "ref_stop_codons": [(27, "TAA")],
        "transcript_end": 88,
        "transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCTCCTAA" + "CCCCTGACCTAGCCC",
        "transcript_length": 58,
        "transcript_exons": [(1, 25), (2, 33)],
        "utr3_length": 15,
    },
)
# Last codon TGG, stop codon TAA, 3'UTR actg ccc: the run TAA|a of 3 A
STOP_TGG_TAA_A = stop_layout(
    "TGG",
    "TAA",
    "actgccc",
    {
        "cds_start": per_strand(23, 17),
        "cds_end": per_strand(73, 67),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTGGTAA",
        "ref_last_codon": "TAA",
        "ref_first_stop_codon": "TAA",
        "ref_stop_codons": [(27, "TAA")],
        "transcript_end": 80,
        "transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCTGGTAA" + "ACTGCCC",
        "transcript_length": 50,
        "transcript_exons": [(1, 25), (2, 25)],
        "utr3_length": 7,
    },
)
# Last codon GCC, stop codon TAA, 3'UTR cc tag ccc
STOP_GCC_TAA_CCTAG = stop_layout(
    "GCC",
    "TAA",
    "cctagccc",
    {
        "cds_start": per_strand(23, 18),
        "cds_end": per_strand(73, 68),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTAA",
        "ref_last_codon": "TAA",
        "ref_first_stop_codon": "TAA",
        "ref_stop_codons": [(27, "TAA")],
        "transcript_end": 81,
        "transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCGCCTAA" + "CCTAGCCC",
        "transcript_length": 51,
        "transcript_exons": [(1, 25), (2, 26)],
        "utr3_length": 8,
    },
)
# Last codon GTA, stop codon TAG, 3'UTR tag cat ccc: a stop codon repeat TAG TAG
STOP_GTA_TAG_TAG = stop_layout(
    "GTA",
    "TAG",
    "tagcatccc",
    {
        "cds_start": per_strand(23, 19),
        "cds_end": per_strand(73, 69),
        "ref_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGTATAG",
        "ref_last_codon": "TAG",
        "ref_first_stop_codon": "TAG",
        "ref_stop_codons": [(27, "TAG")],
        "transcript_end": 82,
        "transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCGTATAG" + "TAGCATCCC",
        "transcript_length": 52,
        "transcript_exons": [(1, 25), (2, 27)],
        "utr3_length": 9,
    },
)

# A transcript whose 3' end lies in a run of A that goes on to the end of the chromosome: TAA, the 3'UTR aaa and a
# 3' flank of 10 A. On the minus strand, the run starts at the start of the chromosome.
#
# 5' [gccacc ATG GCC]|[AAG CTG TAA aaa] AAAAAAAAAA 3'
#     0      6        12      18  21                 tx
#            0        6       12                     CDS
RUN_TO_CHROMOSOME_END = Layout(
    Transcript(("gccaccATGGCC", "AAGCTGTAAaaa"), flanks=("CCCCCCCCCC", "AAAAAAAAAA")),
    {
        **IDS,
        "cds_start": per_strand(16, 13),
        "cds_end": per_strand(51, 48),
        "ref_cds_seq": "ATGGCCAAGCTGTAA",
        "ref_cds_length": 15,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_exons": [(1, 6), (2, 9)],
        "cds_in_transcript": True,
        "start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 12,
        "ref_stop_codon_count": 1,
        "ref_stop_codons": [(12, "TAA")],
        "ref_stop_codon_exons": [2],
        "ref_has_ptc": False,
        "transcript_start": 10,
        "transcript_end": 54,
        "transcript_seq": "GCCACCATGGCCAAGCTGTAAAAA",
        "transcript_length": 24,
        "cds_start_in_transcript": 6,
        "cds_end_in_transcript": 21,
        "transcript_exons": [(1, 12), (2, 12)],
        "utr3_length": 3,
        "utr5_length": 6,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)


# The symbolic alleles and breakends at the SNV site of THREE_EXONS, as (name, ALT as written, reason for no row).
# {ref} is the REF of the record: T on the plus strand, A on the minus strand. Each case also pins SY-05: a VCF with
# only symbolic records gives no row, with all columns and their dtypes.
SYMBOLIC_ALTS = [
    # SY-01: the symbolic alleles of structural variants
    ("symbolic_deletion", "<DEL>", "symbolic allele"),
    ("symbolic_duplication", "<DUP>", "symbolic allele"),
    ("symbolic_insertion", "<INS>", "symbolic allele"),
    ("symbolic_inversion", "<INV>", "symbolic allele"),
    ("symbolic_copy_number_variant", "<CNV>", "symbolic allele"),
    ("symbolic_tandem_duplication", "<DUP:TANDEM>", "symbolic allele"),
    # SY-02: the unspecified allele and a symbolic allele with subtypes
    ("symbolic_unspecified_allele", "<*>", "symbolic allele"),
    ("symbolic_mobile_element_insertion", "<INS:ME:ALU>", "symbolic allele"),
    # SY-03: the 4 breakend forms. t]p]: the reverse complement of the piece left of p is joined after t. ]p]t: the
    # piece left of p is joined before t. t[p[: the piece right of p is joined after t. [p[t: the reverse complement
    # of the piece right of p is joined before t.
    ("breakend_ref_then_closing_brackets", "{ref}]chr2:100]", "breakend"),
    ("breakend_closing_brackets_then_ref", "]chr2:100]{ref}", "breakend"),
    ("breakend_ref_then_opening_brackets", "{ref}[chr2:100[", "breakend"),
    ("breakend_opening_brackets_then_ref", "[chr2:100[{ref}", "breakend"),
    # SY-04: single breakends, also with inserted bases
    ("single_breakend_after_ref", "{ref}.", "breakend"),
    ("single_breakend_before_ref", ".{ref}", "breakend"),
    ("single_breakend_with_inserted_bases", "{ref}TA.", "breakend"),
    # SY-08: a single symbolic allele or breakend with a pipe in its text is no multi-allelic record
    ("symbolic_allele_with_a_pipe_in_its_id", "<INS:ME|ALU>", "symbolic allele"),
    ("breakend_to_a_contig_with_a_pipe_after_ref", "{ref}]gi|123|:100]", "breakend"),
    ("breakend_to_a_contig_with_a_pipe_before_ref", "[gi|123|:100[{ref}", "breakend"),
]
SYMBOLIC_CASES = [
    Case(
        f"{name}_at_a_cds_position_gives_no_row",
        f"{THREE_EXONS_SNV_DRAWING}The VCF record of T>A at tx 23, CDS 19 has the ALT {alt.format(ref='T')}.\n",
        THREE_EXONS,
        Change("GCCT[T>A]CAAG", vcf_alt=alt),
        NoRow(reason),
    )
    for name, alt, reason in SYMBOLIC_ALTS
]

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

CASES = [
    *SYMBOLIC_CASES,
    # VD-05, VD-06, VD-10, SY-06: an SNV, as a padded MNV, and with REF or ALT in lower case
    Case(
        "missense_snv_in_an_internal_exon_gives_one_row_for_each_description",
        THREE_EXONS_SNV_DRAWING
        + """T>A at tx 23, CDS 19: TTC>TAC
Descriptions: T>A; TT>TA and TC>AC (one padding base); TTC>TAC (padding on both sides); T>a, t>A and t>a.
""",
        THREE_EXONS,
        Change("GCCT[T>A]CAAG"),
        {
            **THREE_EXONS_SNV_RECORD,
            **THREE_EXONS_UNCHANGED_STOP,
            "alt": per_strand("A", "T"),
            "alt_cds_seq": "ATGCAACTGGCCGCCGCCTACAAGTGGTAA",
            "alt_transcript_seq": "GACCATGCAACTGGCCGCCGCCTACAAGTGGTAAACG",
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        equivalent=(
            Change("GCC[TT>TA]CAAG"),
            Change("GCCT[TC>AC]AAG"),
            Change("GCC[TTC>TAC]AAG"),
            Change("GCCT[T>A]CAAG", lower_case=("alt",)),
            Change("GCCT[T>A]CAAG", lower_case=("ref",)),
            Change("GCCT[T>A]CAAG", lower_case=("ref", "alt")),
        ),
    ),
    # SY-06: an ALT with N is a sequence, not a symbolic allele
    Case(
        "snv_to_n_puts_n_into_the_alt_cds",
        THREE_EXONS_SNV_DRAWING + "The VCF record of T>A at tx 23, CDS 19 has the ALT N: TTC>TNC.\n",
        THREE_EXONS,
        Change("GCCT[T>A]CAAG", vcf_alt="N"),
        {
            **THREE_EXONS_SNV_RECORD,
            **THREE_EXONS_UNCHANGED_STOP,
            "alt": "N",
            "alt_cds_seq": "ATGCAACTGGCCGCCGCCTNCAAGTGGTAA",
            "alt_transcript_seq": "GACCATGCAACTGGCCGCCGCCTNCAAGTGGTAAACG",
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
    ),
    # VD-09: REF equal to ALT, also if only the case differs
    Case(
        "ref_equal_to_alt_at_a_cds_position_gives_no_row",
        """
        ref 5' [gacc ATG CAA]|[CTG GCC GCC GCC TTC AAG]|[TGG TAA acg] 3'
                                                ^ T>T
        T>T at tx 23, CDS 19: the ALT equals the REF.
        Descriptions: T>T, T>t, TTC>TTC.
        """,
        THREE_EXONS,
        Change("GCCT[T>T]CAAG"),
        NoRow("touches no coding region"),
        equivalent=(Change("GCCT[T>T]CAAG", lower_case=("alt",)), Change("GCC[TTC>TTC]AAG")),
    ),
    # VD-11: the REF check covers the whole REF, not only the padding base
    Case(
        "deletion_whose_ref_mismatches_after_the_padding_base_gives_no_row",
        """
        ref 5' [gacc ATG CAA]|[CTG GCC GCC GCC TTC AAG]|[TGG TAA acg] 3'
        alt 5' [gacc ATG CAA]|[CTG GCC GCC GCC --C AAG]|[TGG TAA acg] 3'
                                               ^^ TT>-
        TT deleted at tx 22-23. The VCF record gives GA as the deleted bases: REF CGA on the plus strand, GTC on
        the minus strand. Its padding base matches the genome, its other 2 bases do not.
        """,
        THREE_EXONS,
        Change("GCC[TT>]CAAG", vcf_ref="GA"),
        NoRow("REF mismatch"),
    ),
    # VD-16, VD-05, VD-06: an in-frame deletion in a codon repeat, left- or right-aligned, shifted, or as a delins
    Case(
        "inframe_deletion_in_a_codon_repeat_gives_one_row_for_each_alignment",
        """
        tx                         13      19
        ref 5' [gacc ATG CAA]|[CTG GCC GCC GCC TTC AAG]|[TGG TAA acg] 3'
        alt 5' [gacc ATG CAA]|[CTG --- GCC GCC TTC AAG]|[TGG TAA acg] 3'
                                   ^^^ GCC>-
        One GCC deleted from the run GCC GCC GCC at tx 13-21.
        Descriptions: GCC at tx 13 (left-aligned), CCG at tx 14, GCC at tx 19 (right-aligned), and the delins
        GGCCG>GG at tx 12.
        """,
        THREE_EXONS,
        Change("CTG[GCC>]GCCGCCTTC"),
        {
            **THREE_EXONS_UNCHANGED_STOP,
            "variant_id": "var1",
            "ref": per_strand("GGCC", "CGGC"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(42, 50),
            "variant_end": per_strand(46, 54),
            "alt_cds_seq": "ATGCAACTGGCCGCCTTCAAGTGGTAA",
            "alt_cds_length": 27,
            "alt_cds_exons": [(1, 6), (2, 15), (3, 6)],
            "alt_first_stop_pos": 24,
            "alt_stop_codons": [(24, "TAA")],
            "alt_transcript_seq": "GACCATGCAACTGGCCGCCTTCAAGTGGTAAACG",
            "alt_transcript_length": 34,
            "alt_transcript_exons": [(1, 10), (2, 15), (3, 9)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(
            Change("CTGGCCGCC[GCC>]TTC"),
            Change("CTGG[CCG>]CCGCCTTC"),
            Change("CT[GGCCG>GG]CCGCCTTC"),
        ),
        ruler=Ruler((13, 19)),
    ),
    # VD-01, VD-12: a 1 nt deletion in a run that starts in the stop codon, described at each of its 3 positions.
    # On the minus strand, the 3'-most description is anchored on a 3'UTR base. The record AC>C lies in the 3'UTR
    # on both strands.
    Case(
        "deletion_of_one_a_in_the_run_from_the_stop_codon_into_the_3utr_gives_one_row_for_each_description",
        """
        tx                                                       41 43
        ref 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC TGG TAA actgccc] 3'
        alt 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC TGG T-A actgccc] 3'
                                                                 ^ A>-
        One A deleted from the run AAa at tx 41-43.
        Descriptions: the A at tx 41, 42 or 43, and AC>C at tx 43. The farthest placement into the 3'UTR deletes
        tx 43, so the CDS stays, and the 3'UTR loses the A. The alt line shows the change text: it deletes the A
        at tx 41.
        """,
        STOP_TGG_TAA_A,
        Change("TGGT[A>]AACTG"),
        {
            **STOP_ALT_CDS_UNCHANGED,
            "variant_id": "var1",
            "ref": per_strand("TA", "TT"),
            "alt": per_strand("T", "T"),
            "variant_start": per_strand(70, 17),
            "variant_end": per_strand(72, 19),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTGGTAA",
            "alt_transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCTGGTAA" + "CTGCCC",
            "alt_transcript_length": 49,
            "alt_transcript_exons": [(1, 25), (2, 24)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("TGGTA[A>]ACTG"), Change("TGGTAA[A>]CTG"), Change("TGGTAA[AC>C]TG")),
        ruler=Ruler((41, 43)),
    ),
    # VD-02, VD-14: TCC inserted right before the stop codon TAA, described at each of its 4 positions. The 3'-most
    # placement is CCT inserted after the T of the stop codon.
    Case(
        "tcc_inserted_right_before_the_stop_codon_gives_one_row_for_each_description",
        """
        tx                                                  37     40
        ref 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC G---CC TAA cccctgacctagccc] 3'
        alt 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC GCCTCC TAA cccctgacctagccc] 3'
                                                             ^^^ ->CCT
        TCC inserted right before the stop codon: alt GCC TCC TAA.
        Descriptions: CCT after tx 37, CTC after tx 38, TCC after tx 39, CCT after tx 40 (inside the stop codon).
        The alt transcript reads TAA at the shifted position of the annotated stop codon: neither a PTC nor a
        stop loss.
        """,
        STOP_GCC_TAA,
        Change("GCCGCCG[>CCT]CCTAA"),
        {
            **STOP_ALT_CDS_UNCHANGED,
            "variant_id": "var1",
            "ref": per_strand("G", "G"),
            "alt": per_strand("GCCT", "GAGG"),
            "variant_start": per_strand(67, 29),
            "variant_end": per_strand(68, 30),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTCCTAA",
            "alt_cds_length": 33,
            "alt_cds_exons": [(1, 12), (2, 21)],
            "alt_first_stop_pos": 30,
            "alt_stop_codons": [(30, "TAA")],
            "alt_transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCGCCTCCTAA" + "CCCCTGACCTAGCCC",
            "alt_transcript_length": 61,
            "alt_transcript_exons": [(1, 25), (2, 36)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(
            Change("GCCGCCGC[>CTC]CTAA"),
            Change("GCCGCCGCC[>TCC]TAA"),
            Change("GCCGCCGCCT[>CCT]AA"),
        ),
        ruler=Ruler((37, 40)),
    ),
    # VD-14: TGGCCC inserted right before the stop codon TAA, described at each of its 4 positions
    Case(
        "tggccc_inserted_right_before_the_stop_codon_gives_one_row_for_each_description",
        """
        tx                                                  37        40
        ref 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC ------TAA cccctgacctagccc] 3'
        alt 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TGGCCCTAA cccctgacctagccc] 3'
                                                                ^^^^^^ ->TGGCCC
        TGGCCC inserted right before the stop codon: alt GCC TGG CCC TAA.
        Descriptions: CCTGGC after tx 37, CTGGCC after tx 38, TGGCCC after tx 39, GGCCCT after tx 40 (inside the
        stop codon). Neither a PTC nor a stop loss.
        """,
        STOP_GCC_TAA,
        Change("GCCGCCGCC[>TGGCCC]TAA"),
        {
            **STOP_ALT_CDS_UNCHANGED,
            "variant_id": "var1",
            "ref": per_strand("C", "A"),
            "alt": per_strand("CTGGCCC", "AGGGCCA"),
            "variant_start": per_strand(69, 27),
            "variant_end": per_strand(70, 28),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTGGCCCTAA",
            "alt_cds_length": 36,
            "alt_cds_exons": [(1, 12), (2, 24)],
            "alt_first_stop_pos": 33,
            "alt_stop_codons": [(33, "TAA")],
            "alt_transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCGCCTGGCCCTAA" + "CCCCTGACCTAGCCC",
            "alt_transcript_length": 64,
            "alt_transcript_exons": [(1, 25), (2, 39)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(
            Change("GCCGCCG[>CCTGGC]CCTAA"),
            Change("GCCGCCGC[>CTGGCC]CTAA"),
            Change("GCCGCCGCCT[>GGCCCT]AA"),
        ),
        ruler=Ruler((37, 40)),
    ),
    # VD-14, VD-16: the last sense codon TCC deleted, described at each of its 4 positions. The 3'-most placement
    # deletes CCT, with the T of the stop codon.
    Case(
        "last_sense_codon_tcc_deleted_before_the_stop_codon_gives_one_row_for_each_description",
        """
        tx                                               35  38
        ref 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC TCC TAA cccctgacctagccc] 3'
        alt 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC --- TAA cccctgacctagccc] 3'
                                                            ^^^ TCC>-
        3 nt deleted from tx 35-40: alt GCC GCC TAA.
        Descriptions: CCT at tx 35, TCC at tx 37 (as the old CTCC>C), CCT at tx 38 (as the old TCCT>T, with the T
        of the stop codon). Neither a PTC nor a stop loss.
        """,
        STOP_TCC_TAA,
        Change("GCCGCC[TCC>]TAA"),
        {
            **STOP_ALT_CDS_UNCHANGED,
            "variant_id": "var1",
            "ref": per_strand("CTCC", "AGGA"),
            "alt": per_strand("C", "A"),
            "variant_start": per_strand(66, 27),
            "variant_end": per_strand(70, 31),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCTAA",
            "alt_cds_length": 27,
            "alt_cds_exons": [(1, 12), (2, 15)],
            "alt_first_stop_pos": 24,
            "alt_stop_codons": [(24, "TAA")],
            "alt_transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCTAA" + "CCCCTGACCTAGCCC",
            "alt_transcript_length": 55,
            "alt_transcript_exons": [(1, 25), (2, 30)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCCG[CCT>]CCTAA"), Change("GCCGCCT[CCT>]AA")),
        ruler=Ruler((35, 38)),
    ),
    # VD-03: a 5 nt deletion across the stop codon, described from the CDS or into the 3'UTR
    Case(
        "deletion_across_the_stop_codon_gives_one_row_for_each_description",
        """
        tx                                                   38   42   46
        ref 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA cctagccc] 3'
        alt 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC GC- --- -ctagccc] 3'
                                                              ^^^^^^^ CTAAC>-
                                                                      *** TAG at the position of the stop codon
        5 nt deleted from tx 38-46: alt GCC GCC TAG ccc.
        Descriptions: CCTAA at tx 38, CTAAC at tx 39 (as the old CCTAAC>C), TAACC at tx 40 (as the old
        CTAACC>C), ACCTA at tx 42. The TAG of the 3'UTR lands at the position of the stop codon: neither a PTC
        nor a stop loss.
        """,
        STOP_GCC_TAA_CCTAG,
        Change("GCCGCCGC[CTAAC>]CTAG"),
        {
            **STOP_ALT_CDS_UNCHANGED,
            "variant_id": "var1",
            "ref": per_strand("CCTAAC", "GGTTAG"),
            "alt": per_strand("C", "G"),
            "variant_start": per_strand(68, 16),
            "variant_end": per_strand(74, 22),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTA",
            "alt_cds_length": 29,
            "alt_cds_exons": [(1, 12), (2, 17)],
            "alt_last_codon": "CTA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_stop_codon_count": 0,
            "alt_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCGC" + "CTAGCCC",
            "alt_transcript_length": 46,
            "alt_transcript_exons": [(1, 25), (2, 21)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(
            Change("GCCGCCG[CCTAA>]CCTAG"),
            Change("GCCGCCGCC[TAACC>]TAG"),
            Change("GCCGCCGCCTA[ACCTA>]G"),
        ),
        marks=(Mark("alt", 40, 43, "*", "TAG at the position of the stop codon"),),
        ruler=Ruler((38, 42, 46)),
    ),
    # VD-04: a deletion from the stop codon into the 3'UTR, with the VCF anchor in the CDS or in the 3'UTR
    Case(
        "deletion_from_the_stop_codon_into_the_3utr_gives_one_row_for_each_anchor",
        """
        tx                                                      40    45
        ref 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA cccctgacctagccc] 3'
        alt 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC T-- --cctgacctagccc] 3'
                                                                 ^^^^^ AACC>-
                                                                        *** TGA, the first stop codon in frame
                                                                <------> annotated_stop_distance = -3
        AAcc deleted at tx 41-44: alt GCC TCC TGA.
        Descriptions: TAACC>T anchored on tx 40 (CDS) and AACCC>C anchored on tx 45 (3'UTR), each on both
        strands. Read in frame, the TGA of the 3'UTR is the first stop codon, 3 nt downstream of the annotated
        one: a stop loss.
        """,
        STOP_GCC_TAA,
        Change("GCCT[AACC>]CCTGA"),
        {
            "variant_id": "var1",
            "ref": per_strand("TAACC", "GGGTT"),
            "alt": per_strand("T", "G"),
            "variant_start": per_strand(70, 22),
            "variant_end": per_strand(75, 27),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCT",
            "alt_cds_length": 28,
            "alt_cds_exons": [(1, 12), (2, 16)],
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
            "alt_transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCGCCT" + "CCTGACCTAGCCC",
            "alt_transcript_length": 54,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TGA",
            "alt_scan_first_stop_pos": 43,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [(43, "TGA")],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -3,
            **NO_RULE,
            "alt_transcript_exons": [(1, 25), (2, 29)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("GCC[TAACC>T]CCTGA"), Change("GCCT[AACCC>C]CTGA")),
        marks=(
            Mark("alt", 43, 46, "*", "TGA, the first stop codon in frame"),
            Span("alt", 40, 43, "annotated_stop_distance = -3"),
        ),
        ruler=Ruler((40, 45)),
    ),
    # VD-07: a delins over the last base of the stop codon and the first base of the 3'UTR. Of its two placements,
    # the one with the deletion at the 3'UTR end counts: A>G turns TAA into TAG, and the 3'UTR loses a c.
    Case(
        "delins_over_the_stop_codon_end_puts_its_deletion_into_the_3utr",
        """
        tx                                                        42
        ref 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAA cccctgacctagccc] 3'
        alt 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC GCC TAG -ccctgacctagccc] 3'
                                                                  ^^^ AC>G
        Ac>G at tx 42-43: alt GCC TAG ccctga...
        """,
        STOP_GCC_TAA,
        Change("GCCTA[AC>G]CCCTGA"),
        {
            **STOP_ALT_CDS_UNCHANGED,
            "variant_id": "var1",
            "ref": per_strand("AC", "GT"),
            "alt": per_strand("G", "C"),
            "variant_start": per_strand(72, 24),
            "variant_end": per_strand(74, 26),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGCCTAG",
            "alt_last_codon": "TAG",
            "alt_first_stop_codon": "TAG",
            "alt_stop_codons": [(27, "TAG")],
            "alt_transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCGCCTAG" + "CCCTGACCTAGCCC",
            "alt_transcript_length": 57,
            "alt_transcript_exons": [(1, 25), (2, 32)],
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((42,)),
    ),
    # VD-08: the delins ATAG>T is not split into an SNV plus an indel. GTA TAG TAG becomes GTT TAG: a stop loss,
    # although the protein stays the same.
    Case(
        "delins_atag_to_t_in_a_stop_codon_repeat_gives_the_same_row_for_each_description",
        """
        ref 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC GTA TAG tagcatccc] 3'
        alt 5' [ccg..10.. ATG GCC GCC GCC]|[CTG ACC GCC GCC GTT --- tagcatccc] 3'
                                                              ^^^^^ ATAG>T
                                                                    *** TAG, the first stop codon in frame
                                                            <-> annotated_stop_distance = -3
        ATAG>T at tx 39-42: alt GTT tag cat ccc.
        Descriptions: ATAG>T, TATAG>TT and ATAGT>TT. The annotated stop codon maps to tx 37 of the alt
        transcript, and the first stop codon in frame is the TAG at tx 40: annotated_stop_distance -3.
        """,
        STOP_GTA_TAG_TAG,
        Change("GT[ATAG>T]TAGCAT"),
        {
            "variant_id": "var1",
            "ref": per_strand("ATAG", "CTAT"),
            "alt": per_strand("T", "A"),
            "variant_start": per_strand(69, 19),
            "variant_end": per_strand(73, 23),
            "alt_cds_seq": "ATGGCCGCCGCCCTGACCGCCGCCGTT",
            "alt_cds_length": 27,
            "alt_cds_exons": [(1, 12), (2, 15)],
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
            "alt_transcript_seq": "CCGCCGCCACCGC" + "ATGGCCGCCGCCCTGACCGCCGCCGTT" + "TAGCATCCC",
            "alt_transcript_length": 49,
            "alt_cds_start_in_transcript": 13,
            "alt_scan_start_codon_pos": 13,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_codon": "TAG",
            "alt_scan_first_stop_pos": 40,
            "alt_scan_stop_codon_count": 1,
            "alt_scan_stop_codons": [(40, "TAG")],
            "alt_scan_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "annotated_stop_distance": -3,
            **NO_RULE,
            "alt_transcript_exons": [(1, 25), (2, 24)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("G[TATAG>TT]TAGCAT"), Change("GT[ATAGT>TT]AGCAT")),
        marks=(
            Mark("alt", 40, 43, "*", "TAG, the first stop codon in frame"),
            Span("alt", 37, 40, "annotated_stop_distance = -3"),
        ),
    ),
    # VD-13: the placements of a deletion in a run stop at the end of the chromosome. Some of them lie inside the
    # transcript and some after its 3' end, so they put the transcript end at different positions.
    Case(
        "deletion_in_a_run_to_the_chromosome_end_makes_the_transcript_end_ambiguous",
        """
        tx                                19 21 24
        ref 5' [gccacc ATG GCC]|[AAG CTG TAA aaa]aaaaaaaaaa 3'
        alt 5' [gccacc ATG GCC]|[AAG CTG T-A aaa]aaaaaaaaaa 3'
                                          ^ A>-
        One A deleted, anywhere from tx 19 to the end of the chromosome: the run goes on into the 3' flank of 10 A.
        Descriptions: the A at tx 19 (stop codon), at tx 21 (3'UTR) and the third A after the transcript end.
        """,
        RUN_TO_CHROMOSOME_END,
        Change("CTGT[A>]AAAA"),
        {
            **UNKNOWN_ALT,
            "variant_id": "var1",
            "ref": per_strand("TA", "TT"),
            "alt": per_strand("T", "T"),
            "variant_start": per_strand(48, 13),
            "variant_end": per_strand(50, 15),
            "unknown_reason": "exon_boundary_ambiguous",
        },
        equivalent=(Change("CTGTAA[A>]AA"), Change("CTGTAAAAAAA[A>]A")),
        ruler=Ruler((19, 21, 24)),
    ),
    # SY-09. "." means that the record has no alternate allele, so it changes no base. The docs skip a record with
    # the ALT "." or "*" with a warning ("Output columns").
    Case(
        "alt_dot_without_an_alternate_allele_gives_no_row",
        """
        ref 5' [gacc ATG GAT G]|[TA AGC TAA gc] 3'
        alt 5' [gacc ATG GAT G]|[TA AGA TAA gc] 3'
                                      ^ C>A
        The VCF record of C>A has the ALT "." (no alternate allele).
        """,
        TWO_EXONS,
        Change("AAG[C>A]TAA", vcf_alt="."),
        NoRow('ALT "." or "*"'),
    ),
    # SY-10. "*" stands for the bases that an overlapping deletion removes, and that deletion has its own record. So
    # the "*" record changes no base itself. The docs skip a record with the ALT "." or "*" with a warning.
    Case(
        "alt_star_of_an_overlapping_deletion_gives_no_row",
        """
        ref 5' [gacc ATG GAT G]|[TA AGC TAA gc] 3'
        alt 5' [gacc ATG GAT G]|[TA AGA TAA gc] 3'
                                      ^ C>A
        The VCF record of C>A has the ALT "*" (allele removed by an overlapping deletion).
        """,
        TWO_EXONS,
        Change("AAG[C>A]TAA", vcf_alt="*"),
        NoRow('ALT "." or "*"'),
    ),
    # VD-15. Two VCF records with the same CHROM, POS, REF and ALT, and the IDs var1 and var2. Each record gives a row
    # with its own ID as variant_id ("Output columns").
    Case(
        "two_identical_records_with_different_ids_give_a_row_each",
        """
        ref 5' [gacc ATG GAT G]|[TA AGC TAA gc] 3'
        alt 5' [gacc ATG GAT G]|[TA AGA TAA gc] 3'
                                      ^ C>A
                                      ^ C>A var2
        C>A: AGC>AGA, in two records: var1 and var2.
        """,
        TWO_EXONS,
        Change("AAG[C>A]TAA"),
        {
            **MISSENSE_ALT_CDS,
            "alt_transcript_seq": "GACCATGGATGTAAGATAAGC",
            "alt_transcript_length": 21,
            "alt_cds_start_in_transcript": 4,
            "alt_transcript_exons": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        more_changes=(Change("AAG[C>A]TAA", vcf_id="var2"),),
        more_rows=({"variant_id": "var2"},),
    ),
]
