"""
Conformance cases of transcripts without start_codon or stop_codon rows, the GENCODE tags
cds_start_NF and cds_end_NF, the Ensembl flavor with its codons from the FASTA, the mitochondrial stop codons, and the
code paths that differ between the strands. A transcript with the strand "." is an error.

The drawings show the transcript 5' to 3', also for the minus strand. The runner renders their layout block from the
case. The line ref is the transcript and the line alt is the transcript with the change. The two lines are aligned,
and `-` fills a gap. `[...]` is an exon and `|` is an exon junction. Upper case is the coding region, in codons of the
annotated frame; lower case is UTR. The coding region is the CDS with its stop codon, as the GFF3 CDS rows hold it.
`..N..` leaves out N bases. `^` marks the change as ref>alt. A mark line puts a character under bases, e.g. `***`
under the PTC, and `<-- label -->` spans a length. A ruler gives tx or CDS positions, as labelled. The comments
above the layouts use the same notation.
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

CDS_START_NF = ("cds_start_NF",)
CDS_END_NF = ("cds_end_NF",)
# The columns of a row that is neither a PTC nor a stop loss, and that has no scan of the alt transcript
NO_FLAGS = {"start_loss": False, "stop_loss": False, **NOT_SCANNED, "unknown_reason": None}


def _set_ensembl_end_phases(rows, phases):
    """Adds ensembl_end_phase to each exon row; ``phases`` maps the exon rank to the phase."""
    for row in rows:
        if row[2] == "exon":
            rank = int(row[8].split("rank=")[1])
            row[8] += f";ensembl_end_phase={phases[rank]}"
    return rows


def exon_1_ends_mid_codon(rows, strand):
    """Ensembl end phases: exon 1 ends mid-codon (2), exon 2 holds the 3' UTR (-1)."""
    return _set_ensembl_end_phases(rows, {1: 2, 2: -1})


def both_exons_end_mid_codon(rows, strand):
    """Ensembl end phases: exon 1 ends mid-codon (2), and so does exon 2, the last coding exon (1)."""
    return _set_ensembl_end_phases(rows, {1: 2, 2: 1})


def without_rank(rows, strand):
    """Removes the rank attribute of the Ensembl exon rows, so that the exon numbers come from the genomic order."""
    for row in rows:
        row[8] = ";".join(part for part in row[8].split(";") if not part.startswith("rank="))
    return rows


def strand_dot(rows, strand):
    """Set the strand of every GFF3 row to ".", which GFF3 uses for a feature without a strand."""
    return [[*row[:6], ".", *row[7:]] for row in rows]


# cds_start_NF: the CDS starts with ATG by chance, and the GFF3 has no start_codon rows
#
# tx      0   3        8             18   23
# ref 5' [ggg ATG AA]|[A CCC GAC TAA ggggg] 3'
NF_ATG = Layout(
    Transcript(("gggATGAA", "ACCCGACTAAggggg"), start_codon=False, tags=CDS_START_NF),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 15),
        "ref_cds_stop": per_strand(48, 50),
        "ref_cds_seq": "ATGAAACCCGACTAA",
        "ref_cds_len": 15,
        "has_start_codon": False,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 5), (2, 10)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": None,
        "ref_start_codon_exon": None,
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
        "likely_misannotated": True,
    },
)

# cds_start_NF: the CDS starts with CTG and has an in-frame ATG at CDS 3, an internal Met. The 3' UTR has a TAG in
# the frame of the CDS.
#
# tx      0   3   6   9   12      18                   66      72  75 78   83
# ref 5' [ggg CTG ATG AAG TGG GAC AAG AAG ..39.. AAG]|[CCC GAC TAA ggctagcc] 3'
# CDS         0   3   6   9       15                   63      69
NF_CTG_EXON_1 = "gggCTGATGAAGTGGGAC" + "AAG" * 16
NF_CTG = Layout(
    Transcript((NF_CTG_EXON_1, "CCCGACTAAggctagcc"), start_codon=False, tags=CDS_START_NF),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 18),
        "ref_cds_stop": per_strand(105, 110),
        "ref_cds_seq": "CTGATGAAGTGGGAC" + "AAG" * 16 + "CCCGACTAA",
        "ref_cds_len": 72,
        "has_start_codon": False,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 63), (2, 9)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": None,
        "ref_start_codon_exon": None,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 69,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(69, "TAA")],
        "ref_stop_codon_exons": [2],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 113,
        "transcript_seq": "GGGCTGATGAAGTGGGAC" + "AAG" * 16 + "CCCGACTAAGGCTAGCC",
        "transcript_length": 83,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 75,
        "transcript_exon_info": [(1, 66), (2, 17)],
        "utr3_length": 8,
        "utr5_length": 3,
        "total_exon_count": 2,
        "likely_misannotated": True,
    },
)

# cds_end_NF: the GFF3 has no stop_codon rows, and the CDS ends in the sense codon TGG at the transcript end
#
# tx      0   3             12         21
# ref 5' [gcc ATG GCC AAG]|[CTG GAC TGG] 3'
# CDS         0             9       15 18
CDS_END_NF_SENSE = Layout(
    Transcript(("gccATGGCCAAG", "CTGGACTGG"), stop_codon=False, tags=CDS_END_NF),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 10),
        "ref_cds_stop": per_strand(51, 48),
        "ref_cds_seq": "ATGGCCAAGCTGGACTGG",
        "ref_cds_len": 18,
        "has_start_codon": True,
        "has_stop_codon": False,
        "cds_frame": 0,
        "ref_cds_info": [(1, 9), (2, 9)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TGG",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 51,
        "transcript_seq": "GCCATGGCCAAGCTGGACTGG",
        "transcript_length": 21,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 21,
        "transcript_exon_info": [(1, 12), (2, 9)],
        "utr3_length": None,
        "utr5_length": 3,
        "total_exon_count": 2,
        "likely_misannotated": True,
    },
)

# The exons of the cds_end_NF trimming cases: the coding region ends in the complete in-frame stop codon TAA
#
# tx      0   3             12      18  21  25
# ref 5' [gcc ATG GCC AAG]|[CTG GAC TAA gccg] 3'
# CDS         0             9       15  18
TAA_EXONS = ("gccATGGCCAAG", "CTGGACTAAgccg")
# The columns of a layout with TAA_EXONS whose coding region keeps the TAA
TAA_KEPT = {
    **IDS,
    "ref_cds_start": per_strand(13, 14),
    "ref_cds_stop": per_strand(51, 52),
    "ref_cds_seq": "ATGGCCAAGCTGGACTAA",
    "ref_cds_len": 18,
    "has_start_codon": True,
    "has_stop_codon": True,
    "cds_frame": 0,
    "ref_cds_info": [(1, 9), (2, 9)],
    "cds_in_transcript": True,
    "ref_start_codon_pos": 0,
    "ref_start_codon_exon": 1,
    "ref_last_codon": "TAA",
    "ref_valid_stop": True,
    "ref_first_stop_codon": "TAA",
    "ref_first_stop_pos": 15,
    "ref_num_stop_codons": 1,
    "ref_all_stop_codons": [(15, "TAA")],
    "ref_stop_codon_exons": [2],
    "ref_is_premature": False,
    "transcript_start": 10,
    "transcript_end": 55,
    "transcript_seq": "GCCATGGCCAAGCTGGACTAAGCCG",
    "transcript_length": 25,
    "cds_start_in_transcript": 3,
    "cds_end_in_transcript": 21,
    "transcript_exon_info": [(1, 12), (2, 13)],
    "utr3_length": 4,
    "utr5_length": 3,
    "total_exon_count": 2,
    "likely_misannotated": False,
}
# cds_end_NF without stop_codon rows: the TAA leaves the coding region and becomes 3' UTR
CDS_END_NF_TRIMMED = Layout(
    Transcript(TAA_EXONS, stop_codon=False, tags=CDS_END_NF),
    {
        **TAA_KEPT,
        "ref_cds_start": per_strand(13, 17),
        "ref_cds_stop": per_strand(48, 52),
        "ref_cds_seq": "ATGGCCAAGCTGGAC",
        "ref_cds_len": 15,
        "has_stop_codon": False,
        "ref_cds_info": [(1, 9), (2, 6)],
        "ref_last_codon": "GAC",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "cds_end_in_transcript": 18,
        "utr3_length": None,
        "likely_misannotated": True,
    },
)
# cds_end_NF with stop_codon rows: the TAA stays the annotated stop codon
CDS_END_NF_WITH_STOP_CODON_ROWS = Layout(Transcript(TAA_EXONS, tags=CDS_END_NF), TAA_KEPT)
# No tag attribute in the GFF3 and no stop_codon rows: the TAA stays in the CDS, but it is no annotated stop codon
NO_TAG_NO_STOP_CODON_ROWS = Layout(
    Transcript(TAA_EXONS, stop_codon=False),
    {
        **TAA_KEPT,
        "has_stop_codon": False,
        "ref_valid_stop": False,
        "ref_is_premature": True,
        "utr3_length": None,
        "likely_misannotated": True,
    },
)
# Ensembl, tagged cds_start_NF and cds_end_NF on the mRNA row: the codons come from the FASTA, not from the tags
ENSEMBL_NF_TAGS = Layout(Transcript(TAA_EXONS, tags=CDS_START_NF + CDS_END_NF, flavor="ensembl"), TAA_KEPT)

# cds_end_NF without stop_codon rows: the stop codon TA|A is split across the intron. Its A in exon 2 is the only
# coding base of exon 2. Trimming removes the whole stop codon, and with it the CDS row of exon 2.
#
# tx      0   3               15   17    22
# ref 5' [gcc ATG GCC AAG CTG TA]|[A gccg] 3'
CDS_END_NF_SPLIT_STOP = Layout(
    Transcript(("gccATGGCCAAGCTGTA", "Agccg"), stop_codon=False, tags=CDS_END_NF),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 37),
        "ref_cds_stop": per_strand(25, 49),
        "ref_cds_seq": "ATGGCCAAGCTG",
        "ref_cds_len": 12,
        "has_start_codon": True,
        "has_stop_codon": False,
        "cds_frame": 0,
        "ref_cds_info": [(1, 12)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "CTG",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 52,
        "transcript_seq": "GCCATGGCCAAGCTGTAAGCCG",
        "transcript_length": 22,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 15,
        "transcript_exon_info": [(1, 17), (2, 5)],
        "utr3_length": None,
        "utr5_length": 3,
        "total_exon_count": 2,
        "likely_misannotated": True,
    },
)

# cds_end_NF: the last 3 CDS bases read TAA out of frame, so they stay in the CDS. The CDS runs to the transcript end.
#
# tx      0   3             12       19
# ref 5' [gcc ATG GCC AAG]|[CTG CTA A] 3'
# CDS         0             9    13  16
CDS_END_NF_TAA_OUT_OF_FRAME = Layout(
    Transcript(("gccATGGCCAAG", "CTGCTAA"), stop_codon=False, tags=CDS_END_NF),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 10),
        "ref_cds_stop": per_strand(49, 46),
        "ref_cds_seq": "ATGGCCAAGCTGCTAA",
        "ref_cds_len": 16,
        "has_start_codon": True,
        "has_stop_codon": False,
        "cds_frame": 0,
        "ref_cds_info": [(1, 9), (2, 7)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 49,
        "transcript_seq": "GCCATGGCCAAGCTGCTAA",
        "transcript_length": 19,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 19,
        "transcript_exon_info": [(1, 12), (2, 7)],
        "utr3_length": None,
        "utr5_length": 3,
        "total_exon_count": 2,
        "likely_misannotated": True,
    },
)

# cds_end_NF: the GFF3 has no stop_codon rows, and the CDS of 8 nt ends in TA. The 3' UTR starts with a, so the
# in-frame TAA starts in the CDS and ends past it.
#
# tx      0    4   7   10 12
# ref 5' [gacc ATG GCC TA agccaggc] 3'
# CDS          0   3   6  8
CDS_END_NF_TAA_ACROSS_THE_CDS_END = Layout(
    Transcript(("gaccATGGCCTAagccaggc",), stop_codon=False, tags=CDS_END_NF),
    {
        **IDS,
        "ref_cds_start": per_strand(14, 18),
        "ref_cds_stop": per_strand(22, 26),
        "ref_cds_seq": "ATGGCCTA",
        "ref_cds_len": 8,
        "has_start_codon": True,
        "has_stop_codon": False,
        "cds_frame": 0,
        "ref_cds_info": [(1, 8)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "CTA",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 30,
        "transcript_seq": "GACCATGGCCTAAGCCAGGC",
        "transcript_length": 20,
        "cds_start_in_transcript": 4,
        "cds_end_in_transcript": 12,
        "transcript_exon_info": [(1, 20)],
        "utr3_length": None,
        "utr5_length": 4,
        "total_exon_count": 1,
        "likely_misannotated": True,
    },
)

# cds_start_NF and cds_end_NF: the CDS covers the whole transcript. Its first base belongs to no complete codon
# (phase 1), and its codons TGG AAG | CTG GAC AAG follow. The CDS row of exon 2 has phase 0.
#
# tx      0 1         7          16
# ref 5' [C TGG AAG]|[CTG GAC AAG] 3'
NF_BOTH = Layout(
    Transcript(("CTGGAAG", "CTGGACAAG"), start_codon=False, stop_codon=False, frame=1, tags=CDS_START_NF + CDS_END_NF),
    {
        **IDS,
        "ref_cds_start": per_strand(10, 10),
        "ref_cds_stop": per_strand(46, 46),
        "ref_cds_seq": "CTGGAAGCTGGACAAG",
        "ref_cds_len": 16,
        "has_start_codon": False,
        "has_stop_codon": False,
        "cds_frame": 1,
        "ref_cds_info": [(1, 7), (2, 9)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": None,
        "ref_start_codon_exon": None,
        "ref_last_codon": "AAG",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 46,
        "transcript_seq": "CTGGAAGCTGGACAAG",
        "transcript_length": 16,
        "cds_start_in_transcript": 0,
        "cds_end_in_transcript": 16,
        "transcript_exon_info": [(1, 7), (2, 9)],
        "utr3_length": None,
        "utr5_length": 0,
        "total_exon_count": 2,
        "likely_misannotated": True,
    },
)

# One exon whose CDS ends in the codon XXX, for the stop codons of the mitochondrial chromosome
#
# tx      0   3               15  18  22
# ref 5' [gcc ATG GCC AAG CTG XXX gccg] 3'
# CDS         0               12  15


def one_exon(last_codon):
    return ("gccATGGCCAAGCTG" + last_codon + "gccg",)


# The columns of a one-exon layout whose CDS ends in a codon that the codon scans do not know as a stop codon
ONE_EXON = {
    **IDS,
    "ref_cds_start": per_strand(13, 14),
    "ref_cds_stop": per_strand(28, 29),
    "ref_cds_len": 15,
    "has_start_codon": True,
    "has_stop_codon": False,
    "cds_frame": 0,
    "ref_cds_info": [(1, 15)],
    "cds_in_transcript": True,
    "ref_start_codon_pos": 0,
    "ref_start_codon_exon": 1,
    "ref_valid_stop": False,
    "ref_first_stop_codon": None,
    "ref_first_stop_pos": None,
    "ref_num_stop_codons": 0,
    "ref_all_stop_codons": [],
    "ref_stop_codon_exons": [],
    "ref_is_premature": False,
    "transcript_start": 10,
    "transcript_end": 32,
    "transcript_length": 22,
    "cds_start_in_transcript": 3,
    "cds_end_in_transcript": 18,
    "transcript_exon_info": [(1, 22)],
    "utr3_length": None,
    "utr5_length": 3,
    "total_exon_count": 1,
    "likely_misannotated": True,
}
AGA = {
    **ONE_EXON,
    "ref_cds_seq": "ATGGCCAAGCTGAGA",
    "ref_last_codon": "AGA",
    "transcript_seq": "GCCATGGCCAAGCTGAGAGCCG",
}
ENSEMBL_AGA_ON_CHR1 = Layout(Transcript(one_exon("AGA"), flavor="ensembl"), AGA)
GENCODE_AGA_ON_CHRM = Layout(
    Transcript(one_exon("AGA"), stop_codon=False, contig="chrM"), {**AGA, "chromosome": "chrM"}
)
ENSEMBL_AGA_ON_MT = Layout(
    Transcript(one_exon("AGA"), flavor="ensembl", contig="MT"),
    {**AGA, "chromosome": "MT", "has_stop_codon": True, "utr3_length": 4},
)
ENSEMBL_AGG_ON_MT = Layout(
    Transcript(one_exon("AGG"), flavor="ensembl", contig="MT"),
    {
        **ONE_EXON,
        "chromosome": "MT",
        "ref_cds_seq": "ATGGCCAAGCTGAGG",
        "has_stop_codon": True,
        "ref_last_codon": "AGG",
        "transcript_seq": "GCCATGGCCAAGCTGAGGGCCG",
        "utr3_length": 4,
    },
)
ENSEMBL_TGA_ON_MT = Layout(
    Transcript(one_exon("TGA"), flavor="ensembl", contig="MT"),
    {
        **ONE_EXON,
        "chromosome": "MT",
        "ref_cds_seq": "ATGGCCAAGCTGTGA",
        "ref_last_codon": "TGA",
        "ref_first_stop_codon": "TGA",
        "ref_first_stop_pos": 12,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(12, "TGA")],
        "ref_stop_codon_exons": [1],
        "ref_is_premature": True,
        "transcript_seq": "GCCATGGCCAAGCTGTGAGCCG",
    },
)

# Ensembl: a CDS of 16 nt whose last 3 bases read TGA out of frame. Exon 1 ends mid-codon (ensembl_end_phase 2), and
# exon 2 holds the 3' UTR.
#
# tx      0   3                14      19  23
# ref 5' [gcc ATG GCC AAG CT]|[G CTG A gccg] 3'
# CDS         0                11 13   16
ENSEMBL_TGA_OUT_OF_FRAME = Layout(
    Transcript(("gccATGGCCAAGCT", "GCTGAgccg"), flavor="ensembl", edit_gff3=exon_1_ends_mid_codon),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 14),
        "ref_cds_stop": per_strand(49, 50),
        "ref_cds_seq": "ATGGCCAAGCTGCTGA",
        "ref_cds_len": 16,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 11), (2, 5)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TGA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 53,
        "transcript_seq": "GCCATGGCCAAGCTGCTGAGCCG",
        "transcript_length": 23,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 19,
        "transcript_exon_info": [(1, 14), (2, 9)],
        "utr3_length": 4,
        "utr5_length": 3,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)

# Ensembl: a CDS of 19 nt that runs to the transcript end and whose last 3 bases read TAA. Exon 1 ends mid-codon
# (ensembl_end_phase 2), and so does exon 2, the last coding exon (ensembl_end_phase 1).
#
# tx      0   3                14         22
# ref 5' [gcc ATG GCC AAG CT]|[G CTG CTA A] 3'
# CDS         0                11     16  19
ENSEMBL_LAST_EXON_ENDS_MID_CODON = Layout(
    Transcript(("gccATGGCCAAGCT", "GCTGCTAA"), flavor="ensembl", edit_gff3=both_exons_end_mid_codon),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 10),
        "ref_cds_stop": per_strand(52, 49),
        "ref_cds_seq": "ATGGCCAAGCTGCTGCTAA",
        "ref_cds_len": 19,
        "has_start_codon": True,
        "has_stop_codon": False,
        "cds_frame": 0,
        "ref_cds_info": [(1, 11), (2, 8)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": False,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": 0,
        "ref_all_stop_codons": [],
        "ref_stop_codon_exons": [],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 52,
        "transcript_seq": "GCCATGGCCAAGCTGCTGCTAA",
        "transcript_length": 22,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 22,
        "transcript_exon_info": [(1, 14), (2, 8)],
        "utr3_length": None,
        "utr5_length": 3,
        "total_exon_count": 2,
        "likely_misannotated": True,
    },
)

# Ensembl without rank attributes: the exon numbers come from the genomic order
#
# tx      0   3         9           16     21   26
# ref 5' [gcc ATG GCC]|[AAG CTG G]|[AC TAA gccgg] 3'
# CDS         0         6           13 15  18
ENSEMBL_WITHOUT_RANK = Layout(
    Transcript(("gccATGGCC", "AAGCTGG", "ACTAAgccgg"), flavor="ensembl", edit_gff3=without_rank),
    {
        **IDS,
        "ref_cds_start": per_strand(13, 15),
        "ref_cds_stop": per_strand(71, 73),
        "ref_cds_seq": "ATGGCCAAGCTGGACTAA",
        "ref_cds_len": 18,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 6), (2, 7), (3, 5)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 15,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(15, "TAA")],
        "ref_stop_codon_exons": [3],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 76,
        "transcript_seq": "GCCATGGCCAAGCTGGACTAAGCCGG",
        "transcript_length": 26,
        "cds_start_in_transcript": 3,
        "cds_end_in_transcript": 21,
        "transcript_exon_info": [(1, 9), (2, 7), (3, 10)],
        "utr3_length": 5,
        "utr5_length": 3,
        "total_exon_count": 3,
        "likely_misannotated": False,
    },
)

# The edges of the coding region: the 5' UTR ends in A before the start codon, and the 3' UTR starts with AAA after
# the stop codon TAA
#
# tx      0      6             15      21  24     31
# ref 5' [ggacca ATG GCC AAG]|[CTG CTG TAA aaagccg] 3'
# CDS            0             9       15  18
EDGES = Layout(
    Transcript(("ggaccaATGGCCAAG", "CTGCTGTAAaaagccg")),
    {
        **IDS,
        "ref_cds_start": per_strand(16, 17),
        "ref_cds_stop": per_strand(54, 55),
        "ref_cds_seq": "ATGGCCAAGCTGCTGTAA",
        "ref_cds_len": 18,
        "has_start_codon": True,
        "has_stop_codon": True,
        "cds_frame": 0,
        "ref_cds_info": [(1, 9), (2, 9)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": 0,
        "ref_start_codon_exon": 1,
        "ref_last_codon": "TAA",
        "ref_valid_stop": True,
        "ref_first_stop_codon": "TAA",
        "ref_first_stop_pos": 15,
        "ref_num_stop_codons": 1,
        "ref_all_stop_codons": [(15, "TAA")],
        "ref_stop_codon_exons": [2],
        "ref_is_premature": False,
        "transcript_start": 10,
        "transcript_end": 61,
        "transcript_seq": "GGACCAATGGCCAAGCTGCTGTAAAAAGCCG",
        "transcript_length": 31,
        "cds_start_in_transcript": 6,
        "cds_end_in_transcript": 24,
        "transcript_exon_info": [(1, 15), (2, 16)],
        "utr3_length": 7,
        "utr5_length": 6,
        "total_exon_count": 2,
        "likely_misannotated": False,
    },
)


# Ensembl: a CDS of 2 nt, too short for a codon check
#
# tx      0     5  7   11
# ref 5' [gccgc AG gccg] 3'
ENSEMBL_CDS_OF_2_NT = Layout(
    Transcript(("gccgcAGgccg",), flavor="ensembl"),
    {
        **IDS,
        "ref_cds_start": per_strand(15, 14),
        "ref_cds_stop": per_strand(17, 16),
        "ref_cds_seq": "AG",
        "ref_cds_len": 2,
        "has_start_codon": False,
        "has_stop_codon": False,
        "cds_frame": 0,
        "ref_cds_info": [(1, 2)],
        "cds_in_transcript": True,
        "ref_start_codon_pos": None,
        "ref_start_codon_exon": None,
        "ref_last_codon": None,
        "ref_valid_stop": None,
        "ref_first_stop_codon": None,
        "ref_first_stop_pos": None,
        "ref_num_stop_codons": None,
        "ref_all_stop_codons": None,
        "ref_stop_codon_exons": None,
        "ref_is_premature": None,
        "transcript_start": 10,
        "transcript_end": 21,
        "transcript_seq": "GCCGCAGGCCG",
        "transcript_length": 11,
        "cds_start_in_transcript": 5,
        "cds_end_in_transcript": 7,
        "transcript_exon_info": [(1, 11)],
        "utr3_length": None,
        "utr5_length": 5,
        "total_exon_count": 1,
        "likely_misannotated": True,
    },
)


# The alt CDS columns of the GCC>GAC missense at CDS 4 in a layout with TAA_EXONS whose coding region keeps the TAA
TAA_KEPT_MISSENSE = {
    "alt_cds_start": per_strand(13, 14),
    "alt_cds_stop": per_strand(51, 52),
    "alt_cds_seq": "ATGGACAAGCTGGACTAA",
    "alt_cds_len": 18,
    "alt_cds_info": [(1, 9), (2, 9)],
    "alt_start_codon_pos": 0,
    "alt_start_codon_exon": 1,
    "alt_last_codon": "TAA",
    "alt_first_stop_codon": "TAA",
    "alt_first_stop_pos": 15,
    "alt_num_stop_codons": 1,
    "alt_all_stop_codons": [(15, "TAA")],
    "alt_stop_codon_exons": [2],
    "alt_transcript_seq": "GCCATGGACAAGCTGGACTAAGCCG",
    "alt_transcript_length": 25,
    "alt_cds_start_in_transcript": 3,
}
# The record of the GCC>GAC missense at CDS 4 in a layout with TAA_EXONS
TAA_EXONS_MISSENSE_RECORD = {
    "variant_id": "var1",
    "ref": per_strand("C", "G"),
    "alt": per_strand("A", "T"),
    "start_variant": per_strand(17, 47),
    "end_variant": per_strand(18, 48),
}
# The alt CDS columns of the GCC>GAC missense at CDS 4 in a one-exon layout (one_exon)
ONE_EXON_MISSENSE = {
    "variant_id": "var1",
    "ref": per_strand("C", "G"),
    "alt": per_strand("A", "T"),
    "start_variant": per_strand(17, 24),
    "end_variant": per_strand(18, 25),
    "alt_cds_start": per_strand(13, 14),
    "alt_cds_stop": per_strand(28, 29),
    "alt_cds_len": 15,
    "alt_cds_info": [(1, 15)],
    "alt_start_codon_pos": 0,
    "alt_start_codon_exon": 1,
    "alt_valid_stop": False,
    "alt_first_stop_codon": None,
    "alt_first_stop_pos": None,
    "alt_num_stop_codons": 0,
    "alt_all_stop_codons": [],
    "alt_stop_codon_exons": [],
    "alt_is_premature": False,
    "alt_transcript_length": 22,
    "alt_cds_start_in_transcript": 3,
    **NO_FLAGS,
    **NO_PTC_FEATURES,
    "stop_codon_distance": None,
    **NO_RULE,
}
# The GCC>GAC missense in a one-exon layout whose CDS ends in AGA
AGA_MISSENSE = {
    **ONE_EXON_MISSENSE,
    "alt_cds_seq": "ATGGACAAGCTGAGA",
    "alt_last_codon": "AGA",
    "alt_transcript_seq": "GCCATGGACAAGCTGAGAGCCG",
}

CASES = [
    # NA-01, NF-02, NA-10
    Case(
        "atg_to_acg_at_the_start_of_a_cds_start_nf_cds_is_no_start_loss",
        """
        cds_start_NF: the GFF3 has no start_codon rows, so the ATG at CDS 0 is no annotated start codon. ATG>ACG
        changes it, but the transcript has no start codon to lose.

        tx      0   3        8             18   23
        ref 5' [ggg ATG AA]|[A CCC GAC TAA ggggg] 3'
        alt 5' [ggg ACG AA]|[A CCC GAC TAA ggggg] 3'
                     ^ T>C
        """,
        NF_ATG,
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
            "alt_transcript_seq": "GGGACGAAACCCGACTAAGGGGG",
            "alt_transcript_length": 23,
            "alt_cds_start_in_transcript": 3,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 3, 8, 18, 23)),
    ),
    # NF-03, NA-10
    Case(
        "lost_internal_atg_of_a_cds_start_nf_cds_is_no_start_loss",
        """
        cds_start_NF: the CDS starts with CTG, and its in-frame ATG at CDS 3 is an internal Met. ATG>ACG removes that
        ATG, but the transcript has no start codon to lose.

        CDS         0   3                        63      69
        ref 5' [ggg CTG ATG AAG TGG ..48.. AAG]|[CCC GAC TAA ggctagcc] 3'
        alt 5' [ggg CTG ACG AAG TGG ..48.. AAG]|[CCC GAC TAA ggctagcc] 3'
                         ^ T>C
        """,
        NF_CTG,
        Change("CTGA[T>C]GAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "start_variant": per_strand(17, 105),
            "end_variant": per_strand(18, 106),
            "alt_cds_start": per_strand(13, 18),
            "alt_cds_stop": per_strand(105, 110),
            "alt_cds_seq": "CTGACGAAGTGGGAC" + "AAG" * 16 + "CCCGACTAA",
            "alt_cds_len": 72,
            "alt_cds_info": [(1, 63), (2, 9)],
            "alt_start_codon_pos": None,
            "alt_start_codon_exon": None,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 69,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(69, "TAA")],
            "alt_stop_codon_exons": [2],
            "alt_is_premature": False,
            "alt_transcript_seq": "GGGCTGACGAAGTGGGAC" + "AAG" * 16 + "CCCGACTAAGGCTAGCC",
            "alt_transcript_length": 83,
            "alt_cds_start_in_transcript": 3,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 3, 63, 69), "CDS"),
    ),
    # NF-04, NA-10
    Case(
        "ptc_of_a_cds_start_nf_cds_has_no_distance_to_the_start_codon",
        """
        cds_start_NF: TGG>TAG at CDS 9 is a PTC 9 nt downstream of the CDS start and 6 nt downstream of the internal
        Met at CDS 3. The true start codon lies upstream of the CDS, at an unknown distance, so ptc_to_start_codon is
        null and the start-proximal rule is False. The PTC lies 54 nt upstream of the only exon junction.

        CDS         0   3       9
        ref 5' [ggg CTG ATG AAG TGG GAC AAG ..39.. AAG AAG]|[CCC GAC TAA ggctagcc] 3'
        alt 5' [ggg CTG ATG AAG TAG GAC AAG ..39.. AAG AAG]|[CCC GAC TAA ggctagcc] 3'
                                 ^ G>A
                                *** PTC
                                <-- ptc_to_intron = 54 -->
        """,
        NF_CTG,
        Change("GAAGT[G>A]GGAC"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(23, 99),
            "end_variant": per_strand(24, 100),
            "alt_cds_start": per_strand(13, 18),
            "alt_cds_stop": per_strand(105, 110),
            "alt_cds_seq": "CTGATGAAGTAGGAC" + "AAG" * 16 + "CCCGACTAA",
            "alt_cds_len": 72,
            "alt_cds_info": [(1, 63), (2, 9)],
            "alt_start_codon_pos": None,
            "alt_start_codon_exon": None,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 9,
            "alt_num_stop_codons": 2,
            "alt_all_stop_codons": [(9, "TAG"), (69, "TAA")],
            "alt_stop_codon_exons": [1, 2],
            "alt_is_premature": True,
            "alt_transcript_seq": "GGGCTGATGAAGTAGGAC" + "AAG" * 16 + "CCCGACTAAGGCTAGCC",
            "alt_transcript_length": 83,
            "alt_cds_start_in_transcript": 3,
            **NO_FLAGS,
            "upstream_exon_count": 0,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": None,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 66,
            "stop_codon_distance": 60,
            "ptc_to_intron": 54,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_annotated_start",
        },
        marks=(Mark("alt", 12, 15, "*", "PTC"), Span("alt", 12, 66, "ptc_to_intron = 54")),
        ruler=Ruler((0, 3, 9), "CDS"),
    ),
    # NF-05, NA-10, STR-08
    Case(
        "stop_loss_in_a_cds_start_nf_cds_reads_through_in_the_cds_frame_without_a_start_codon",
        """
        cds_start_NF: TAA>CAA loses the stop codon. The scan reads on from the CDS start at tx 3, in the frame of the
        CDS, to the TAG at tx 78 in the 3' UTR. It has no start codon: not the CTG at tx 3, and not the internal Met
        at tx 6. On the minus strand, the 3' UTR of 8 nt lies on the genomic left; a scan from tx 8 would read
        another frame.

        tx      0   3   6                    66      72     78   83
        ref 5' [ggg CTG ATG AAG ..51.. AAG]|[CCC GAC TAA ggctagcc] 3'
        alt 5' [ggg CTG ATG AAG ..51.. AAG]|[CCC GAC CAA ggctagcc] 3'
                                                     ^ T>C
        """,
        NF_CTG,
        Change("CCCGAC[T>C]AAGGC"),
        {
            "variant_id": "var1",
            "ref": per_strand("T", "A"),
            "alt": per_strand("C", "G"),
            "start_variant": per_strand(102, 20),
            "end_variant": per_strand(103, 21),
            "alt_cds_start": per_strand(13, 18),
            "alt_cds_stop": per_strand(105, 110),
            "alt_cds_seq": "CTGATGAAGTGGGAC" + "AAG" * 16 + "CCCGACCAA",
            "alt_cds_len": 72,
            "alt_cds_info": [(1, 63), (2, 9)],
            "alt_start_codon_pos": None,
            "alt_start_codon_exon": None,
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
            "alt_transcript_seq": "GGGCTGATGAAGTGGGAC" + "AAG" * 16 + "CCCGACCAAGGCTAGCC",
            "alt_transcript_length": 83,
            "alt_cds_start_in_transcript": 3,
            "transcript_start_codon_pos": None,
            "transcript_start_codon_exon": None,
            "transcript_last_codon": "GCC",
            "transcript_valid_stop": False,
            "transcript_first_stop_codon": "TAG",
            "transcript_first_stop_pos": 78,
            "transcript_num_stop_codons": 1,
            "transcript_all_stop_codons": [(78, "TAG")],
            "transcript_stop_codon_exons": [2],
            "unknown_reason": None,
            **NO_PTC_FEATURES,
            "stop_codon_distance": -6,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 3, 6, 66, 72, 78, 83)),
    ),
    # NA-03, NA-09, NF-13
    Case(
        "stop_gained_in_the_last_codon_of_a_cds_end_nf_cds_is_a_ptc",
        """
        cds_end_NF: the GFF3 has no stop_codon rows, and the CDS ends in the sense codon TGG at the transcript end.
        Without an annotated stop codon, a stop codon inside the alt CDS is a PTC: TGG>TAG in the last codon.

        tx      0   3             12      18 21
        ref 5' [gcc ATG GCC AAG]|[CTG GAC TGG] 3'
        alt 5' [gcc ATG GCC AAG]|[CTG GAC TAG] 3'
                                           ^ G>A
                                          *** PTC
        """,
        CDS_END_NF_SENSE,
        Change("GACT[G>A]GCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(49, 11),
            "end_variant": per_strand(50, 12),
            "alt_cds_start": per_strand(13, 10),
            "alt_cds_stop": per_strand(51, 48),
            "alt_cds_seq": "ATGGCCAAGCTGGACTAG",
            "alt_cds_len": 18,
            "alt_cds_info": [(1, 9), (2, 9)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAG",
            "alt_valid_stop": False,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 15,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(15, "TAG")],
            "alt_stop_codon_exons": [2],
            "alt_is_premature": True,
            "alt_transcript_seq": "GCCATGGCCAAGCTGGACTAG",
            "alt_transcript_length": 21,
            "alt_cds_start_in_transcript": 3,
            **NO_FLAGS,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 15,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 9,
            "stop_codon_distance": None,
            "ptc_to_intron": 3,
            **NO_RULE,
            "nmd_last_exon_rule": True,
            "nmd_start_proximal_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_annotated_stop",
        },
        marks=(Mark("alt", 18, 21, "*", "PTC"),),
        ruler=Ruler((0, 3, 12, 18, 21)),
    ),
    # Pins "Without an annotated stop codon (`has_stop_codon` False), a first stop codon inside the alt CDS is a PTC,
    # and one past its end is neither" for a stop codon across the CDS end (alt_is_premature: "Without an annotated
    # stop codon: whether it lies inside the alt CDS"). Its first 2 nt lie in the alt CDS and its last nt lies past
    # the end, so it does not lie inside the alt CDS and is no PTC. In the closest cases, the stop codon lies wholly
    # inside the alt CDS (NA-03, a PTC) or wholly past its end (ST-23 in cases_stop, neither).
    Case(
        "stop_codon_that_straddles_the_end_of_a_cds_end_nf_cds_is_no_ptc",
        """
        cds_end_NF: the GFF3 has no stop_codon rows, so has_stop_codon is False. The CDS of 8 nt ends at tx 12 in
        TA, and the 3' UTR starts with a. C>A changes GCC to GAC. The first in-frame stop codon of the alt
        transcript is the TAA at tx 10 to 12 (`*`). Its first 2 nt lie in the alt CDS, and its last nt lies past
        the end. So it does not lie inside the alt CDS: alt_is_premature is False, and stop_loss is False.

        tx      0    4          12      20
        ref 5' [gacc ATG GCC TA agccaggc] 3'
        alt 5' [gacc ATG GAC TA agccaggc] 3'
                          ^ C>A
                             **** first in-frame stop codon
        """,
        CDS_END_NF_TAA_ACROSS_THE_CDS_END,
        Change("ATGG[C>A]CTA"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(18, 21),
            "end_variant": per_strand(19, 22),
            "alt_cds_start": per_strand(14, 18),
            "alt_cds_stop": per_strand(22, 26),
            "alt_cds_seq": "ATGGACTA",
            "alt_cds_len": 8,
            "alt_cds_info": [(1, 8)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "CTA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": 0,
            "alt_all_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_is_premature": False,
            "alt_transcript_seq": "GACCATGGACTAAGCCAGGC",
            "alt_transcript_length": 20,
            "alt_cds_start_in_transcript": 4,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        marks=(Mark("alt", 10, 13, "*", "first in-frame stop codon"),),
        ruler=Ruler((0, 4, 12, 20)),
    ),
    # NA-09
    Case(
        "missense_in_the_last_codon_of_a_cds_end_nf_cds_is_no_stop_loss",
        """
        cds_end_NF: the GFF3 has no stop_codon rows, and the CDS ends in the sense codon TGG at the transcript end.
        TGG>TGC changes the last codon. Without an annotated stop codon, there is no stop codon to lose.

        tx      0   3             12      18 21
        ref 5' [gcc ATG GCC AAG]|[CTG GAC TGG] 3'
        alt 5' [gcc ATG GCC AAG]|[CTG GAC TGC] 3'
                                            ^ G>C
        """,
        CDS_END_NF_SENSE,
        Change("GACTG[G>C]CC"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("C", "G"),
            "start_variant": per_strand(50, 10),
            "end_variant": per_strand(51, 11),
            "alt_cds_start": per_strand(13, 10),
            "alt_cds_stop": per_strand(51, 48),
            "alt_cds_seq": "ATGGCCAAGCTGGACTGC",
            "alt_cds_len": 18,
            "alt_cds_info": [(1, 9), (2, 9)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TGC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": 0,
            "alt_all_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_is_premature": False,
            "alt_transcript_seq": "GCCATGGCCAAGCTGGACTGC",
            "alt_transcript_length": 21,
            "alt_cds_start_in_transcript": 3,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 3, 12, 18, 21)),
    ),
    # NF-07
    Case(
        "cds_end_nf_coding_region_without_stop_codon_rows_loses_its_complete_stop_codon",
        """
        cds_end_NF: the GFF3 has no stop_codon rows, and its CDS rows end in the complete in-frame stop codon TAA.
        These 3 nt leave the CDS and become 3' UTR. GCC>GAC is a missense. The TAA lies past the end of the alt CDS,
        so it is neither a PTC nor the annotated stop codon.

        tx      0   3             12      18  21  25
        ref 5' [gcc ATG GCC AAG]|[CTG GAC TAA gccg] 3'
        alt 5' [gcc ATG GAC AAG]|[CTG GAC TAA gccg] 3'
                         ^ C>A
        """,
        CDS_END_NF_TRIMMED,
        Change("ATGG[C>A]CAAG"),
        {
            **TAA_EXONS_MISSENSE_RECORD,
            **TAA_KEPT_MISSENSE,
            "alt_cds_start": per_strand(13, 17),
            "alt_cds_stop": per_strand(48, 52),
            "alt_cds_seq": "ATGGACAAGCTGGAC",
            "alt_cds_len": 15,
            "alt_cds_info": [(1, 9), (2, 6)],
            "alt_last_codon": "GAC",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": 0,
            "alt_all_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_is_premature": False,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 3, 12, 18, 21, 25)),
    ),
    # NF-08
    Case(
        "cds_end_nf_split_stop_codon_leaves_the_cds_with_the_cds_row_of_its_last_exon",
        """
        cds_end_NF: the GFF3 has no stop_codon rows, and its CDS rows end in the complete in-frame stop codon TA|A,
        split across the intron. The A is the only coding base of exon 2. The stop codon leaves the CDS, so exon 2 has
        no CDS row any more, and exon 1 loses TA. GCC>GAC is a missense.

        tx      0   3               15   17    22
        ref 5' [gcc ATG GCC AAG CTG TA]|[A gccg] 3'
        alt 5' [gcc ATG GAC AAG CTG TA]|[A gccg] 3'
                         ^ C>A
        """,
        CDS_END_NF_SPLIT_STOP,
        Change("ATGG[C>A]CAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(17, 44),
            "end_variant": per_strand(18, 45),
            "alt_cds_start": per_strand(13, 37),
            "alt_cds_stop": per_strand(25, 49),
            "alt_cds_seq": "ATGGACAAGCTG",
            "alt_cds_len": 12,
            "alt_cds_info": [(1, 12)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "CTG",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": 0,
            "alt_all_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_is_premature": False,
            "alt_transcript_seq": "GCCATGGACAAGCTGTAAGCCG",
            "alt_transcript_length": 22,
            "alt_cds_start_in_transcript": 3,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 3, 15, 17, 22)),
    ),
    # NF-09
    Case(
        "cds_end_nf_cds_that_ends_in_a_stop_codon_out_of_frame_keeps_it",
        """
        cds_end_NF: the GFF3 has no stop_codon rows. The last 3 CDS bases read TAA, but out of frame: the CDS of 16 nt
        holds 5 codons and one more base. So they stay in the CDS. The CDS runs to the transcript end. GCC>GAC is a
        missense, and the alt transcript has no in-frame stop codon.

        CDS         0             9    13  16
        ref 5' [gcc ATG GCC AAG]|[CTG CTA A] 3'
        alt 5' [gcc ATG GAC AAG]|[CTG CTA A] 3'
                         ^ C>A
        """,
        CDS_END_NF_TAA_OUT_OF_FRAME,
        Change("ATGG[C>A]CAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(17, 41),
            "end_variant": per_strand(18, 42),
            "alt_cds_start": per_strand(13, 10),
            "alt_cds_stop": per_strand(49, 46),
            "alt_cds_seq": "ATGGACAAGCTGCTAA",
            "alt_cds_len": 16,
            "alt_cds_info": [(1, 9), (2, 7)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": 0,
            "alt_all_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_is_premature": False,
            "alt_transcript_seq": "GCCATGGACAAGCTGCTAA",
            "alt_transcript_length": 19,
            "alt_cds_start_in_transcript": 3,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 9, 13, 16), "CDS"),
    ),
    # NF-10
    Case(
        "cds_end_nf_with_stop_codon_rows_keeps_its_stop_codon",
        """
        cds_end_NF with stop_codon rows on the TAA: the TAA stays the annotated stop codon. GCC>GAC is a missense.

        tx      0   3             12      18  21  25
        ref 5' [gcc ATG GCC AAG]|[CTG GAC TAA gccg] 3'
        alt 5' [gcc ATG GAC AAG]|[CTG GAC TAA gccg] 3'
                         ^ C>A
        """,
        CDS_END_NF_WITH_STOP_CODON_ROWS,
        Change("ATGG[C>A]CAAG"),
        {
            **TAA_EXONS_MISSENSE_RECORD,
            **TAA_KEPT_MISSENSE,
            "alt_valid_stop": True,
            "alt_is_premature": False,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 3, 12, 18, 21, 25)),
    ),
    # NF-11, NA-09
    Case(
        "missense_without_stop_codon_rows_and_tags_is_a_ptc_at_the_last_codon",
        """
        The GFF3 has neither stop_codon rows nor a tag attribute, so the CDS keeps the TAA (`*`) at its end, and it is
        no annotated stop codon. Without an annotated stop codon, every in-frame stop codon is premature, also this
        TAA. So the missense GCC>GAC gives a PTC row: the TAA lies inside the alt CDS.

        tx      0   3             12      18  21  25
        ref 5' [gcc ATG GCC AAG]|[CTG GAC TAA gccg] 3'
        alt 5' [gcc ATG GAC AAG]|[CTG GAC TAA gccg] 3'
                         ^ C>A
                                          *** PTC
                                          <------> ptc_to_intron = 7
        """,
        NO_TAG_NO_STOP_CODON_ROWS,
        Change("ATGG[C>A]CAAG"),
        {
            **TAA_EXONS_MISSENSE_RECORD,
            **TAA_KEPT_MISSENSE,
            "alt_valid_stop": False,
            "alt_is_premature": True,
            **NO_FLAGS,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 15,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 13,
            "stop_codon_distance": None,
            "ptc_to_intron": 7,
            **NO_RULE,
            "nmd_last_exon_rule": True,
            "nmd_start_proximal_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ref_ptc",
        },
        marks=(Mark("alt", 18, 21, "*", "PTC"), Span("alt", 18, 25, "ptc_to_intron = 7")),
        ruler=Ruler((0, 3, 12, 18, 21, 25)),
    ),
    # NF-15, NF-13, NF-14, NF-01 (phase 1), STR-03
    Case(
        "ptc_in_a_cds_start_nf_and_cds_end_nf_cds_with_phase_1",
        """
        cds_start_NF and cds_end_NF: the GFF3 has neither start_codon nor stop_codon rows, and the CDS covers the
        whole transcript. Its first base belongs to no complete codon (phase 1 on the CDS row of exon 1; the CDS row
        of exon 2 has phase 0). TGG>TAG at CDS 1 is an in-frame PTC, 6 nt upstream of the only exon junction. Read in
        phase 0, there is no stop codon.

        tx      0 1         7          16
        ref 5' [C TGG AAG]|[CTG GAC AAG] 3'
        alt 5' [C TAG AAG]|[CTG GAC AAG] 3'
                   ^ G>A
                  *** PTC
                  <-----> ptc_to_intron = 6
        """,
        NF_BOTH,
        Change("CT[G>A]GAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("G", "C"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(12, 43),
            "end_variant": per_strand(13, 44),
            "alt_cds_start": per_strand(10, 10),
            "alt_cds_stop": per_strand(46, 46),
            "alt_cds_seq": "CTAGAAGCTGGACAAG",
            "alt_cds_len": 16,
            "alt_cds_info": [(1, 7), (2, 9)],
            "alt_start_codon_pos": None,
            "alt_start_codon_exon": None,
            "alt_last_codon": "AAG",
            "alt_valid_stop": False,
            "alt_first_stop_codon": "TAG",
            "alt_first_stop_pos": 1,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(1, "TAG")],
            "alt_stop_codon_exons": [1],
            "alt_is_premature": True,
            "alt_transcript_seq": "CTAGAAGCTGGACAAG",
            "alt_transcript_length": 16,
            "alt_cds_start_in_transcript": 0,
            **NO_FLAGS,
            "upstream_exon_count": 0,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": None,
            "ptc_less_than_150nt_to_start": False,
            "ptc_exon_length": 7,
            "stop_codon_distance": None,
            "ptc_to_intron": 6,
            **NO_RULE,
            "nmd_50nt_penultimate_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_annotated_stop",
        },
        marks=(Mark("alt", 1, 4, "*", "PTC"), Span("alt", 1, 7, "ptc_to_intron = 6")),
        ruler=Ruler((0, 1, 7, 16)),
    ),
    # NF-06, NF-12, NA-05 (TAA)
    Case(
        "ensembl_cds_start_nf_and_cds_end_nf_tags_do_not_change_the_codons_from_the_fasta",
        """
        Ensembl: the GFF3 has no start_codon or stop_codon rows, and the mRNA row has the tags cds_start_NF and
        cds_end_NF. The CDS starts with ATG in phase 0 and ends in TAA, so it has a start and a stop codon, and the
        TAA stays in the CDS. GCC>GAC is a missense.

        tx      0   3             12      18  21  25
        ref 5' [gcc ATG GCC AAG]|[CTG GAC TAA gccg] 3'
        alt 5' [gcc ATG GAC AAG]|[CTG GAC TAA gccg] 3'
                         ^ C>A
        """,
        ENSEMBL_NF_TAGS,
        Change("ATGG[C>A]CAAG"),
        {
            **TAA_EXONS_MISSENSE_RECORD,
            **TAA_KEPT_MISSENSE,
            "alt_valid_stop": True,
            "alt_is_premature": False,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 3, 12, 18, 21, 25)),
    ),
    # NA-05 (AGA on chr1)
    Case(
        "ensembl_cds_ending_in_aga_on_chr1_has_no_stop_codon",
        """
        Ensembl: the last 3 CDS bases AGA are no stop codon outside the mitochondrial chromosome. GCC>GAC is a
        missense.

        tx      0   3               15  18  22
        ref 5' [gcc ATG GCC AAG CTG AGA gccg] 3'
        alt 5' [gcc ATG GAC AAG CTG AGA gccg] 3'
                         ^ C>A
        """,
        ENSEMBL_AGA_ON_CHR1,
        Change("ATGG[C>A]CAAG"),
        {**AGA_MISSENSE, "alt_transcript_exon_info": SAME_EXONS, "nmd_model_status": "no_ptc"},
        ruler=Ruler((0, 3, 15, 18, 22)),
    ),
    # NA-06 (AGA)
    Case(
        "ensembl_cds_ending_in_aga_on_mt_has_a_stop_codon_that_the_codon_scans_do_not_know",
        """
        Ensembl on MT: the last 3 CDS bases AGA are a stop codon of the mitochondrial code, so the CDS has an
        annotated stop codon. The codon scans know only TAA, TAG and TGA, so ref_valid_stop is False, the transcript
        does not stop at the annotated stop codon, and the row keeps the flags from the CDS. GCC>GAC is a missense.

        tx      0   3               15  18  22
        ref 5' [gcc ATG GCC AAG CTG AGA gccg] 3'
        alt 5' [gcc ATG GAC AAG CTG AGA gccg] 3'
                         ^ C>A
        """,
        ENSEMBL_AGA_ON_MT,
        Change("ATGG[C>A]CAAG"),
        {**AGA_MISSENSE, "alt_transcript_exon_info": SAME_EXONS, "nmd_model_status": "no_ptc"},
        ruler=Ruler((0, 3, 15, 18, 22)),
    ),
    # NA-06 (AGG)
    Case(
        "ensembl_cds_ending_in_agg_on_mt_has_a_stop_codon_that_the_codon_scans_do_not_know",
        """
        Ensembl on MT: the last 3 CDS bases AGG are a stop codon of the mitochondrial code, as AGA is. GCC>GAC is a
        missense.

        tx      0   3               15  18  22
        ref 5' [gcc ATG GCC AAG CTG AGG gccg] 3'
        alt 5' [gcc ATG GAC AAG CTG AGG gccg] 3'
                         ^ C>A
        """,
        ENSEMBL_AGG_ON_MT,
        Change("ATGG[C>A]CAAG"),
        {
            **ONE_EXON_MISSENSE,
            "alt_cds_seq": "ATGGACAAGCTGAGG",
            "alt_last_codon": "AGG",
            "alt_transcript_seq": "GCCATGGACAAGCTGAGGGCCG",
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 3, 15, 18, 22)),
    ),
    # NA-06 (TGA)
    Case(
        "ensembl_cds_ending_in_tga_on_mt_has_no_stop_codon_and_its_tga_is_a_ptc",
        """
        Ensembl on MT: TGA is no stop codon of the mitochondrial code, so the CDS has no annotated stop codon. The
        codon scans know TGA as a stop codon, and without an annotated stop codon, every in-frame stop codon is
        premature. So the missense GCC>GAC gives a PTC row: the TGA (`*`) lies inside the alt CDS.

        tx      0   3               15  18  22
        ref 5' [gcc ATG GCC AAG CTG TGA gccg] 3'
        alt 5' [gcc ATG GAC AAG CTG TGA gccg] 3'
                         ^ C>A
                                    *** PTC
                                    <------> ptc_to_intron = 7
        """,
        ENSEMBL_TGA_ON_MT,
        Change("ATGG[C>A]CAAG"),
        {
            **ONE_EXON_MISSENSE,
            "alt_cds_seq": "ATGGACAAGCTGTGA",
            "alt_last_codon": "TGA",
            "alt_first_stop_codon": "TGA",
            "alt_first_stop_pos": 12,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(12, "TGA")],
            "alt_stop_codon_exons": [1],
            "alt_is_premature": True,
            "alt_transcript_seq": "GCCATGGACAAGCTGTGAGCCG",
            "upstream_exon_count": 0,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 12,
            "ptc_less_than_150nt_to_start": True,
            "ptc_exon_length": 22,
            "ptc_to_intron": 7,
            "nmd_last_exon_rule": True,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": True,
            "nmd_escape": True,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "ref_ptc",
        },
        marks=(Mark("alt", 15, 18, "*", "PTC"), Span("alt", 15, 22, "ptc_to_intron = 7")),
        ruler=Ruler((0, 3, 15, 18, 22)),
    ),
    # NA-07
    Case(
        "gencode_cds_ending_in_aga_on_chrm_without_stop_codon_rows_has_no_stop_codon",
        """
        GENCODE on chrM without stop_codon rows: the CDS has no annotated stop codon, although its last 3 bases AGA
        are a stop codon of the mitochondrial code. GCC>GAC is a missense.

        tx      0   3               15  18  22
        ref 5' [gcc ATG GCC AAG CTG AGA gccg] 3'
        alt 5' [gcc ATG GAC AAG CTG AGA gccg] 3'
                         ^ C>A
        """,
        GENCODE_AGA_ON_CHRM,
        Change("ATGG[C>A]CAAG"),
        {**AGA_MISSENSE, "alt_transcript_exon_info": SAME_EXONS, "nmd_model_status": "no_ptc"},
        ruler=Ruler((0, 3, 15, 18, 22)),
    ),
    # NA-05 (TGA out of frame)
    Case(
        "ensembl_cds_ending_in_tga_out_of_frame_has_a_stop_codon",
        """
        Ensembl: the last 3 CDS bases read TGA, so the CDS has an annotated stop codon, although the CDS of 16 nt
        holds 5 codons and one more base. The frame is not checked. Exon 1 ends mid-codon (ensembl_end_phase 2), which
        does not matter, since it is not the last coding exon. The transcript does not stop at the TGA, so the row
        keeps the flags from the CDS. GCC>GAC is a missense.

        CDS         0                11 13   16
        ref 5' [gcc ATG GCC AAG CT]|[G CTG A gccg] 3'
        alt 5' [gcc ATG GAC AAG CT]|[G CTG A gccg] 3'
                         ^ C>A
        """,
        ENSEMBL_TGA_OUT_OF_FRAME,
        Change("ATGG[C>A]CAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(17, 45),
            "end_variant": per_strand(18, 46),
            "alt_cds_start": per_strand(13, 14),
            "alt_cds_stop": per_strand(49, 50),
            "alt_cds_seq": "ATGGACAAGCTGCTGA",
            "alt_cds_len": 16,
            "alt_cds_info": [(1, 11), (2, 5)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TGA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": 0,
            "alt_all_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_is_premature": False,
            "alt_transcript_seq": "GCCATGGACAAGCTGCTGAGCCG",
            "alt_transcript_length": 23,
            "alt_cds_start_in_transcript": 3,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 11, 13, 16), "CDS"),
    ),
    # NA-05 (ensembl_end_phase 1 on the last coding exon)
    Case(
        "ensembl_cds_whose_last_coding_exon_ends_mid_codon_has_no_stop_codon",
        """
        Ensembl: the last 3 CDS bases read TAA, but the last coding exon ends mid-codon (ensembl_end_phase 1), so the
        CDS has no annotated stop codon. The CDS runs to the transcript end. GCC>GAC is a missense, and the alt
        transcript has no in-frame stop codon.

        CDS         0                11     16  19
        ref 5' [gcc ATG GCC AAG CT]|[G CTG CTA A] 3'
        alt 5' [gcc ATG GAC AAG CT]|[G CTG CTA A] 3'
                         ^ C>A
        """,
        ENSEMBL_LAST_EXON_ENDS_MID_CODON,
        Change("ATGG[C>A]CAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(17, 44),
            "end_variant": per_strand(18, 45),
            "alt_cds_start": per_strand(13, 10),
            "alt_cds_stop": per_strand(52, 49),
            "alt_cds_seq": "ATGGACAAGCTGCTGCTAA",
            "alt_cds_len": 19,
            "alt_cds_info": [(1, 11), (2, 8)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": False,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": 0,
            "alt_all_stop_codons": [],
            "alt_stop_codon_exons": [],
            "alt_is_premature": False,
            "alt_transcript_seq": "GCCATGGACAAGCTGCTGCTAA",
            "alt_transcript_length": 22,
            "alt_cds_start_in_transcript": 3,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 11, 16, 19), "CDS"),
    ),
    # NA-08 (fewer than 3 CDS bases)
    Case(
        "ensembl_cds_of_2_nt_has_neither_start_nor_stop_codon",
        """
        Ensembl: the CDS has 2 nt, too few for a codon check, so it has neither a start nor a stop codon. The codon
        columns of both CDS are null. AG>CG changes the first CDS base. Without an annotated stop codon, the row does
        not keep the flags from the CDS, and the alt transcript has no in-frame stop codon.

        tx      0     5  7   11
        ref 5' [gccgc AG gccg] 3'
        alt 5' [gccgc CG gccg] 3'
                      ^ A>C
        """,
        ENSEMBL_CDS_OF_2_NT,
        Change("gccgc[A>C]Ggccg"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("C", "G"),
            "start_variant": per_strand(15, 15),
            "end_variant": per_strand(16, 16),
            "alt_cds_start": per_strand(15, 14),
            "alt_cds_stop": per_strand(17, 16),
            "alt_cds_seq": "CG",
            "alt_cds_len": 2,
            "alt_cds_info": [(1, 2)],
            "alt_start_codon_pos": None,
            "alt_start_codon_exon": None,
            "alt_last_codon": None,
            "alt_valid_stop": None,
            "alt_first_stop_codon": None,
            "alt_first_stop_pos": None,
            "alt_num_stop_codons": None,
            "alt_all_stop_codons": None,
            "alt_stop_codon_exons": None,
            "alt_is_premature": False,
            "alt_transcript_seq": "GCCGCCGGCCG",
            "alt_transcript_length": 11,
            "alt_cds_start_in_transcript": 5,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": None,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 5, 7, 11)),
    ),
    # STR-10
    Case(
        "ensembl_exon_numbers_without_rank_follow_the_direction_of_transcription",
        """
        Ensembl without rank attributes: the exon numbers come from the genomic order, in the direction of
        transcription. On the minus strand, exon 1 has the largest Start. GCC>GAC is a missense.

        tx      0   3         9           16     21   26
        ref 5' [gcc ATG GCC]|[AAG CTG G]|[AC TAA gccgg] 3'
        alt 5' [gcc ATG GAC]|[AAG CTG G]|[AC TAA gccgg] 3'
                         ^ C>A
        """,
        ENSEMBL_WITHOUT_RANK,
        Change("ATGG[C>A]CGTAAG"),
        {
            "variant_id": "var1",
            "ref": per_strand("C", "G"),
            "alt": per_strand("A", "T"),
            "start_variant": per_strand(17, 68),
            "end_variant": per_strand(18, 69),
            "alt_cds_start": per_strand(13, 15),
            "alt_cds_stop": per_strand(71, 73),
            "alt_cds_seq": "ATGGACAAGCTGGACTAA",
            "alt_cds_len": 18,
            "alt_cds_info": [(1, 6), (2, 7), (3, 5)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 15,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(15, "TAA")],
            "alt_stop_codon_exons": [3],
            "alt_is_premature": False,
            "alt_transcript_seq": "GCCATGGACAAGCTGGACTAAGCCGG",
            "alt_transcript_length": 26,
            "alt_cds_start_in_transcript": 3,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": SAME_EXONS,
            "nmd_model_status": "no_ptc",
        },
        ruler=Ruler((0, 3, 9, 16, 21, 26)),
    ),
    # STR-06
    Case(
        "deletion_of_an_a_before_the_start_codon_shortens_the_5utr",
        """
        Deleting one A of the run aA in ccaATG: only the placement in the 5' UTR keeps an ATG at the start, so the
        5' UTR loses the A, and the CDS stays as it is. On the minus strand, the 5' UTR lies on the genomic right.

        tx      0      6             15      21  24     31
        ref 5' [ggacca ATG GCC AAG]|[CTG CTG TAA aaagccg] 3'
        alt 5' [ggacc- ATG GCC AAG]|[CTG CTG TAA aaagccg] 3'
                     ^ a>-
        """,
        EDGES,
        Change("ggacc[a>]ATGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("CA", "TT"),
            "alt": per_strand("C", "T"),
            "start_variant": per_strand(14, 54),
            "end_variant": per_strand(16, 56),
            "alt_cds_start": per_strand(16, 17),
            "alt_cds_stop": per_strand(54, 55),
            "alt_cds_seq": "ATGGCCAAGCTGCTGTAA",
            "alt_cds_len": 18,
            "alt_cds_info": [(1, 9), (2, 9)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 15,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(15, "TAA")],
            "alt_stop_codon_exons": [2],
            "alt_is_premature": False,
            "alt_transcript_seq": "GGACCATGGCCAAGCTGCTGTAAAAAGCCG",
            "alt_transcript_length": 30,
            "alt_cds_start_in_transcript": 5,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": [(1, 14), (2, 16)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("ggacca[A>]TGGCC"),),
        ruler=Ruler((0, 6, 15, 21, 24, 31)),
    ),
    # STR-05, STR-09
    Case(
        "deletion_of_an_a_in_the_run_after_the_stop_codon_shortens_the_3utr",
        """
        Deleting one A of the run AAaaa in TAAaaa: the stop codon edge takes the placement shifted farthest into the
        3' UTR, so the 3' UTR loses an A, and the CDS stays as it is. The alt line draws the 5'-most placement, in the
        stop codon. On the minus strand, the left-normalized VCF record lies in the 3' UTR, and its other placements
        reach into the stop codon.

        tx      0      6             15      21  24     31
        ref 5' [ggacca ATG GCC AAG]|[CTG CTG TAA aaagccg] 3'
        alt 5' [ggacca ATG GCC AAG]|[CTG CTG T-A aaagccg] 3'
                                              ^ A>-
        """,
        EDGES,
        Change("CTGCTGT[A>]AAAAGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("TA", "TT"),
            "alt": per_strand("T", "T"),
            "start_variant": per_strand(51, 17),
            "end_variant": per_strand(53, 19),
            "alt_cds_start": per_strand(16, 17),
            "alt_cds_stop": per_strand(54, 55),
            "alt_cds_seq": "ATGGCCAAGCTGCTGTAA",
            "alt_cds_len": 18,
            "alt_cds_info": [(1, 9), (2, 9)],
            "alt_start_codon_pos": 0,
            "alt_start_codon_exon": 1,
            "alt_last_codon": "TAA",
            "alt_valid_stop": True,
            "alt_first_stop_codon": "TAA",
            "alt_first_stop_pos": 15,
            "alt_num_stop_codons": 1,
            "alt_all_stop_codons": [(15, "TAA")],
            "alt_stop_codon_exons": [2],
            "alt_is_premature": False,
            "alt_transcript_seq": "GGACCAATGGCCAAGCTGCTGTAAAAGCCG",
            "alt_transcript_length": 30,
            "alt_cds_start_in_transcript": 6,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": [(1, 15), (2, 15)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("CTGCTGTAAAA[A>]GCC"),),
        ruler=Ruler((0, 6, 15, 21, 24, 31)),
    ),
    # STR-04
    Case(
        "insertion_that_repeats_the_start_codon_adds_a_met",
        """
        Inserting ATG at the start codon gives ATGATG. A scanning ribosome starts at the first ATG, so the coding
        region starts at the inserted ATG and gains a Met. The alt line draws the insertion before the annotated
        start codon, so in lower case. The insertion has 4 placements, from before the ATG to after it. On the minus
        strand, the first ATG is the 3'-most CAT in the genome.

        tx      0         6             15      21  24     31
        ref 5' [ggacca ---ATG GCC AAG]|[CTG CTG TAA aaagccg] 3'
        alt 5' [ggacca atgATG GCC AAG]|[CTG CTG TAA aaagccg] 3'
                       ^^^ ->ATG
        """,
        EDGES,
        Change("ggacca[>ATG]ATGGCC"),
        {
            "variant_id": "var1",
            "ref": per_strand("A", "T"),
            "alt": per_strand("AATG", "TCAT"),
            "start_variant": per_strand(15, 54),
            "end_variant": per_strand(16, 55),
            "alt_cds_start": per_strand(16, 17),
            "alt_cds_stop": per_strand(54, 55),
            "alt_cds_seq": "ATGATGGCCAAGCTGCTGTAA",
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
            "alt_transcript_seq": "GGACCAATGATGGCCAAGCTGCTGTAAAAAGCCG",
            "alt_transcript_length": 34,
            "alt_cds_start_in_transcript": 6,
            **NO_FLAGS,
            **NO_PTC_FEATURES,
            "stop_codon_distance": 0,
            **NO_RULE,
            "alt_transcript_exon_info": [(1, 18), (2, 16)],
            "nmd_model_status": "no_ptc",
        },
        equivalent=(Change("ggaccaATG[>ATG]GCCAAG"), Change("ggaccaA[>TGA]TGGCC")),
        ruler=Ruler((0, 6, 15, 21, 24, 31)),
    ),
    # STR-11
    Case(
        "transcript_with_strand_dot_is_an_error_that_names_the_transcript",
        """
        Every GFF3 row of the transcript has the strand "." instead of + or -. annotate() rejects the GFF3 with a
        ValueError that names the transcript tx1.

        ref 5' [ggacca ATG GCC AAG]|[CTG CTG TAA aaagccg] 3'
        alt 5' [ggacca ATG GAC AAG]|[CTG CTG TAA aaagccg] 3'
                            ^ C>A
        """,
        Layout(Transcript(EDGES.transcript.exons, edit_gff3=strand_dot), EDGES.ref),
        Change("ATGG[C>A]CAAG"),
        Raises(ValueError, match="tx1"),
    ),
]
