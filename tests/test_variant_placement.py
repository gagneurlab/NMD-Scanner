"""
The layout drawings show a synthetic genome in transcript orientation, 5' to 3', and are not to scale. `.` is a flank
or an intron, `[...]` is an exon, `u` is UTR, `=` is CDS and `s` is the stop codon. The numbers under a layout are
layout positions: 0-based positions in transcript orientation, also on the minus strand.
"""

import pandas as pd
import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

from nmd_scanner.extra_features import add_nmd_features, evaluate_nmd_escape_rules
from nmd_scanner.rules import extract_ptc
from nmd_scanner.variant_placement import (
    EXON_BOUNDARY_AMBIGUOUS,
    SPLICE_SITE_DESTROYED,
    Placement,
    ReferenceSequence,
    equivalent_placements,
    place_in_transcript,
    trim_alleles,
    variant_placements,
)


class _Reference:
    """ReferenceSequence stand-in on a string."""

    def __init__(self, sequence):
        self.sequence = sequence

    def base(self, position):
        return self.sequence[position] if 0 <= position < len(self.sequence) else None

    def bases(self, start, end):
        return self.sequence[max(start, 0) : max(end, 0)]


def test_reference_sequence_keeps_at_most_max_chunks(tmp_path):
    # 4 chunks of 1024 bases. Reads across a chunk edge, and reads of an evicted chunk, give the right bases.
    sequence = "ACGTTGCA" * 512
    (tmp_path / "genome.fa").write_text(f">chrR\n{sequence.lower()}\n")
    reference = ReferenceSequence(Fasta(str(tmp_path / "genome.fa")), "chrR", max_chunks=2)
    for start, end in [(1020, 1030), (3000, 3100), (0, 5), (2040, 2050), (1020, 1030), (4090, 5000)]:
        assert reference.bases(start, end) == sequence[start:end]
        assert len(reference._chunks) <= 2
    assert reference.base(4095) == sequence[4095]
    assert reference.base(4096) is None


def test_trim_alleles():
    # VCF padding base of an indel
    assert trim_alleles(10, "CAG", "C") == (11, "AG", "")
    assert trim_alleles(10, "C", "CTT") == (11, "", "TT")
    # an MNV that changes only its first base
    assert trim_alleles(10, "TG", "CG") == (10, "T", "C")
    # shared bases on both sides
    assert trim_alleles(10, "ATCATG", "ATTTG") == (12, "CA", "T")
    assert trim_alleles(10, "A", "A") == (11, "", "")


def test_equivalent_placements_deletion_in_repeat():
    reference = _Reference("GGCTGTGTAC")
    # deleting TG at 3, GT at 4, TG at 5 or GT at 6 gives GGCTGTAC
    assert equivalent_placements(3, "TG", "", reference) == [Placement(s, s + 2, "") for s in [3, 4, 5, 6]]
    assert equivalent_placements(6, "GT", "", reference) == [Placement(s, s + 2, "") for s in [3, 4, 5, 6]]


def test_equivalent_placements_insertion_in_repeat():
    reference = _Reference("CCAGAGTT")
    # inserting AG anywhere in AGAG gives CCAGAGAGTT
    assert equivalent_placements(4, "", "AG", reference) == [
        Placement(2, 2, "AG"),
        Placement(3, 3, "GA"),
        Placement(4, 4, "AG"),
        Placement(5, 5, "GA"),
        Placement(6, 6, "AG"),
    ]


def test_equivalent_placements_at_chromosome_ends():
    reference = _Reference("AAAC")
    assert equivalent_placements(1, "A", "", reference) == [Placement(s, s + 1, "") for s in [0, 1, 2]]
    reference = _Reference("CAAA")
    assert equivalent_placements(1, "A", "", reference) == [Placement(s, s + 1, "") for s in [1, 2, 3]]


def test_equivalent_placements_substitution():
    reference = _Reference("AAAAAA")
    assert equivalent_placements(2, "AA", "CC", reference) == [Placement(2, 4, "CC")]
    # A delins whose REF and ALT differ in length is matched from the right and from the left
    assert equivalent_placements(2, "AA", "C", reference) == [
        Placement(2, 4, "C", match_left=False),
        Placement(2, 4, "C"),
    ]


def test_variant_placements_without_change():
    # REF and ALT differ only in case
    assert variant_placements(2, "g", "G", _Reference("ACGTACGT")) == []


# Genome for place_in_transcript, plus strand: exon [2, 10), intron, exon [30, 40)
#                   0         1         2         3
#                   0123456789012345678901234567890123456789
_GENOME = _Reference("CCATGGCTCTGTGTAAGCCCCCCTTTTCAGAGTGAACGTTGGCC")
_EXONS = [(2, 10), (30, 40)]


def _effect(variant_start, ref, alt, coding_rows=((2, 10), (30, 40)), exons=_EXONS):
    placements = variant_placements(variant_start, ref, alt, _GENOME)
    coding_region = (min(start for start, _ in coding_rows), max(end for _, end in coding_rows))
    return place_in_transcript(placements, list(coding_rows), exons, _GENOME, "+", coding_region)


def test_place_in_transcript_exonic_snv():
    effect = _effect(9, "T", "C")
    assert effect.unknown_reason is None
    assert effect.alt_coding == {(2, 10): "ATGGCTCC", (30, 40): "AGTGAACGTT"}


def test_place_in_transcript_untouched():
    # intron +5, and an intron variant that does not shift into the splice dinucleotide
    assert _effect(14, "A", "G") is None
    assert _effect(18, "C", "") is None


def test_place_in_transcript_splice_dinucleotide_snv():
    assert _effect(10, "G", "A").unknown_reason == SPLICE_SITE_DESTROYED
    assert _effect(29, "G", "C").unknown_reason == SPLICE_SITE_DESTROYED


def test_place_in_transcript_deletion_over_boundary_destroys_splice_site():
    # exon -1 to intron +5: TGTGTA leaves CTCAG..., no GT after any position
    assert _effect(9, "TGTGTA", "").unknown_reason == SPLICE_SITE_DESTROYED


def test_place_in_transcript_ambiguous_insertion():
    # inserting GT in CT|GTGT: the exon keeps its end or gains GT
    assert _effect(10, "", "GT").unknown_reason == EXON_BOUNDARY_AMBIGUOUS


def test_place_in_transcript_transcript_end_has_no_splice_site():
    # the coding region reaches the transcript start at 2; a deletion over it leaves the exon shorter
    effect = _effect(1, "CA", "")
    assert effect.unknown_reason is None
    assert effect.alt_coding[2, 10] == "TGGCTCT"


def test_place_in_transcript_changes_over_boundary_that_keep_the_donor():
    # CT|GTGT to CA|GTCT: a substitution maps base for base, and GT stays after the exon
    effect = _effect(9, "TGTG", "AGTC")
    assert effect.unknown_reason is None
    assert effect.alt_coding[2, 10] == "ATGGCTCA"
    # CT|GTGT to CTCC|GTAA: GT stays after the boundary if the delins is matched from the left, not from the right
    effect = _effect(9, "TGTG", "CGTAA")
    assert effect.unknown_reason is None
    assert effect.alt_coding[2, 10] == "ATGGCTCC"
    # TG to ATG trims to T to AT, which ends at the boundary
    effect = _effect(9, "TG", "ATG")
    assert effect.unknown_reason is None
    assert effect.alt_coding[2, 10] == "ATGGCTCAT"
    # CT|GT to CA|TT changes the donor
    assert _effect(9, "TG", "AT").unknown_reason == SPLICE_SITE_DESTROYED


def _start_codon_effect(strand, utr5, position, ref, alt):
    """
    place_in_transcript for a variant near the start codon of a one-exon transcript. In transcript orientation, the
    layout is: flank 0-9, exon 10-37 (5'UTR 10-15, coding row 16-27, 3'UTR 28-37), flank 38-47.

    :param utr5: the 6 bases of the 5'UTR
    :param position, ref, alt: the change in transcript orientation at a layout position; ref is empty for an
        insertion before ``position``
    :return: the alt bases of the coding row in transcript orientation
    """
    layout = "CCCCCCCCCC" + utr5 + "ATGCTGCTGTAA" + "GGCCGGCCGG" + "CCCCCCCCCC"
    assert layout[position : position + len(ref)] == ref
    length = len(layout)
    if strand == "+":
        genome, start, coding_row, exon = layout, position, (16, 28), (10, 38)
    else:
        genome = str(Seq(layout).reverse_complement())
        start = length - position - len(ref)
        ref, alt = str(Seq(ref).reverse_complement()), str(Seq(alt).reverse_complement())
        coding_row, exon = (length - 28, length - 16), (length - 38, length - 10)
    reference = _Reference(genome)
    placements = variant_placements(start, ref, alt, reference)
    effect = place_in_transcript(placements, [coding_row], [exon], reference, strand, coding_row)
    assert effect.unknown_reason is None
    alt_coding = effect.alt_coding[coding_row]
    return alt_coding if strand == "+" else str(Seq(alt_coding).reverse_complement())


@pytest.mark.parametrize(
    "utr5, position, ref, alt, alt_coding",
    [
        # ATG to ATGATG: the 5'-most placement of the insertion gives the first ATG, so the coding region gains a Met
        ("GGACCC", 19, "", "ATG", "ATGATGCTGCTGTAA"),
        # A|ATG minus one A: only the placement in the 5'UTR keeps an ATG at the start
        ("GGACCA", 16, "A", "", "ATGCTGCTGTAA"),
        # T|ATG minus TA: no placement keeps an ATG at the start, so the placement farthest into the 5'UTR counts
        ("GGACCT", 15, "TA", "", "TGCTGCTGTAA"),
    ],
    ids=["atg_duplication", "atg_on_the_utr_side", "no_atg"],
)
def test_place_in_transcript_start_codon(strand, utr5, position, ref, alt, alt_coding):
    assert _start_codon_effect(strand, utr5, position, ref, alt) == alt_coding


# Synthetic transcript in transcript orientation (5' to 3'), with layout positions:
# flank 0-9, exon 1 10-21 (5'UTR 10-13, CDS 14-21), intron 1 22-41, exon 2 42-57 (CDS), intron 2 58-77,
# exon 3 78-95 (CDS 78-83, stop codon 84-86, 3'UTR 87-95), flank 96-105
#
# 5' ....[uuuu========]....[================]....[======sssuuuuuuuuu].... 3'
#    0   10  14       22   42               58   78    84 87       96   106
_FLANK = "CCCCCCCCCC"
_UTR5 = "GACC"
_EXON1_CDS = "ATGGCTCT"
_INTRON1 = "GTGTAAGCCCCCCTTTTCAG"
_EXON2 = "AGTGAACGTTGGAAGC"
_INTRON2 = "GTAAGTCCCCCCTTTCCTAG"
_EXON3_CDS = "CTGCGT"
_STOP = "TAA"
_UTR3 = "AAAGCTGCC"
_LAYOUT = _FLANK + _UTR5 + _EXON1_CDS + _INTRON1 + _EXON2 + _INTRON2 + _EXON3_CDS + _STOP + _UTR3 + _FLANK
_REF_CDS = _EXON1_CDS + _EXON2 + _EXON3_CDS + _STOP
_DONOR1, _ACCEPTOR2 = 22, 42  # end of exon 1, start of exon 2
_STOP_END = 87
# (feature, exon number, start, end) in layout positions. A CDS row includes the stop codon.
_ROWS = [
    ("exon", 1, 10, 22),
    ("exon", 2, 42, 58),
    ("exon", 3, 78, 96),
    ("CDS", 1, 14, 22),
    ("CDS", 2, 42, 58),
    ("CDS", 3, 78, 87),
]

# A second transcript, whose coding region is exon 2: flank 0-9, exon 1 10-19 (5'UTR), intron 1 20-39,
# exon 2 40-63 (CDS 40-60, stop codon 61-63), intron 2 64-83, exon 3 84-95 (3'UTR), flank 96-105
#
# 5' ....[uuuuuuuuuu]....[=====================sss]....[uuuuuuuuuuuu].... 3'
#    0   10         20   40                    61 64   84           96   106
_EDGE_CDS = "ATGGCCAAGCTGCTGAAGCTGTAA"
_EDGE_INTRON = "GTAAGTCCCCCCTTTTTCAG"
_EDGE_LAYOUT = _FLANK + "GCCACCGCAG" + _EDGE_INTRON + _EDGE_CDS + _EDGE_INTRON + "GCCCCCCCCCCC" + _FLANK
_EDGE_ROWS = [
    ("exon", 1, 10, 20),
    ("exon", 2, 40, 64),
    ("exon", 3, 84, 96),
    ("CDS", 2, 40, 64),
]


def _run(tmp_path, strand, position, ref, alt, layout=_LAYOUT, rows=_ROWS):
    """
    Run extract_ptc for one variant in a synthetic transcript, by default the first one.

    :param position, ref, alt: the change in transcript orientation at a layout position; ref is empty for an
        insertion before ``position``. The VCF record is on the plus strand, with a padding base on its left for
        an indel.
    :param rows: (feature, exon number, start, end) of the exon and CDS rows in layout positions. The CDS rows
        include the stop codon.
    :return: the result rows (DataFrame)
    """
    assert layout[position : position + len(ref)] == ref
    length = len(layout)
    if strand == "+":
        genome = layout

        def to_genome(start, end):
            return start, end

    else:
        genome = str(Seq(layout).reverse_complement())

        def to_genome(start, end):
            return length - end, length - start

        ref, alt = str(Seq(ref).reverse_complement()), str(Seq(alt).reverse_complement())

    start, _ = to_genome(position, position + len(ref))
    if len(ref) != len(alt):
        start, ref, alt = start - 1, genome[start - 1] + ref, genome[start - 1] + alt

    chrom = f"chr_{tmp_path.name}"
    (tmp_path / "genome.fa").write_text(f">{chrom}\n{genome}\n")
    fasta = Fasta(str(tmp_path / "genome.fa"))

    annotation = pd.DataFrame(
        [
            {
                "Chromosome": chrom,
                "Start": to_genome(s, e)[0],
                "End": to_genome(s, e)[1],
                "Strand": strand,
                "Feature": feature,
                "exon_number": str(number),
                "transcript_id": "tx",
                "gene_id": "gene",
            }
            for feature, number, s, e in rows
        ]
    )
    vcf = pd.DataFrame(
        [{"Chromosome": chrom, "Start": start, "End": start + len(ref), "ID": "var", "Ref": ref, "Alt": alt}]
    )
    coding = annotation[annotation["Feature"] == "CDS"].assign(has_stop_codon=True)
    return extract_ptc(coding, vcf, fasta, annotation[annotation["Feature"] == "exon"])


def _single_row(results):
    assert len(results) == 1
    return results.iloc[0]


@pytest.fixture(params=["+", "-"])
def strand(request):
    return request.param


def test_exonic_snv(tmp_path, strand):
    row = _single_row(_run(tmp_path, strand, 15, "T", "C"))
    assert pd.isna(row["unknown_reason"])
    assert row["ref_cds_seq"] == _REF_CDS
    assert row["alt_cds_seq"] == "ACGGCTCT" + _REF_CDS[8:]


def test_mnv_over_donor_that_changes_only_the_exon_base(tmp_path, strand):
    # exon -1 T>C, intron +1 G stays
    row = _single_row(_run(tmp_path, strand, _DONOR1 - 1, "TG", "CG"))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == "ATGGCTCC" + _REF_CDS[8:]


def test_deletion_with_an_equivalent_intronic_placement(tmp_path, strand):
    # CT|GTGT: deleting TG over the boundary equals deleting the intron's GT, so the exon stays as it is
    row = _single_row(_run(tmp_path, strand, _DONOR1 - 1, "TG", ""))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == _REF_CDS


def test_acceptor_deletion_left_aligned_into_the_intron(tmp_path, strand):
    # CAG|AGT: deleting the intron's AG equals deleting the exon's AG; only the latter keeps the acceptor
    row = _single_row(_run(tmp_path, strand, _ACCEPTOR2 - 2, "AG", ""))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == _EXON1_CDS + _EXON2[2:] + _EXON3_CDS + _STOP


@pytest.mark.parametrize(
    "position, inserted, alt_cds",
    [
        (_ACCEPTOR2, "T", _EXON1_CDS + "T" + _EXON2 + _EXON3_CDS + _STOP),
        (_DONOR1, "A", _EXON1_CDS + "A" + _EXON2 + _EXON3_CDS + _STOP),
    ],
    ids=["acceptor", "donor"],
)
def test_insertion_between_intron_and_exon(tmp_path, strand, position, inserted, alt_cds):
    # The inserted bases go into the exon. On the plus strand, the acceptor lies on the left genomic side of the
    # exon and the donor on the right; on the minus strand it is the other way round.
    row = _single_row(_run(tmp_path, strand, position, "", inserted))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == alt_cds


@pytest.mark.parametrize("position, alt", [(_DONOR1, "A"), (_DONOR1 + 1, "C"), (_ACCEPTOR2 - 1, "C")])
def test_splice_dinucleotide_snv(tmp_path, strand, position, alt):
    row = _single_row(_run(tmp_path, strand, position, _LAYOUT[position], alt))
    assert row["unknown_reason"] == SPLICE_SITE_DESTROYED
    assert row["ref_cds_seq"] == _REF_CDS
    assert pd.isna(row["alt_cds_seq"])
    assert pd.isna(row["start_loss"]) and pd.isna(row["stop_loss"])

    row = row.to_dict()
    features = add_nmd_features(row)
    assert features["ptc_less_than_150nt_to_start"] is None
    assert features["total_exon_count"] == 3
    assert all(value is None for value in evaluate_nmd_escape_rules(row).values())


def test_intron_snv_outside_the_dinucleotide_has_no_row(tmp_path, strand):
    assert _run(tmp_path, strand, _DONOR1 + 2, "G", "A").empty


def test_ambiguous_acceptor_insertion(tmp_path, strand):
    # CAG|AG + AG: the exon may start at either AG
    row = _single_row(_run(tmp_path, strand, _ACCEPTOR2, "", "AG"))
    assert row["unknown_reason"] == EXON_BOUNDARY_AMBIGUOUS
    assert pd.isna(row["alt_cds_seq"])


def test_deletion_in_a_run_at_the_stop_codon(tmp_path, strand):
    # TAA|AAA: deleting one A anywhere in the run leaves TAA in place. On the minus strand, the left-normalized
    # VCF placement lies fully in the 3'UTR, so the VCF interval does not touch the coding region.
    for position in [_STOP_END - 2, _STOP_END + 2]:
        row = _single_row(_run(tmp_path, strand, position, "A", ""))
        assert pd.isna(row["unknown_reason"])
        assert row["alt_cds_seq"] == _REF_CDS
        assert row["stop_loss"] == False
        # The 3'UTR loses the A
        assert row["alt_transcript_seq"] == _UTR5 + _REF_CDS + _UTR3[1:]


def test_insertion_after_the_stop_codon_has_no_row(tmp_path, strand):
    assert _run(tmp_path, strand, _STOP_END, "", "C").empty


def test_insertion_before_the_start_codon_has_no_row(tmp_path, strand):
    # GACC|ATG to GACCC|ATG: the inserted base does not start an ATG, so it goes into the 5'UTR
    assert _run(tmp_path, strand, 14, "", "C").empty


def test_insertion_that_repeats_the_start_codon(tmp_path, strand):
    # ATG to ATGATG: a scanning ribosome starts at the first ATG, so the coding region gains a Met
    row = _single_row(_run(tmp_path, strand, 17, "", "ATG"))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == "ATG" + _REF_CDS
    assert row["start_loss"] == False


def test_deletion_in_the_start_codon(tmp_path, strand):
    # ATG to AG: no placement keeps an ATG at the start
    row = _single_row(_run(tmp_path, strand, 15, "T", ""))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == "AGGCTCT" + _REF_CDS[8:]
    assert row["start_loss"] == True


# In _EDGE_LAYOUT, the coding region starts at the start of exon 2 and its stop codon ends exon 2. An insertion at
# such an exon edge goes into the exon, as the splice site says, but into its UTR part, as the coding region edge
# rule says. So these insertions touch no coding base and give no row.
@pytest.mark.parametrize("position", [40, 64], ids=["before_the_start_codon", "after_the_stop_codon"])
def test_insertion_at_an_exon_edge_of_the_coding_region_has_no_row(tmp_path, strand, position):
    assert _run(tmp_path, strand, position, "", "CC", _EDGE_LAYOUT, _EDGE_ROWS).empty


def test_insertion_that_repeats_the_start_codon_at_an_exon_start(tmp_path, strand):
    # CAG|ATG to CAG|ATGATG: the 5'-most placement gives the first ATG in the exon, so the coding region gains a Met
    row = _single_row(_run(tmp_path, strand, 43, "", "ATG", _EDGE_LAYOUT, _EDGE_ROWS))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == "ATG" + _EDGE_CDS
    assert row["start_loss"] == False


def test_insertion_in_the_stop_codon_run_at_an_exon_end(tmp_path, strand):
    # TAA|GT to TAAAA|GT: the placement right after the stop codon puts the inserted bases into the 3'UTR part
    row = _single_row(_run(tmp_path, strand, 62, "", "AA", _EDGE_LAYOUT, _EDGE_ROWS))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == _EDGE_CDS
    assert row["stop_loss"] == False


def _donor_layout(exon1_end):
    """_LAYOUT with exon 1 ending in ``exon1_end`` and intron 1 starting with GTAAGT."""
    layout = _FLANK + _UTR5 + "ATGGC" + exon1_end + _EDGE_INTRON + _LAYOUT[_ACCEPTOR2:]
    assert len(layout) == len(_LAYOUT)
    return layout


# A delins whose REF and ALT differ in length is matched from the left and from the right. The splice site rule and
# the coding region edge rules treat these two like equivalent placements.
@pytest.mark.parametrize(
    "exon1_end, position, ref, alt",
    [
        # CTG|GTA to CGGTTTT: GT follows GG, but neither matching keeps GT after the boundary
        ("CTG", 20, "TGGTA", "GGTTTT"),
        # CAG|GTA to GGGTTTT: the same with GGG
        ("CAG", 19, "CAGGTA", "GGGTTTT"),
    ],
    ids=["TGGTA_GGTTTT", "CAGGTA_GGGTTTT"],
)
def test_delins_over_donor_that_destroys_it(tmp_path, strand, exon1_end, position, ref, alt):
    row = _single_row(_run(tmp_path, strand, position, ref, alt, _donor_layout(exon1_end)))
    assert row["unknown_reason"] == SPLICE_SITE_DESTROYED


@pytest.mark.parametrize(
    "alt, exon1_cds",
    [
        # CTCT|GTGTAA to CTCCC|GTATAA: only matching from the right keeps GT after the boundary
        ("CCGTA", "ATGGCTCCC"),
        # CTCT|GTGTAA to CTCA|GTCCTAA: only matching from the left keeps GT after the boundary
        ("AGTCC", "ATGGCTCA"),
    ],
    ids=["from_the_right", "from_the_left"],
)
def test_delins_over_donor_that_one_matching_keeps(tmp_path, strand, alt, exon1_cds):
    row = _single_row(_run(tmp_path, strand, _DONOR1 - 1, "TGTG", alt))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == exon1_cds + _REF_CDS[8:]


def test_delins_over_donor_that_both_matchings_keep(tmp_path, strand):
    # CT|GTGT to CTAGTGTC: GT follows A if matched from the left, and AGT if matched from the right
    row = _single_row(_run(tmp_path, strand, _DONOR1 - 1, "TGTG", "AGTGTC"))
    assert row["unknown_reason"] == EXON_BOUNDARY_AMBIGUOUS


def test_delins_over_the_start_codon_edge(tmp_path, strand):
    # GACC|ATGG to GACGATGTGG: only matching from the left puts an ATG at the start codon edge, so it starts there
    row = _single_row(_run(tmp_path, strand, 13, "CA", "GATG"))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == "ATG" + _REF_CDS[1:]
    # The G that replaces the last 5'UTR base stays in the 5'UTR
    assert row["alt_transcript_seq"] == "GACG" + "ATG" + _REF_CDS[1:] + _UTR3
    assert row["alt_cds_start_in_transcript"] == 4


def test_delins_over_the_stop_codon_edge(tmp_path, strand):
    # TA|A|AAA to TA|CCC|AAA: the matching with its length change in the 3'UTR counts, so C replaces the last stop
    # codon base and CC goes into the 3'UTR
    row = _single_row(_run(tmp_path, strand, _STOP_END - 1, "AA", "CCC"))
    assert row["alt_cds_seq"] == _REF_CDS[:-1] + "C"
    assert row["stop_loss"] == True
    assert row["alt_transcript_seq"] == _UTR5 + _REF_CDS[:-1] + "C" + "CC" + _UTR3[1:]


def test_deletion_that_shortens_the_5utr(tmp_path, strand):
    """
    GGACCA|ATG minus one A: only the placement in the 5'UTR keeps an ATG at the start, so the 5'UTR is 1 nt shorter
    and the coding region is unchanged. The alt CDS starts 1 nt earlier in the alt transcript.

    ref 5' ....[uuuuuu=========sssuuuuuuuuuu].... 3'
               10    16       25 28        38
       tx      0     6        15 18        28
    alt 5' ....[uuuuu=========sssuuuuuuuuuu]..... 3'
       tx      0    5        14 17        27
    """
    layout = _FLANK + "GGACCA" + "ATGCTGCTGTAA" + "GGCCGGCCGG" + _FLANK
    rows = [("exon", 1, 10, 38), ("CDS", 1, 16, 28)]
    row = _single_row(_run(tmp_path, strand, 16, "A", "", layout, rows))
    assert row["alt_cds_seq"] == "ATGCTGCTGTAA"
    assert row["start_loss"] == False
    assert row["alt_transcript_seq"] == "GGACC" + "ATGCTGCTGTAA" + "GGCCGGCCGG"
    assert row["cds_start_in_transcript"] == 6
    assert row["alt_cds_start_in_transcript"] == 5


def test_delins_over_the_stop_codon_and_the_donor(tmp_path, strand):
    # CTG TA|A|GTA to CTG TG|TAAGT|T: only matching from the right keeps the donor GT, 4 bases after the old exon end.
    # The coding region ends at the exon end, so the transcript holds every alt exon base.
    row = _single_row(_run(tmp_path, strand, 62, "AAGTA", "GTAAGTT", _EDGE_LAYOUT, _EDGE_ROWS))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == _EDGE_CDS[:-2] + "GTAA"
    assert row["alt_transcript_seq"] == "GCCACCGCAG" + _EDGE_CDS[:-2] + "GTAA" + "GCCCCCCCCCCC"


# A third transcript with a short coding exon 2: flank 0-9, exon 1 10-29 (5'UTR 10-14, CDS 15-29), intron 1 30-49,
# exon 2 50-61 (CDS), intron 2 62-81, exon 3 82-99 (CDS 82-87, stop codon 88-90, 3'UTR 91-99), flank 100-109
_SHORT_EXON1_CDS = "ATGGCCAAGCTGCTG"
_SHORT_EXON2 = "CAGCAGCTGCTG"
_SHORT_EXON3_CDS = "AAGCTGTAA"
#
# 5' ....[uuuuu===============]....[============]....[======sssuuuuuuuuu].... 3'
#    0   10   15              30   50           62   82    88 91       100  110
_SHORT_LAYOUT = (
    _FLANK + "GCCAC" + _SHORT_EXON1_CDS + _EDGE_INTRON + _SHORT_EXON2 + _EDGE_INTRON + "AAGCTGTAAGCCCCCCCC" + _FLANK
)
_SHORT_ROWS = [
    ("exon", 1, 10, 30),
    ("exon", 2, 50, 62),
    ("exon", 3, 82, 100),
    ("CDS", 1, 15, 30),
    ("CDS", 2, 50, 62),
    ("CDS", 3, 82, 91),
]


# A placement counts only if it keeps every splice site it maps. Taking each splice site from another placement
# can make two exons overlap.
@pytest.mark.parametrize(
    "position, ref, alt",
    [
        # G|intron 1|CA to 19 C, AG, CCC, GT, 20 C: matched from the left, only the acceptor keeps its AG. Matched
        # from the right, only the donor keeps its GT, 3 bases after that acceptor, so CCC would be in both exons.
        (29, "G" + _EDGE_INTRON + "CA", "C" * 19 + "AGCCCGT" + "C" * 20),
        # TCAG|exon 2|GTA to CCAGTC: matched from the left, only the acceptor keeps its AG. Matched from the right,
        # only the donor keeps its GT, 1 base before that acceptor, so exon 2 would have a negative length.
        (46, "TCAG" + _SHORT_EXON2 + "GTA", "CCAGTC"),
    ],
    ids=["overlapping_exons", "exon_of_negative_length"],
)
def test_delins_whose_matchings_each_keep_one_splice_site(tmp_path, strand, position, ref, alt):
    row = _single_row(_run(tmp_path, strand, position, ref, alt, _SHORT_LAYOUT, _SHORT_ROWS))
    assert row["unknown_reason"] == SPLICE_SITE_DESTROYED
    assert pd.isna(row["alt_cds_seq"])


def test_deletion_of_a_whole_short_exon(tmp_path, strand):
    # CAG|exon 2|GT to CAG|GT. Deleting G and exon 2 without its last G gives the same sequence but loses the
    # acceptor AG, so only the placement of exon 2 counts. It keeps both splice sites and leaves exon 2 empty.
    row = _single_row(_run(tmp_path, strand, 50, _SHORT_EXON2, "", _SHORT_LAYOUT, _SHORT_ROWS))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == _SHORT_EXON1_CDS + _SHORT_EXON3_CDS


def test_delins_over_a_short_exon_that_one_matching_keeps(tmp_path, strand):
    # CAG|exon 2|G to TAG|AGCAGCTGCTGC|GTC: matched from the left, both splice sites stay. Matched from the right,
    # the acceptor moves 2 bases to TAGAG|, but the donor loses its GT, so that matching does not count.
    alt = "TAGAGCAGCTGCTGCGTC"
    row = _single_row(_run(tmp_path, strand, 47, "CAG" + _SHORT_EXON2 + "G", alt, _SHORT_LAYOUT, _SHORT_ROWS))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == _SHORT_EXON1_CDS + "AGCAGCTGCTGC" + _SHORT_EXON3_CDS


def test_delins_over_the_start_codon_and_a_donor_uses_the_matching_that_keeps_the_donor(tmp_path, strand):
    # GACC|ATGGCTCT|GTGTAAG to 12 C, GTGCAG, 17 C: only matching from the left keeps the donor GT. The coding region
    # starts where that matching puts the start codon, so 8 C stay in exon 1. Matched from the right, the start
    # codon would lie after the donor.
    alt = "C" * 12 + "GTGCAG" + "C" * 17
    row = _single_row(_run(tmp_path, strand, 10, _UTR5 + _EXON1_CDS + _INTRON1[:7], alt))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == "C" * 8 + _REF_CDS[8:]


@pytest.mark.parametrize(
    "layout, rows, position, ref, alt, alt_cds, utr5, utr3",
    [
        # TAA|AA to CCC|TGAC
        (_LAYOUT, _ROWS, _STOP_END - 3, "TAAAA", "CCCTGAC", _REF_CDS[:-3] + "CCC", _UTR5, "TGAC" + _UTR3[2:]),
        # TAA|G to CCC|TGA
        (
            _SHORT_LAYOUT,
            _SHORT_ROWS,
            88,
            "TAAG",
            "CCCTGA",
            _SHORT_EXON1_CDS + _SHORT_EXON2 + "AAGCTG" + "CCC",
            "GCCAC",
            "TGA" + "CCCCCCCC",
        ),
    ],
    ids=["utr3_run", "utr3_g"],
)
def test_delins_from_the_stop_codon_into_the_3utr(
    tmp_path, strand, layout, rows, position, ref, alt, alt_cds, utr5, utr3
):
    # Matching from the left gives the coding region CCC in place of the stop codon, and the 3'UTR starts with TGA.
    # Read on in frame, TGA is the next codon: a stop loss with a stop codon right after the coding region.
    row = _single_row(_run(tmp_path, strand, position, ref, alt, layout, rows))
    assert row["alt_cds_seq"] == alt_cds
    assert row["alt_transcript_seq"] == utr5 + alt_cds + utr3
    assert row["stop_loss"] == True
    assert row["alt_is_premature"] == False
    assert row["transcript_first_stop_pos"] == len(utr5) + len(alt_cds)
