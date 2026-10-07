"""
The layout drawings show a synthetic genome in transcript orientation, 5' to 3', and are not to scale. `.` is a flank
or an intron, `[...]` is an exon, `u` is UTR, `=` is CDS and `s` is the stop codon. The numbers under a layout are
layout positions: 0-based positions in transcript orientation, also on the minus strand.
"""

import pandas as pd
import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

from nmd_scanner.rules import extract_ptc
from nmd_scanner.variant_placement import (
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
# (feature, exon number, start, end) in layout positions. A CDS row includes the stop codon.
_ROWS = [
    ("exon", 1, 10, 22),
    ("exon", 2, 42, 58),
    ("exon", 3, 78, 96),
    ("CDS", 1, 14, 22),
    ("CDS", 2, 42, 58),
    ("CDS", 3, 78, 87),
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

        ref = str(Seq(ref).reverse_complement())
        alt = str(Seq(alt).reverse_complement())

    start, _ = to_genome(position, position + len(ref))
    if len(ref) != len(alt):
        ref = genome[start - 1] + ref
        alt = genome[start - 1] + alt
        start -= 1

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
                "Frame": "0",
            }
            for feature, number, s, e in rows
        ]
    )
    vcf = pd.DataFrame(
        [{"Chromosome": chrom, "Start": start, "End": start + len(ref), "ID": "var", "Ref": ref, "Alt": alt}]
    )
    coding = annotation[annotation["Feature"] == "CDS"].assign(has_start_codon=True, has_stop_codon=True)
    return extract_ptc(coding, vcf, fasta, annotation[annotation["Feature"] == "exon"])


def _single_row(results):
    assert len(results) == 1
    return results.iloc[0]


@pytest.fixture(params=["+", "-"])
def strand(request):
    return request.param


def test_deletion_in_the_start_codon(tmp_path, strand):
    # ATG to AG: no placement keeps an ATG at the start
    row = _single_row(_run(tmp_path, strand, 15, "T", ""))
    assert pd.isna(row["unknown_reason"])
    assert row["alt_cds_seq"] == "AGGCTCT" + _REF_CDS[8:]
    assert row["start_loss"] == True
