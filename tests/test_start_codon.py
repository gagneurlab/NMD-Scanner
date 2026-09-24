"""
Tests of the start codon and of start loss, through annotate() on one synthetic transcript and a GFF3 file.

The transcript drawings show a transcript 5' to 3' and are not to scale. `u` is UTR, `=` is CDS, `[...]` is an exon,
`|` between two exons is an exon junction, and `x` marks the changed bases. The numbers under a transcript are
transcript positions (tx). A drawing is in transcript orientation, also on the minus strand.
"""

import pandas as pd
import pytest
from Bio.Seq import Seq

from nmd_scanner.cli import annotate


class SyntheticTranscript:
    """
    A transcript on its own synthetic chromosome, for runs of annotate() on both strands.

    The transcript is given in transcript orientation (5' to 3'). Its exons are separated by introns, and the
    chromosome holds them on the given strand. The annotation is a GENCODE GFF3: its CDS rows include the stop codon,
    and start_codon and stop_codon rows mark the two codons. A transcript without an annotated start codon is tagged
    cds_start_NF.
    """

    flank = "CCCCCCCCCC"
    intron = "GTAAGTCCCCCCCCTTTCAG"

    def __init__(self, tmp_path, strand, exons, cds_start, stop_codon, frame=0, start_codon=True):
        """
        :param tmp_path: directory for the input files. Its name names the chromosome, because catch_sequence
            caches sequences by chromosome across tests.
        :param strand: "+" or "-"
        :param exons: exon sequences in transcript order
        :param cds_start: transcript position of the first CDS base
        :param stop_codon: transcript position of the first base of the stop codon, which ends the CDS
        :param frame: phase of the 5'-most CDS row: the number of bases before the first complete codon
        :param start_codon: whether the first 3 CDS bases are annotated as start codon
        """
        self.tmp_path = tmp_path
        self.strand = strand
        self.chrom = f"chr_{tmp_path.name}"
        layout = self.flank + self.intron.join(exons) + self.flank
        self.genome = layout if strand == "+" else str(Seq(layout).reverse_complement())

        # 0-based genomic position of each transcript position
        positions = []
        offset = len(self.flank)
        for exon in exons:
            positions.extend(range(offset, offset + len(exon)))
            offset += len(exon) + len(self.intron)
        self.positions = positions if strand == "+" else [len(layout) - 1 - p for p in positions]

        tag = "" if start_codon else ";tag=cds_start_NF"
        types = "gene_type=protein_coding;transcript_type=protein_coding"
        self.lines = [
            self._line("gene", 0, len(positions), ".", "ID=g1;gene_id=g1;gene_type=protein_coding"),
            self._line(
                "transcript", 0, len(positions), ".", f"ID=tx1;Parent=g1;gene_id=g1;transcript_id=tx1;{types}{tag}"
            ),
        ]
        features = [("CDS", cds_start, stop_codon + 3), ("stop_codon", stop_codon, stop_codon + 3)]
        if start_codon:
            features.append(("start_codon", cds_start, cds_start + 3))
        exon_start = 0
        phase = frame
        for number, exon in enumerate(exons, start=1):
            exon_end = exon_start + len(exon)
            parts = [("exon", exon_start, exon_end)]
            parts += [(feature, max(exon_start, start), min(exon_end, end)) for feature, start, end in features]
            for feature, start, end in parts:
                if start >= end:
                    continue
                attributes = (
                    f"ID={feature}:tx1:{number};Parent=tx1;gene_id=g1;transcript_id=tx1;{types};exon_number={number}"
                )
                self.lines.append(self._line(feature, start, end, phase if feature == "CDS" else ".", attributes + tag))
                if feature == "CDS":
                    # phase of the next CDS row: the bases that its first codon still needs
                    phase = (3 - (end - start - phase) % 3) % 3
            exon_start = exon_end

    def _line(self, feature, start, end, phase, attributes):
        """GFF3 line of a feature at transcript positions start to end (0-based, half-open)."""
        genomic = self.positions[start:end]
        columns = [self.chrom, "test", feature, min(genomic) + 1, max(genomic) + 1, ".", self.strand, phase, attributes]
        return "\t".join(str(column) for column in columns)

    def variant(self, position, ref, alt):
        """
        VCF record (POS, REF, ALT) of a change in transcript orientation within one exon: ref at transcript position
        ``position`` becomes alt. An indel gets a padding base on its left in the genome.
        """
        if self.strand == "+":
            start = self.positions[position]
        else:
            start = self.positions[position + len(ref) - 1]
            ref, alt = str(Seq(ref).reverse_complement()), str(Seq(alt).reverse_complement())
        assert self.genome[start : start + len(ref)] == ref
        if len(ref) != len(alt):
            start, ref, alt = start - 1, self.genome[start - 1] + ref, self.genome[start - 1] + alt
        return start + 1, ref, alt

    def run(self, position, ref, alt):
        """Run annotate() on this transcript and the variant (see ``variant``); return the single result row."""
        pos, vcf_ref, vcf_alt = self.variant(position, ref, alt)
        (self.tmp_path / "genome.fa").write_text(f">{self.chrom}\n{self.genome}\n")
        (self.tmp_path / "tx.gff3").write_text("##gff-version 3\n" + "\n".join(self.lines) + "\n")
        (self.tmp_path / "variant.vcf").write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            f"{self.chrom}\t{pos}\tvar1\t{vcf_ref}\t{vcf_alt}\t.\t.\t.\n"
        )
        results = annotate(
            str(self.tmp_path / "variant.vcf"), str(self.tmp_path / "tx.gff3"), str(self.tmp_path / "genome.fa")
        )
        assert len(results) == 1
        return results.iloc[0]


def _assert_values(row, expected):
    """Assert the expected value of each column; None means null."""
    actual = {
        column: None if pd.api.types.is_scalar(row[column]) and pd.isna(row[column]) else row[column]
        for column in expected
    }
    assert actual == expected


@pytest.fixture(params=["+", "-"], ids=["plus", "minus"])
def strand(request):
    return request.param


def test_stop_loss_in_a_cds_without_a_leading_atg_reads_through_the_3utr(tmp_path, strand):
    """
    A cds_start_NF CDS without a leading ATG has no start codon to lose. No codon of its frame is ATG, in the ref and in
    the alt CDS. The variant TAA>CAA at `x` loses the stop codon, and the scan reads on in the frame of the CDS. Its
    first stop codon `s` is the TAG at t21 in the 3' UTR.

          8 nt         20 nt
    5' [uuu=====]|[=======xxxuuusssuuuu] 3'
    tx  0  3       8      15 18 21     28
    """
    cds = "CTGAAACCCGAC"
    tx = SyntheticTranscript(
        tmp_path, strand, ["GGG" + cds[:5], cds[5:] + "TAA" + "GGGTAGCCGG"], 3, 15, start_codon=False
    )

    row = tx.run(15, "T", "C")

    _assert_values(
        row,
        {
            "ref_start_codon_pos": None,
            "alt_start_codon_pos": None,
            "start_loss": False,
            "stop_loss": True,
            "alt_is_premature": False,
            "transcript_start_codon_pos": None,
            "transcript_first_stop_codon": "TAG",
            "transcript_first_stop_pos": 21,
            "transcript_all_stop_codons": [(21, "TAG")],
            "stop_codon_distance": -6,
        },
    )


def test_lost_internal_atg_is_no_start_loss(tmp_path, strand):
    """
    A cds_start_NF CDS starts with CTG, and its first in-frame ATG lies at t9, CDS position 6. The variant ATG>ACG at
    `x`, t10, removes that ATG. It is an internal Met, so the CDS has no start codon to lose. `s` is the stop codon.

          8 nt         20 nt
    5' [uuu=====]|[==x=======sssuuuuuuu] 3'
    tx  0  3       8 10      18 21     28
    """
    cds = "CTGAAAATGCCCGAC"
    tx = SyntheticTranscript(tmp_path, strand, ["GGG" + cds[:5], cds[5:] + "TAA" + "GGGGGGG"], 3, 18, start_codon=False)

    row = tx.run(10, "T", "C")

    _assert_values(
        row,
        {
            "ref_start_codon_pos": 6,
            "alt_start_codon_pos": None,
            "start_loss": False,
            "stop_loss": False,
            "alt_is_premature": False,
            "transcript_start_codon_pos": None,
            "transcript_first_stop_pos": None,
        },
    )
