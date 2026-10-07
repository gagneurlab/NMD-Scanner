"""
Tests of the reading frame, the start codon and start loss, through annotate() on one synthetic transcript and a GFF3
file.

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
        # The layout is the chromosome in transcript orientation: the flanks, the exons and the introns
        self.layout = self.flank + self.intron.join(exons) + self.flank
        self.genome = self.layout if strand == "+" else str(Seq(self.layout).reverse_complement())

        # Layout position and 0-based genomic position of each transcript position
        positions = []
        offset = len(self.flank)
        for exon in exons:
            positions.extend(range(offset, offset + len(exon)))
            offset += len(exon) + len(self.intron)
        self.layout_positions = positions
        self.positions = positions if strand == "+" else [len(self.layout) - 1 - p for p in positions]

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
        return self.layout_variant(self.layout_positions[position], ref, alt)

    def layout_variant(self, start, ref, alt):
        """
        VCF record (POS, REF, ALT) of a change in transcript orientation: ref at layout position ``start`` becomes
        alt. The change can reach into the introns. An indel gets a padding base on its left in the genome.
        """
        assert self.layout[start : start + len(ref)] == ref
        if self.strand == "-":
            start = len(self.layout) - start - len(ref)
            ref = str(Seq(ref).reverse_complement())
            alt = str(Seq(alt).reverse_complement())
        if len(ref) != len(alt):
            ref = self.genome[start - 1] + ref
            alt = self.genome[start - 1] + alt
            start -= 1
        return start + 1, ref, alt

    def run(self, position, ref, alt):
        """Run annotate() on this transcript and the variant (see ``variant``); return the single result row."""
        return self.run_record(*self.variant(position, ref, alt))

    def run_record(self, pos, vcf_ref, vcf_alt):
        """Run annotate() on this transcript and the VCF record; return the single result row."""
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


def test_a_stop_codon_as_annotated_start_codon_gives_no_ptc_distance(tmp_path, strand):
    """
    The annotated start codon is TAG (`*`), a stop codon. Translation cannot start on a stop codon, so this start
    codon is a misannotation. The missense GCC>GAC at t7 (`x`) leaves it unchanged. TAG is the first in-frame stop
    codon of the alt CDS, so the PTC is the start codon itself. ptc_to_start_codon is null, and so are
    ptc_less_than_150nt_to_start and the start-proximal rule.

    5' [uuu***=x=======sssuuuuu] 3'
    tx  0  3   7       15 18
    CDS    0           12
    """
    tx = SyntheticTranscript(tmp_path, strand, ["GGG" + "TAGGCCAAGCTG" + "TAA" + "GGGGG"], 3, 15)

    row = tx.run(7, "C", "A")

    _assert_values(
        row,
        {
            "alt_first_stop_pos": 0,
            "alt_has_ptc": True,
            "start_loss": False,
            "stop_loss": False,
            "ptc_to_start_codon": None,
            "ptc_less_than_150nt_to_start": None,
            "nmd_start_proximal_rule": None,
            "likely_misannotated": False,
        },
    )


@pytest.mark.parametrize(("start_codon", "start_pos", "start_exon"), [(False, None, None)], ids=["cds_start_nf"])
def test_stop_loss_scan_starts_at_the_annotated_start_codon(tmp_path, strand, start_codon, start_pos, start_exon):
    """
    The CDS is CTG AAA ATG CCC, then the stop codon TAA at t14 (`s`). CTG at t2 (`x`) is the annotated start codon.
    TAA>CAA at t14 loses the stop codon. The scan reads on from t2 to the TAG at t20 (`t`) in the 3' UTR. Its start
    codon is the annotated CTG at t2, not the in-frame ATG at t8 (`a`).

         7 nt          18 nt
    5' [uuxxx==]|[=a=====sssuuutttuu] 3'
    tx  0 2       7      14    20   25
    """
    cds = "CTGAAAATGCCC"
    exons = ["GG" + cds[:5], cds[5:] + "TAA" + "GGGTAGCC"]
    tx = SyntheticTranscript(tmp_path, strand, exons, 2, 14, start_codon=start_codon)

    row = tx.run(14, "T", "C")

    _assert_values(
        row,
        {
            "start_loss": False,
            "stop_loss": True,
            "alt_scan_start_codon_pos": start_pos,
            "alt_scan_start_codon_exon": start_exon,
            "alt_scan_first_stop_codon": "TAG",
            "alt_scan_first_stop_pos": 20,
            "annotated_stop_distance": -6,
        },
    )


def test_ptc_of_the_scan_after_a_start_codon_deletion_lies_in_the_exon_of_the_alt_transcript(tmp_path, strand):
    """
    The deletion of the T of the start codon ATG at t5 (`x`) is a start loss, and it shortens exon 1 to 10 nt in the
    alt transcript. The scan takes the ATG at t7 (`a`), out of frame, and its first stop codon is the TAA at t10 (`*`).
    In the alt transcript, t10 is the first base of exon 2, the last exon. So the PTC exon is exon 2, and ptc_to_exon_end
    runs to the transcript end. The ref exon lengths would put t10 into exon 1.

    ref 5' [uuuu=x=====]|[========uu] 3'
        tx 0    4        11         21
    alt 5' [uuuu===aaa]|[***=====uu] 3'
        tx 0    4  7    10         20
                             *-------->|  ptc_to_exon_end = 10
                   <---->  ptc_to_start_codon = 3
    """
    tx = SyntheticTranscript(tmp_path, strand, ["GACCATGGATG", "TAAGCTAAGC"], 4, 16)

    row = tx.run(5, "T", "")

    _assert_values(
        row,
        {
            "start_loss": True,
            "stop_loss": False,
            "alt_has_ptc": True,
            "alt_transcript_seq": "GACCAGGATGTAAGCTAAGC",
            "transcript_exons": [{"exon_number": 1, "length": 11}, {"exon_number": 2, "length": 10}],
            "alt_transcript_exons": [{"exon_number": 1, "length": 10}, {"exon_number": 2, "length": 10}],
            "alt_scan_start_codon_pos": 7,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_pos": 10,
            "alt_scan_stop_codons": [{"position": 10, "codon": "TAA"}],
            "alt_scan_stop_codon_exons": [2],
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 3,
            "ptc_exon_length": 10,
            "ptc_to_exon_end": 10,
            "annotated_stop_distance": 5,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_escape": True,
        },
    )


def test_scan_after_a_deletion_of_5utr_and_start_codon_bases_takes_the_alt_exon_lengths(tmp_path, strand):
    """
    The deletion CAT at t3 to t5 (`x`) takes the last 5' UTR base C and the AT of the start codon: a start loss that
    shortens the 5' UTR by 1 nt and the CDS by 2 nt. Exon 1 has 8 nt in the alt transcript. The scan takes the ATG at
    t5 (`a`), and its first stop codon is the TAA at t8 (`*`), the first base of exon 2 in the alt transcript. Without
    the 5' UTR change, exon 1 would have 9 nt and hold t8.

    ref 5' [uuuxxx=====]|[========uu] 3'
        tx 0   3         11         21
    alt 5' [uuu==aaa]|[***=====uu] 3'
        tx 0   3 5    8          18
    """
    tx = SyntheticTranscript(tmp_path, strand, ["GACCATGGATG", "TAAGCTAAGC"], 4, 16)

    row = tx.run(3, "CAT", "")

    _assert_values(
        row,
        {
            "start_loss": True,
            "alt_has_ptc": True,
            "alt_transcript_seq": "GACGGATGTAAGCTAAGC",
            "alt_cds_start_in_transcript": 3,
            "transcript_exons": [{"exon_number": 1, "length": 11}, {"exon_number": 2, "length": 10}],
            "alt_transcript_exons": [{"exon_number": 1, "length": 8}, {"exon_number": 2, "length": 10}],
            "alt_scan_start_codon_pos": 5,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_pos": 8,
            "alt_scan_stop_codon_exons": [2],
            "total_exon_count": 2,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 3,
            "ptc_exon_length": 10,
            "ptc_to_exon_end": 10,
            "nmd_last_exon_rule": True,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_single_exon_rule": False,
        },
    )


def test_scan_after_a_deletion_of_the_stop_codon_and_3utr_bases_takes_the_alt_exon_lengths(tmp_path, strand):
    """
    The deletion TAACCA at t10 to t15 (`x`) takes the stop codon TAA and 3 nt of the 3' UTR of exon 2: a stop loss
    that shortens exon 2 to 1 nt in the alt transcript. The scan reads on in the frame of the CDS to the TGA at t13
    (`*`), which lies in exon 3 of the alt transcript. Without the 3' UTR change, exon 2 would have 4 nt and hold t13.

    ref 5' [uuuu======]|[xxxxxxu]|[uuuuuuu] 3'
        tx 0    4       10        17       24
    alt 5' [uuuu======]|[u]|[uu***uu] 3'
        tx 0    4       10  11 13   18
    """
    tx = SyntheticTranscript(tmp_path, strand, ["GACCATGAAG", "TAACCAC", "CATGACC"], 4, 10)

    row = tx.run(10, "TAACCA", "")

    _assert_values(
        row,
        {
            "stop_loss": True,
            "alt_has_ptc": False,
            "alt_transcript_seq": "GACCATGAAGCCATGACC",
            "transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 7},
                {"exon_number": 3, "length": 7},
            ],
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 1},
                {"exon_number": 3, "length": 7},
            ],
            "alt_scan_start_codon_pos": 4,
            "alt_scan_start_codon_exon": 1,
            "alt_scan_first_stop_pos": 13,
            "alt_scan_stop_codons": [{"position": 13, "codon": "TGA"}],
            "alt_scan_stop_codon_exons": [3],
            "annotated_stop_distance": -9,
        },
    )


def test_ptc_features_take_the_3utr_length_change_in_the_ptc_exon(tmp_path, strand):
    """
    The delins CCCGGGTTTGCCTAACCA>TGAGGGTTTGCCTAACCTTC at t160 to t177 (`x`) reaches from the CDS over the stop codon
    into the 3' UTR of exon 2. Its bases are matched from the left, so its length change of +2 nt lies at the UTR end:
    the CDS keeps its length and gets the PTC TGA at t160 (`*`), and the 3' UTR of exon 2 becomes 2 nt longer
    (CCA>CCTTC). In the alt transcript, exon 2 has 55 nt, and the PTC lies 52 nt upstream of its end, the last exon
    junction: no 50 nt rule. With the ref exon lengths, the junction would lie 50 nt downstream of the PTC.

    ref 5' [uuuu===...===]|[===xxxxxxxxxxxxxxxxxxuuuuuuuuuu]|[uuuuuuuuuu] 3'
        tx 0    4         157 160            172  178         210        220
    alt 5' [uuuu===...===]|[===*=============uuuuuuuuuuuuuuu]|[uuuuuuuuuu] 3'
        tx 0    4         157 160            172               212        222
                              *------------------------------>|  ptc_to_exon_end = 52
    """
    exons = ["GACC" + "ATG" + "GCA" * 50, "GCACCCGGGTTTGCCTAACCA" + "CCCA" * 8, "CCACCACCAC"]
    tx = SyntheticTranscript(tmp_path, strand, exons, 4, 172)

    row = tx.run(160, "CCCGGGTTTGCCTAACCA", "TGAGGGTTTGCCTAACCTTC")

    _assert_values(
        row,
        {
            "start_loss": False,
            "stop_loss": False,
            "alt_has_ptc": True,
            "alt_first_stop_pos": 156,
            "alt_cds_start_in_transcript": 4,
            "transcript_exons": [
                {"exon_number": 1, "length": 157},
                {"exon_number": 2, "length": 53},
                {"exon_number": 3, "length": 10},
            ],
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 157},
                {"exon_number": 2, "length": 55},
                {"exon_number": 3, "length": 10},
            ],
            "total_exon_count": 3,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 156,
            "ptc_exon_length": 55,
            "ptc_to_exon_end": 52,
            "annotated_stop_distance": 12,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": False,
        },
    )


def test_a_deleted_exon_is_no_upstream_exon_of_the_ptc(tmp_path, strand):
    """
    The deletion (`x`) takes exon 2 with the intron bases next to it, from the sixth base of intron 1 to the fourth
    base of intron 2. It keeps the splice sites: AG stays before the deletion and GT after it. So the alt transcript
    skips exon 2, which has length 0 there. The CDS loses its 7 nt in exon 2, and the frameshift gives the PTC TAG at
    t13 (`*`) in exon 3. The deleted exon is not in the mRNA, so the PTC has 1 upstream exon, exon 1. The exon numbers
    of the ref transcript would give 2.

    ref 5' [uuuu======]|[xxxxxxx]|[============]|[==uuuuuuuuu] 3'
        tx 0    4      10        17             29            40
    alt 5' [uuuu======]|[===*========]|[==uuuuuuuuu] 3'
        tx 0    4      10  13         22            33
                           *-------->|  ptc_to_exon_end = 9
    """
    tx = SyntheticTranscript(tmp_path, strand, ["GACCATGGCA", "GCAGCAG", "GCCTAGCCGCCG", "CATAAGCCACC"], 4, 31)
    # From the sixth base of intron 1, after the last base of exon 1 at t9, to the fourth base of intron 2, before
    # the first base of exon 3 at t17
    start = tx.layout_positions[9] + 1 + 5
    end = tx.layout_positions[17] - len(tx.intron) + 4

    row = tx.run_record(*tx.layout_variant(start, tx.layout[start:end], ""))

    _assert_values(
        row,
        {
            "start_loss": False,
            "stop_loss": False,
            "alt_has_ptc": True,
            "alt_cds_seq": "ATGGCAGCCTAGCCGCCGCATAA",
            "alt_cds_exons": [
                {"exon_number": 1, "length": 6},
                {"exon_number": 2, "length": 0},
                {"exon_number": 3, "length": 12},
                {"exon_number": 4, "length": 5},
            ],
            "alt_first_stop_pos": 9,
            "alt_stop_codon_exons": [3],
            "transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 7},
                {"exon_number": 3, "length": 12},
                {"exon_number": 4, "length": 11},
            ],
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 10},
                {"exon_number": 2, "length": 0},
                {"exon_number": 3, "length": 12},
                {"exon_number": 4, "length": 11},
            ],
            "total_exon_count": 4,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 9,
            "ptc_exon_length": 12,
            "ptc_to_exon_end": 9,
            "annotated_stop_distance": 11,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": True,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
        },
    )


def test_long_exon_rule_takes_the_length_of_the_ptc_exon_in_the_alt_transcript(tmp_path, strand):
    """
    The deletion of the G at t157 (`x`), in exon 2 of 408 nt, gives the PTC TAA at t157 (`*`), 153 nt downstream of
    the start codon. In the alt transcript, exon 2 has 407 nt: no long exon rule, and the PTC escapes NMD by no rule.
    With the length of exon 2 in the ref transcript, the long exon rule would fire.

    ref 5' [uuuu===...===]|[===x=====...=====]|[======uuuuu] 3'
        tx 0    4         154 157               562         573
    alt 5' [uuuu===...===]|[===*====...=====]|[======uuuuu] 3'
        tx 0    4         154 157              561         572
                              <----------------->  ptc_exon_length = 407
                              *---------------->|  ptc_to_exon_end = 404
    """
    exons = ["GACC" + "ATG" + "GCA" * 49, "GCAGTAAGC" + "GCA" * 133, "GCATAACCACC"]
    tx = SyntheticTranscript(tmp_path, strand, exons, 4, 565)

    row = tx.run(157, "G", "")

    _assert_values(
        row,
        {
            "start_loss": False,
            "stop_loss": False,
            "alt_has_ptc": True,
            "alt_first_stop_pos": 153,
            "transcript_exons": [
                {"exon_number": 1, "length": 154},
                {"exon_number": 2, "length": 408},
                {"exon_number": 3, "length": 11},
            ],
            "alt_transcript_exons": [
                {"exon_number": 1, "length": 154},
                {"exon_number": 2, "length": 407},
                {"exon_number": 3, "length": 11},
            ],
            "total_exon_count": 3,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 153,
            "ptc_exon_length": 407,
            "ptc_to_exon_end": 404,
            "annotated_stop_distance": 407,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": False,
        },
    )
