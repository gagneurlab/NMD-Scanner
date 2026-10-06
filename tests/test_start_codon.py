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
            ref, alt = str(Seq(ref).reverse_complement()), str(Seq(alt).reverse_complement())
        if len(ref) != len(alt):
            start, ref, alt = start - 1, self.genome[start - 1] + ref, self.genome[start - 1] + alt
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
    `x`, t10, removes that ATG. It is an internal Met, so the CDS has no start codon to lose, and neither CDS has a
    start codon position. `s` is the stop codon.

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
            "ref_start_codon_pos": None,
            "alt_start_codon_pos": None,
            "start_loss": False,
            "stop_loss": False,
            "alt_is_premature": False,
            "transcript_start_codon_pos": None,
            "transcript_first_stop_pos": None,
        },
    )


@pytest.mark.parametrize("frame", [0, 1, 2])
@pytest.mark.parametrize(
    ("variant", "expected"),
    [
        (
            (5, "C", "T"),
            {"alt_is_premature": False, "alt_first_stop_pos": 15, "alt_all_stop_codons": [(15, "TAA")], "distance": 0},
        ),
        (
            (12, "C", "T"),
            {
                "alt_is_premature": True,
                "alt_first_stop_pos": 9,
                "alt_all_stop_codons": [(9, "TAA"), (15, "TAA")],
                "distance": 6,
            },
        ),
    ],
    ids=["out_of_frame_stop", "in_frame_stop"],
)
def test_codon_scan_starts_at_the_first_complete_codon(tmp_path, strand, frame, variant, expected):
    """
    The CDS of a cds_start_NF transcript starts with `frame` bases (`f`) that belong to no complete codon: the GFF3
    phase of its CDS row. The codons TGC AAA CCC CAA GGC and the stop codon TAA follow. The drawing shows frame 1, with
    the first complete codon at t4. TGC>TGT at t6 (`a`) puts a TAA out of frame. CAA>TAA at t13 (`b`) is an in-frame
    PTC. The positions in the test are those of frame 0, shifted by `frame`.

          8 nt           19 nt
    5' [uuuf==a=]|[=====b=====sssuuuuu] 3'
    tx  0  3       8    13    19      27
    """
    cds = "AC"[:frame] + "TGCAAACCCCAAGGC"
    exons = ["GGG" + cds[:5], cds[5:] + "TAA" + "GGGGG"]
    tx = SyntheticTranscript(tmp_path, strand, exons, 3, 3 + len(cds), frame=frame, start_codon=False)
    position, ref, alt = variant

    row = tx.run(position + frame, ref, alt)

    _assert_values(
        row,
        {
            "cds_frame": frame,
            "ref_start_codon_pos": None,
            "ref_all_stop_codons": [(15 + frame, "TAA")],
            "alt_first_stop_pos": expected["alt_first_stop_pos"] + frame,
            "alt_all_stop_codons": [(pos + frame, codon) for pos, codon in expected["alt_all_stop_codons"]],
            "alt_is_premature": expected["alt_is_premature"],
            "start_loss": False,
            "stop_loss": False,
            "stop_codon_distance": expected["distance"],
        },
    )


def test_stop_loss_reads_through_the_3utr_in_the_cds_frame(tmp_path, strand):
    """
    The CDS of a cds_start_NF transcript has phase 1: its first base A belongs to no complete codon (`f`). The codons
    TGC AAA CCC and the stop codon TAA (`x`) follow, from t4. TAA>CAA at t13 loses the stop codon. The scan reads on in
    the frame of t4, to the TGA at t25 (`s`) in the 3' UTR. Read from t3, the alt transcript has no stop codon.

          8 nt              22 nt
    5' [uuuf====]|[=====xxxuuuuuuuuusssuu] 3'
    tx  0  3       8    13 16       25   30
    """
    cds = "ATGCAAACCC"
    exons = ["GGG" + cds[:5], cds[5:] + "TAA" + "GGGCCCGGGTGACC"]
    tx = SyntheticTranscript(tmp_path, strand, exons, 3, 13, frame=1, start_codon=False)

    row = tx.run(13, "T", "C")

    _assert_values(
        row,
        {
            "cds_frame": 1,
            "ref_start_codon_pos": None,
            "alt_first_stop_pos": None,
            "alt_is_premature": False,
            "start_loss": False,
            "stop_loss": True,
            "transcript_start_codon_pos": None,
            "transcript_first_stop_codon": "TGA",
            "transcript_first_stop_pos": 25,
            "transcript_all_stop_codons": [(25, "TGA")],
            "stop_codon_distance": -12,
        },
    )


@pytest.mark.parametrize(
    ("cds", "start_codon", "variant", "expected"),
    [
        (
            "CTGAAACCCGAC",
            True,
            (4, "T", "C"),
            {"ref_start_codon_pos": 0, "alt_start_codon_pos": None, "start_loss": True, "scanned_stops": 0},
        ),
        (
            "CTGAAACCCGAC",
            True,
            (10, "C", "G"),
            {"ref_start_codon_pos": 0, "alt_start_codon_pos": 0, "start_loss": False, "scanned_stops": None},
        ),
        (
            "ATGAAACCCGAC",
            True,
            (4, "T", "C"),
            {"ref_start_codon_pos": 0, "alt_start_codon_pos": None, "start_loss": True, "scanned_stops": 0},
        ),
        (
            "ATGAAACCCGAC",
            False,
            (4, "T", "C"),
            {"ref_start_codon_pos": None, "alt_start_codon_pos": None, "start_loss": False, "scanned_stops": None},
        ),
    ],
    ids=["ctg_to_ccg", "missense_after_ctg", "atg_to_acg", "atg_to_acg_without_start_codon_rows"],
)
def test_start_loss_is_a_change_of_the_annotated_start_codon(tmp_path, strand, cds, start_codon, variant, expected):
    """
    The CDS starts with the start codon CTG or ATG at t3 (`x`), which the GFF3 marks with a start_codon row or not.
    CTG>CCG and ATG>ACG change its second base, at t4. The missense CCC>CGC changes t10 in exon 2. `s` is the stop
    codon. After a start loss, the scan looks for an ATG from t3 on and finds none, so it reads no stop codon.

          8 nt         15 nt
    5' [uuuxxx==]|[=======sssuuuuu] 3'
    tx  0  3       8      15 18   23
    """
    tx = SyntheticTranscript(
        tmp_path, strand, ["GGG" + cds[:5], cds[5:] + "TAA" + "GGGGG"], 3, 15, start_codon=start_codon
    )

    row = tx.run(*variant)

    _assert_values(
        row,
        {
            "has_start_codon": start_codon,
            "ref_start_codon_pos": expected["ref_start_codon_pos"],
            "alt_start_codon_pos": expected["alt_start_codon_pos"],
            "start_loss": expected["start_loss"],
            "stop_loss": False,
            "transcript_start_codon_pos": None,
            "transcript_num_stop_codons": expected["scanned_stops"],
        },
    )


@pytest.mark.parametrize(
    ("start_codon", "expected"),
    [
        (
            True,
            {
                "ref_start_codon_pos": 0,
                "alt_start_codon_pos": 0,
                "ptc_to_start_codon": 159,
                "likely_misannotated": False,
            },
        ),
        (
            False,
            {
                "ref_start_codon_pos": None,
                "alt_start_codon_pos": None,
                "ptc_to_start_codon": None,
                "likely_misannotated": True,
            },
        ),
    ],
    ids=["start_codon", "cds_start_nf"],
)
def test_ptc_distance_is_measured_from_the_annotated_start_codon(tmp_path, strand, start_codon, expected):
    """
    The CDS is CTG, 9 times AAA, ATG, 42 times AAA, TGG, GAC, then the stop codon TAA. CTG at t3 (`x`) is the annotated
    start codon, or the transcript has none, as one tagged cds_start_NF. The in-frame ATG at CDS position 30 (`a`) is
    an internal Met. TGG>TAG at t163 is a PTC (`*`) at CDS position 159: 159 nt from the start codon CTG, but only 129
    nt from the ATG, < 150. Without an annotated start codon, the start lies upstream of the CDS, at an unknown
    distance.

              103 nt                73 nt
    5' [uuuxxx===a=========]|[=====*===sssuuuuu] 3'
    tx  0  3                  103  162 168     176
    CDS    0     30                159 165
           <----------------------->  ptc_to_start_codon = 159, not < 150
    """
    cds = "CTG" + "AAA" * 9 + "ATG" + "AAA" * 42 + "TGG" + "GAC"
    exons = ["GGG" + cds[:100], cds[100:] + "TAA" + "GGGGG"]
    tx = SyntheticTranscript(tmp_path, strand, exons, 3, 3 + len(cds), start_codon=start_codon)

    row = tx.run(163, "G", "A")

    _assert_values(
        row,
        {
            **expected,
            "alt_is_premature": True,
            "alt_first_stop_pos": 159,
            "start_loss": False,
            "ptc_less_than_150nt_to_start": False,
            "nmd_start_proximal_rule": False,
        },
    )


@pytest.mark.parametrize(
    ("start_codon", "start_pos", "start_exon"), [(True, 2, 1), (False, None, None)], ids=["start_codon", "cds_start_nf"]
)
def test_stop_loss_scan_starts_at_the_annotated_start_codon(tmp_path, strand, start_codon, start_pos, start_exon):
    """
    The CDS is CTG AAA ATG CCC, then the stop codon TAA at t14 (`s`). CTG at t2 (`x`) is the annotated start codon,
    or the transcript has none, as one tagged cds_start_NF. TAA>CAA at t14 loses the stop codon. The scan reads on from
    t2 to the TAG at t20 (`t`) in the 3' UTR. Its start codon is the annotated CTG at t2, not the in-frame ATG at t8
    (`a`). Without an annotated start codon, it has none.

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
            "transcript_start_codon_pos": start_pos,
            "transcript_start_codon_exon": start_exon,
            "transcript_first_stop_codon": "TAG",
            "transcript_first_stop_pos": 20,
            "stop_codon_distance": -6,
        },
    )


def test_start_proximal_rule_after_a_start_loss_measures_from_the_atg_of_the_scan(tmp_path, strand):
    """
    The CDS is ATG GCC ATG CAG, 36 times AAA, GAC, then the stop codon TAA at t126 (`s`). One MNV changes the start
    codon ATG>ACG at t3 (`x`) and CAG>TAG at t12 (`*`). The start loss makes the scan look for the next ATG from t3 on.
    It finds the in-frame ATG at t9 (`a`), whose first stop codon is the TAG at t12. This TAG is the PTC: 3 nt
    downstream of the ATG of the scan, < 150. Measured from the CDS start, it would be 9 nt.

              63 nt               60 nt          11 nt
    5' [uuuxxx===aaa***=======]|[=============]|[====sssuuuuu] 3'
    tx  0  3     9  12          63              123  126   134
                 <-->  ptc_to_start_codon = 3
    """
    cds = "ATGGCCATGCAG" + "AAA" * 36 + "GAC"
    exons = ["GGG" + cds[:60], cds[60:120], cds[120:] + "TAA" + "GGGGG"]
    tx = SyntheticTranscript(tmp_path, strand, exons, 3, 126)

    row = tx.run(4, "TGGCCATGC", "CGGCCATGT")

    _assert_values(
        row,
        {
            "start_loss": True,
            "alt_start_codon_pos": None,
            "alt_is_premature": True,
            "transcript_start_codon_pos": 9,
            "transcript_first_stop_pos": 12,
            "ptc_to_start_codon": 3,
            "ptc_less_than_150nt_to_start": True,
            "nmd_start_proximal_rule": True,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
        },
    )


def test_start_loss_classifies_the_first_stop_codon_of_the_rescued_orf(tmp_path, strand):
    """
    The CDS is ATG GAT GCC CTA, 50 times AAA, GAC, then the stop codon TAA at t168 (`s`). ATG>ACG at t4 (`x`) is a
    start loss. Translation starts at the next ATG, at t7 (`a`), which is out of frame with the annotated stop codon.
    Its ORF ATG CCC TAA ends at the TAA at t13 (`*`). Upstream of the annotated stop codon, this stop codon is a PTC,
    as after a frameshift. The PTC features measure it in alt CDS coordinates: CDS position 10 in exon 1. The frame of
    the CDS start reads on to the annotated stop codon, at alt_first_stop_pos 165.

              63 nt                60 nt           53 nt
    5' [uuuxxx=aaa===***=======]|[=============]|[====sssuuuuu] 3'
    tx  0  3   7     13           63              123  168   176
               <---->  ptc_to_start_codon = 6
                     <-------------------------------->  stop_codon_distance = 155
                     *--------->|  ptc_to_intron = 50
    """
    cds = "ATGGATGCCCTA" + "AAA" * 50 + "GAC"
    exons = ["GGG" + cds[:60], cds[60:120], cds[120:] + "TAA" + "GGGGG"]
    tx = SyntheticTranscript(tmp_path, strand, exons, 3, 168)

    row = tx.run(4, "T", "C")

    _assert_values(
        row,
        {
            "start_loss": True,
            "stop_loss": False,
            "alt_is_premature": True,
            "alt_first_stop_pos": 165,
            "transcript_start_codon_pos": 7,
            "transcript_first_stop_pos": 13,
            "stop_codon_distance": 155,
            "upstream_exon_count": 0,
            "downstream_exon_count": 2,
            "ptc_exon_length": 63,
            "ptc_to_intron": 50,
            "ptc_to_start_codon": 6,
            "ptc_less_than_150nt_to_start": True,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
        },
    )


def test_start_loss_with_a_frameshift_reads_the_frame_of_the_next_atg(tmp_path, strand):
    """
    The CDS is ATG GTA AGC ATG GCC AAA GAC, then the stop codon TAA. The deletion of the T of the start codon at t4 is a
    start loss and a frameshift. The drawing shows the alt transcript. Read from the alt CDS start at t3, the frame
    meets the TAA at t6 (`*`), across the exon junction. But translation starts at the next ATG, at t11 (`a`), which is
    the Met at CDS position 9 of the reference, in frame with the annotated stop codon. Its ORF ends at the annotated
    stop codon, at t23 (`s`). So the row is neither a PTC nor a stop loss.

           7 nt              24 nt
    5' [uuu===*]|[**==aaa======sssuuuuu] 3'
    tx  0  3  6   7   11       23     31
    """
    cds = "ATGGTAAGCATGGCCAAAGAC"
    tx = SyntheticTranscript(tmp_path, strand, ["GGG" + cds[:5], cds[5:] + "TAA" + "GGGGG"], 3, 24)

    row = tx.run(4, "T", "")

    _assert_values(
        row,
        {
            "start_loss": True,
            "stop_loss": False,
            "alt_is_premature": False,
            "alt_first_stop_pos": 3,
            "transcript_start_codon_pos": 11,
            "transcript_first_stop_pos": 23,
            "stop_codon_distance": 0,
            "ptc_to_start_codon": None,
            "nmd_start_proximal_rule": False,
            "nmd_escape": False,
        },
    )


@pytest.mark.parametrize(
    ("utr3", "start_pos", "stop_pos"),
    [("GGATGCCCTGAGG", 20, 26), ("GGACGCCCTGAGG", None, None)],
    ids=["atg_in_the_3utr", "no_atg"],
)
def test_start_loss_without_an_atg_upstream_of_the_stop_codon_is_neither(tmp_path, strand, utr3, start_pos, stop_pos):
    """
    The CDS is ATG CCC AAA GAC, then the stop codon TAA at t15 (`s`). ATG>ACG at t4 (`x`) is a start loss, and no ATG
    follows in the CDS. In the first case, the next ATG lies in the 3' UTR, at t20 (`a`), and its ORF ends at the TGA
    at t26 (`t`). In the second case, the 3' UTR has ACG there, and the scan finds no ATG. Either way, no ORF overlaps
    the CDS: the row is neither a PTC nor a stop loss, and stop_codon_distance is null.

          8 nt              23 nt
    5' [uuuxxx==]|[=======sssuuaaauuutttuu] 3'
    tx  0  3       8      15   20    26   31
    """
    cds = "ATGCCCAAAGAC"
    tx = SyntheticTranscript(tmp_path, strand, ["GGG" + cds[:5], cds[5:] + "TAA" + utr3], 3, 15)

    row = tx.run(4, "T", "C")

    _assert_values(
        row,
        {
            "start_loss": True,
            "stop_loss": False,
            "alt_is_premature": False,
            "transcript_start_codon_pos": start_pos,
            "transcript_first_stop_pos": stop_pos,
            "stop_codon_distance": None,
            "ptc_to_start_codon": None,
            "nmd_escape": False,
        },
    )


def test_ptc_of_the_scan_after_a_start_codon_deletion_lies_in_the_exon_of_the_alt_transcript(tmp_path, strand):
    """
    The deletion of the T of the start codon ATG at t5 (`x`) is a start loss, and it shortens exon 1 to 10 nt in the
    alt transcript. The scan takes the ATG at t7 (`a`), out of frame, and its first stop codon is the TAA at t10 (`*`).
    In the alt transcript, t10 is the first base of exon 2, the last exon. So the PTC exon is exon 2, and ptc_to_intron
    runs to the transcript end. The ref exon lengths would put t10 into exon 1.

    ref 5' [uuuu=x=====]|[========uu] 3'
        tx 0    4        11         21
    alt 5' [uuuu===aaa]|[***=====uu] 3'
        tx 0    4  7    10         20
                             *-------->|  ptc_to_intron = 10
                   <---->  ptc_to_start_codon = 3
    """
    tx = SyntheticTranscript(tmp_path, strand, ["GACCATGGATG", "TAAGCTAAGC"], 4, 16)

    row = tx.run(5, "T", "")

    _assert_values(
        row,
        {
            "start_loss": True,
            "stop_loss": False,
            "alt_is_premature": True,
            "alt_transcript_seq": "GACCAGGATGTAAGCTAAGC",
            "transcript_exon_info": [(1, 11), (2, 10)],
            "alt_transcript_exon_info": [(1, 10), (2, 10)],
            "transcript_start_codon_pos": 7,
            "transcript_start_codon_exon": 1,
            "transcript_first_stop_pos": 10,
            "transcript_all_stop_codons": [(10, "TAA")],
            "transcript_stop_codon_exons": [2],
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 3,
            "ptc_exon_length": 10,
            "ptc_to_intron": 10,
            "stop_codon_distance": 5,
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
            "alt_is_premature": True,
            "alt_transcript_seq": "GACGGATGTAAGCTAAGC",
            "alt_cds_start_in_transcript": 3,
            "transcript_exon_info": [(1, 11), (2, 10)],
            "alt_transcript_exon_info": [(1, 8), (2, 10)],
            "transcript_start_codon_pos": 5,
            "transcript_start_codon_exon": 1,
            "transcript_first_stop_pos": 8,
            "transcript_stop_codon_exons": [2],
            "total_exon_count": 2,
            "upstream_exon_count": 1,
            "downstream_exon_count": 0,
            "ptc_to_start_codon": 3,
            "ptc_exon_length": 10,
            "ptc_to_intron": 10,
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
            "alt_is_premature": False,
            "alt_transcript_seq": "GACCATGAAGCCATGACC",
            "transcript_exon_info": [(1, 10), (2, 7), (3, 7)],
            "alt_transcript_exon_info": [(1, 10), (2, 1), (3, 7)],
            "transcript_start_codon_pos": 4,
            "transcript_start_codon_exon": 1,
            "transcript_first_stop_pos": 13,
            "transcript_all_stop_codons": [(13, "TGA")],
            "transcript_stop_codon_exons": [3],
            "stop_codon_distance": -9,
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
                              *------------------------------>|  ptc_to_intron = 52
    """
    exons = ["GACC" + "ATG" + "GCA" * 50, "GCACCCGGGTTTGCCTAACCA" + "CCCA" * 8, "CCACCACCAC"]
    tx = SyntheticTranscript(tmp_path, strand, exons, 4, 172)

    row = tx.run(160, "CCCGGGTTTGCCTAACCA", "TGAGGGTTTGCCTAACCTTC")

    _assert_values(
        row,
        {
            "start_loss": False,
            "stop_loss": False,
            "alt_is_premature": True,
            "alt_first_stop_pos": 156,
            "alt_cds_start_in_transcript": 4,
            "transcript_exon_info": [(1, 157), (2, 53), (3, 10)],
            "alt_transcript_exon_info": [(1, 157), (2, 55), (3, 10)],
            "total_exon_count": 3,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 156,
            "ptc_exon_length": 53,
            "ptc_to_intron": 52,
            "stop_codon_distance": 12,
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
                           *-------->|  ptc_to_intron = 9
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
            "alt_is_premature": True,
            "alt_cds_seq": "ATGGCAGCCTAGCCGCCGCATAA",
            "alt_cds_info": [(1, 6), (2, 0), (3, 12), (4, 5)],
            "alt_first_stop_pos": 9,
            "alt_stop_codon_exons": [3],
            "transcript_exon_info": [(1, 10), (2, 7), (3, 12), (4, 11)],
            "alt_transcript_exon_info": [(1, 10), (2, 0), (3, 12), (4, 11)],
            "total_exon_count": 4,
            "upstream_exon_count": 1,
            "downstream_exon_count": 1,
            "ptc_to_start_codon": 9,
            "ptc_exon_length": 12,
            "ptc_to_intron": 9,
            "stop_codon_distance": 11,
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": True,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": True,
            "nmd_single_exon_rule": False,
            "nmd_escape": True,
        },
    )
