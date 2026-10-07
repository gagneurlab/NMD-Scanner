"""
Runner of the conformance cases.

A case pins one constellation: one transcript layout and one variant. A few cases have more variants or none. The
runner writes a GFF3 (GENCODE or Ensembl flavor), a FASTA and a VCF file for the case, runs annotate() on them and
compares the result with the expected value of every output column. A case on input that annotate() rejects expects
an error instead.
Each case runs on both strands. The layout is given 5' to 3', in transcript orientation, and on the minus strand the
chromosome holds its reverse complement. A column in genomic coordinates has one expected value per strand.

The expected values come in two parts. A Layout gives the columns that depend on the transcript only (REF_COLUMNS),
and a Case gives the other columns (CASE_COLUMNS). Together they give every output column.

Each case has a drawing. The drawing holds the layout block that render_case() draws from the case: the transcript,
the change, and the marks, spans and rulers of the case. Hand-written sentences before or after the block explain the
constellation. A meta-test in test_conformance.py checks that each drawing holds its block.
"""

import dataclasses
import re
from collections.abc import Callable
from typing import Any, NamedTuple

import pandas as pd
import pytest
from Bio.Seq import Seq

from nmd_scanner import annotate
from nmd_scanner.schema import KIND_DTYPES, OUTPUT_COLUMN_KINDS

CHROMOSOME = "chr1"
FLANK = "CCCCCCCCCC"
INTRON = "GTAAGTCCCCCCCCTTTCAG"
TYPES = "gene_type=protein_coding;transcript_type=protein_coding"
FLAVORS = ("gencode", "ensembl")

# The columns that depend on the transcript only. A Layout gives their expected values.
REF_COLUMNS = (
    "transcript_id",
    "gene_id",
    "chrom",
    "strand",
    "cds_start",
    "cds_end",
    "ref_cds_seq",
    "ref_cds_length",
    "has_start_codon",
    "has_stop_codon",
    "cds_frame",
    "ref_cds_exons",
    "cds_in_transcript",
    "start_codon_exon",
    "ref_last_codon",
    "ref_valid_stop",
    "ref_first_stop_codon",
    "ref_first_stop_pos",
    "ref_stop_codon_count",
    "ref_stop_codons",
    "ref_stop_codon_exons",
    "ref_has_ptc",
    "transcript_start",
    "transcript_end",
    "transcript_seq",
    "transcript_length",
    "cds_start_in_transcript",
    "cds_end_in_transcript",
    "transcript_exons",
    "utr3_length",
    "utr5_length",
    "total_exon_count",
    "likely_misannotated",
)
# The columns that depend on the variant. A Case gives their expected values.
CASE_COLUMNS = tuple(column for column in OUTPUT_COLUMN_KINDS if column not in REF_COLUMNS)
# The columns that echo the VCF record
RECORD_COLUMNS = ("ref", "alt", "start", "end")


class PerStrand(NamedTuple):
    """The expected value of a column on the plus and on the minus strand."""

    plus: Any
    minus: Any


def per_strand(plus, minus):
    return PerStrand(plus, minus)


# Expected values that many cases share
IDS = {"transcript_id": "tx1", "gene_id": "g1", "chrom": CHROMOSOME, "strand": per_strand("+", "-")}
NOT_SCANNED = {
    "alt_scan_start_codon_pos": None,
    "alt_scan_start_codon_exon": None,
    "alt_scan_first_stop_codon": None,
    "alt_scan_first_stop_pos": None,
    "alt_scan_stop_codon_count": None,
    "alt_scan_stop_codons": None,
    "alt_scan_stop_codon_exons": None,
}
NO_PTC_FEATURES = {
    "upstream_exon_count": None,
    "downstream_exon_count": None,
    "ptc_to_start_codon": None,
    "ptc_less_than_150nt_to_start": None,
    "ptc_exon_length": None,
    "ptc_to_exon_end": None,
}
# The NMD rules of a row that is not a PTC row
NO_RULE = {
    "nmd_last_exon_rule": None,
    "nmd_50nt_penultimate_rule": None,
    "nmd_long_exon_rule": None,
    "nmd_start_proximal_rule": None,
    "nmd_single_exon_rule": None,
    "nmd_escape": None,
}
# The NMD rules of a PTC row whose PTC meets no rule
RULES_FALSE = dict.fromkeys(NO_RULE, False)
# A row with unknown_reason: the alt transcript is unknown, so every column of the alt side is null, and the
# model cannot score the row
UNKNOWN_ALT = {
    "alt_cds_seq": None,
    "alt_cds_length": None,
    "alt_cds_exons": None,
    "alt_last_codon": None,
    "alt_valid_stop": None,
    "alt_first_stop_codon": None,
    "alt_first_stop_pos": None,
    "alt_stop_codon_count": None,
    "alt_stop_codons": None,
    "alt_stop_codon_exons": None,
    "alt_has_ptc": None,
    "start_loss": None,
    "stop_loss": None,
    "alt_transcript_seq": None,
    "alt_transcript_length": None,
    "alt_cds_start_in_transcript": None,
    "alt_transcript_exons": None,
    **NOT_SCANNED,
    "upstream_exon_count": None,
    "downstream_exon_count": None,
    "ptc_to_start_codon": None,
    "ptc_less_than_150nt_to_start": None,
    "ptc_exon_length": None,
    "annotated_stop_distance": None,
    "ptc_to_exon_end": None,
    "nmd_last_exon_rule": None,
    "nmd_50nt_penultimate_rule": None,
    "nmd_long_exon_rule": None,
    "nmd_start_proximal_rule": None,
    "nmd_single_exon_rule": None,
    "nmd_escape": None,
    "nmd_model_status": "unknown_effect",
}


class _SameExons:
    """The type of SAME_EXONS."""

    def __repr__(self):
        return "SAME_EXONS"


# The expected alt_transcript_exons of a variant that keeps the length of every exon: the transcript_exons of
# the row
SAME_EXONS = _SameExons()


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


@dataclasses.dataclass(frozen=True)
class Transcript:
    """
    A transcript on its own chromosome, given 5' to 3' in transcript orientation.

    :param exons: exon sequences in transcript order. Upper case bases form the coding region: the CDS with its stop
        codon, as a GFF3 CDS row holds it. Lower case bases are UTR.
    :param introns: intron sequences, INTRON by default
    :param start_codon: whether the GFF3 has start_codon rows, on the first 3 coding bases
    :param stop_codon: whether the GFF3 has stop_codon rows, on the last 3 coding bases
    :param frame: GFF3 phase of the 5'-most CDS row: the number of bases before the first complete codon
    :param tags: tags of the transcript, e.g. cds_start_NF. GENCODE repeats them on every row of the transcript;
        Ensembl has them on the transcript row only.
    :param exon_rows: whether the GFF3 has the exon rows
    :param flanks: the chromosome bases before the 5' end and after the 3' end of the transcript
    :param contig: the name of the chromosome, e.g. chrM for the mitochondrial stop codons
    :param flavor: "gencode" or "ensembl". An Ensembl GFF3 has biotype attributes, an mRNA row, the ID/Parent
        hierarchy with gene: and transcript: prefixes, exon numbers as rank, and never start_codon or stop_codon rows.
    :param edit_gff3: a function that takes the GFF3 rows (lists of the 9 columns, as strings) and the strand, and
        returns the rows to write. For annotations that the other fields cannot express, e.g. a CDS outside the exons.
    :param fasta_contig: the name of the chromosome in the FASTA, if it differs from contig, e.g. for a FASTA that
        lacks the chromosome of the transcript
    """

    exons: tuple[str, ...]
    introns: tuple[str, ...] | None = None
    start_codon: bool = True
    stop_codon: bool = True
    frame: int = 0
    tags: tuple[str, ...] = ()
    exon_rows: bool = True
    flanks: tuple[str, str] = (FLANK, FLANK)
    contig: str = CHROMOSOME
    flavor: str = "gencode"
    edit_gff3: Callable[[list[list[str]], str], list[list[str]]] | None = None
    fasta_contig: str | None = None

    def __post_init__(self):
        assert self.flavor in FLAVORS, self.flavor
        coding = [i for i, base in enumerate("".join(self.exons)) if base.isupper()]
        assert coding == list(range(coding[0], coding[-1] + 1)), "the coding bases must form one block"
        introns = self.introns if self.introns is not None else (INTRON,) * (len(self.exons) - 1)
        assert len(introns) == len(self.exons) - 1

    @property
    def layout(self):
        """The chromosome in transcript orientation: flank, exons and introns, flank."""
        introns = self.introns if self.introns is not None else (INTRON,) * (len(self.exons) - 1)
        parts = [self.flanks[0]]
        for exon, intron in zip(self.exons, [*introns, self.flanks[1]]):
            parts += [exon, intron]
        return "".join(parts)

    def exon_spans(self):
        """(start, end) of each exon in the layout, in transcript order."""
        introns = self.introns if self.introns is not None else (INTRON,) * (len(self.exons) - 1)
        spans = []
        start = len(self.flanks[0])
        for exon, intron in zip(self.exons, [*introns, ""]):
            spans.append((start, start + len(exon)))
            start += len(exon) + len(intron)
        return spans

    def chromosome(self, strand):
        layout = self.layout.upper()
        return layout if strand == "+" else reverse_complement(layout)

    def to_genome(self, start, end, strand):
        """0-based half-open genomic interval of the layout interval start to end."""
        length = len(self.layout)
        return (start, end) if strand == "+" else (length - end, length - start)

    def gff3(self, strand):
        """The lines of the GFF3 with the gene, the transcript and its exon, CDS and codon rows."""
        spans = self.exon_spans()
        ensembl = self.flavor == "ensembl"
        tag = f";tag={','.join(self.tags)}" if self.tags else ""
        if ensembl:
            gene = "ID=gene:g1;gene_id=g1;biotype=protein_coding"
            transcript_row = ("mRNA", "ID=transcript:tx1;Parent=gene:g1;transcript_id=tx1;biotype=protein_coding" + tag)
        else:
            gene = "ID=g1;gene_id=g1;" + TYPES
            transcript_row = ("transcript", "ID=tx1;Parent=g1;gene_id=g1;transcript_id=tx1;" + TYPES + tag)
        transcript = ";gene_id=g1;transcript_id=tx1;" + TYPES

        def row(feature, start, end, phase, attributes):
            genomic_start, genomic_end = self.to_genome(start, end, strand)
            columns = [self.contig, "test", feature, genomic_start + 1, genomic_end, ".", strand, phase, attributes]
            return [str(column) for column in columns]

        rows = [
            row("gene", spans[0][0], spans[-1][1], ".", gene),
            row(transcript_row[0], spans[0][0], spans[-1][1], ".", transcript_row[1]),
        ]
        # The coding bases in layout positions, 5' to 3'
        coding = [
            start + i for (start, _), exon in zip(spans, self.exons) for i, base in enumerate(exon) if base.isupper()
        ]
        features = [("CDS", coding)]
        if self.start_codon and not ensembl:
            features.append(("start_codon", coding[:3]))
        if self.stop_codon and not ensembl:
            features.append(("stop_codon", coding[-3:]))
        phase = self.frame
        for number, (start, end) in enumerate(spans, start=1):
            if ensembl:
                exon_attributes = f"Parent=transcript:tx1;Name=e{number};exon_id=e{number};rank={number}"
            else:
                attributes = f"Parent=tx1{transcript};exon_number={number}{tag}"
                exon_attributes = f"ID=exon:tx1:{number};" + attributes
            if self.exon_rows:
                rows.append(row("exon", start, end, ".", exon_attributes))
            for feature, positions in features:
                part = [position for position in positions if start <= position < end]
                if not part:
                    continue
                if ensembl:
                    feature_attributes = "ID=CDS:p1;Parent=transcript:tx1;protein_id=p1"
                else:
                    feature_attributes = f"ID={feature}:tx1:{number};" + attributes
                rows.append(row(feature, part[0], part[-1] + 1, phase if feature == "CDS" else ".", feature_attributes))
                if feature == "CDS":
                    # the phase of the next CDS row: the bases that its first codon still needs
                    phase = (3 - (len(part) - phase) % 3) % 3
        if self.edit_gff3 is not None:
            rows = self.edit_gff3(rows, strand)
        return ["\t".join(columns) for columns in rows]


@dataclasses.dataclass(frozen=True)
class Change:
    """
    A variant, given as a change in transcript orientation: "left[ref>alt]right". left, ref and right must occur
    once in the layout of the transcript, and ref becomes alt. ref is empty for an insertion, and alt for a deletion.

    The VCF record is on the plus strand. An insertion or a deletion gets a padding base on its left in the genome,
    which on the minus strand is the base 3' of the change in the transcript.

    :param vcf_ref: the REF of the VCF record in transcript orientation, if it differs from ref, e.g. for a reference
        mismatch
    :param vcf_alt: the ALT of the VCF record as written, e.g. a symbolic allele. {ref} stands for the REF of the
        record, e.g. for a breakend.
    :param vcf_id: the ID of the VCF record
    :param lower_case: the alleles that the VCF record writes in lower case: "ref", "alt" or both
    """

    text: str
    vcf_ref: str | None = None
    vcf_alt: str | None = None
    vcf_id: str = "var1"
    lower_case: tuple[str, ...] = ()

    def parse(self):
        match = re.fullmatch(r"([ACGTacgt]*)\[([ACGTacgt]*)>([ACGTacgt]*)\]([ACGTacgt]*)", self.text)
        assert match, f"not a change: {self.text}"
        return tuple(part.upper() for part in match.groups())

    def locate(self, transcript):
        """The layout position of the first ref base, or of the base after an insertion."""
        left, ref, _, right = self.parse()
        layout = transcript.layout.upper()
        starts = [m.start() for m in re.finditer(f"(?={left + ref + right})", layout)]
        assert len(starts) == 1, f"{self.text} occurs {len(starts)} times in the layout"
        return starts[0] + len(left)

    def record(self, transcript, strand):
        """(POS, REF, ALT) of the VCF record."""
        _, ref, alt, _ = self.parse()
        start = self.locate(transcript)
        if self.vcf_ref is not None:
            ref = self.vcf_ref.upper()
        genome = transcript.chromosome(strand)
        start, end = transcript.to_genome(start, start + len(ref), strand)
        if strand == "-":
            ref = reverse_complement(ref)
            alt = reverse_complement(alt)
        if not ref or not alt:
            start -= 1
            ref = genome[start] + ref
            alt = genome[start] + alt
        if "ref" in self.lower_case:
            ref = ref.lower()
        if "alt" in self.lower_case:
            alt = alt.lower()
        if self.vcf_alt is not None:
            alt = self.vcf_alt.format(ref=ref)
        return start + 1, ref, alt


@dataclasses.dataclass(frozen=True)
class Mark:
    """
    A mark under the bases start to end of the ref or the alt line of a drawing: char under each base, then the label.
    E.g. the PTC: Mark("alt", 44, 47, "*", "PTC").

    :param line: "ref" or "alt": the line whose transcript positions start and end are
    :param start: the transcript position of the first marked base, 0-based
    :param end: the transcript position after the last marked base
    """

    line: str
    start: int
    end: int
    char: str = "*"
    label: str = ""


@dataclasses.dataclass(frozen=True)
class Span:
    """
    An arrow <---> under the bases start to end of the ref or the alt line of a drawing, with a label, e.g.
    Span("alt", 44, 56, "ptc_to_exon_end = 12"). The arrow covers end - start bases: < stands under the first one and
    > under the last one.

    :param line: "ref" or "alt": the line whose transcript positions start and end are
    """

    line: str
    start: int
    end: int
    label: str


@dataclasses.dataclass(frozen=True)
class Ruler:
    """
    Numbers above the bases at the given positions of the ref or the alt line of a drawing. A number stands left-aligned
    over its base, and the transcript length over the column after the last base.

    :param unit: "tx" for transcript positions, "CDS" for positions from the first coding base of the line, or
        "layout" for layout positions of the ref line (from the first base of the 5' flank, also in introns)
    """

    positions: tuple[int, ...]
    unit: str = "tx"
    line: str = "ref"


@dataclasses.dataclass(frozen=True)
class Layout:
    """A transcript and the expected values of the columns that depend on it only (REF_COLUMNS)."""

    transcript: Transcript
    ref: dict

    def __post_init__(self):
        assert sorted(self.ref) == sorted(REF_COLUMNS), _column_difference(self.ref, REF_COLUMNS)


# Why a variant gives no row ("Technical Notes.md", section "Output columns"), and a VCF without a record
NO_ROW_REASONS = (
    "touches no coding region",
    "symbolic allele",
    "breakend",
    'ALT "." or "*"',
    "REF mismatch",
    "no record",
)


@dataclasses.dataclass(frozen=True)
class NoRow:
    """The expected result of a variant that gives no row, with the reason (one of NO_ROW_REASONS)."""

    reason: str

    def __post_init__(self):
        assert self.reason in NO_ROW_REASONS, self.reason


@dataclasses.dataclass(frozen=True)
class Raises:
    """
    The expected result of input that annotate() rejects: an exception of the type exception, whose message matches
    the regular expression match, as pytest.raises checks it. match is per_strand for a message with genomic
    coordinates.
    """

    exception: type[Exception]
    match: str | PerStrand


def names(*parts):
    """A pattern for a message that holds each of the parts, in any order. Each part starts and ends a word."""
    return "(?s)" + "".join(rf"(?=.*\b{re.escape(part)}\b)" for part in parts)


@dataclasses.dataclass(frozen=True)
class Case:
    """
    One constellation: a transcript layout and a variant, with the expected value of every column of CASE_COLUMNS.

    :param name: the constellation, as the test id
    :param drawing: the transcript and the variant, 5' to 3': the block that render_case() draws, and sentences that
        explain the constellation
    :param change: the variant, or None for a VCF without a record
    :param expected: the expected values, NoRow if the variant gives no row, or Raises for an error of annotate()
    :param equivalent: other descriptions of the same variant. Each gives the same result, except for the columns that
        echo the VCF record.
    :param bug: the id and a short description of a known bug that makes the case fail. The case then runs as a
        strict xfail and keeps the documented or decided result as the expected one.
    :param more_changes: further variants in the same VCF, each with its own vcf_id
    :param more_rows: the expected values of further rows, e.g. of a second transcript or of more_changes. Each row
        gives the columns whose values differ from the first row. The runner compares the rows in the order of
        transcript_id, then variant_id.
    :param reassign_exons: the reassign_exons argument of annotate()
    :param marks: the Mark, Span and Ruler lines of the drawing, under the change
    :param ruler: the Ruler line of the drawing above the ref line
    """

    name: str
    drawing: str
    layout: Layout
    change: Change | None
    expected: dict | NoRow | Raises
    equivalent: tuple[Change, ...] = ()
    bug: str | None = None
    more_changes: tuple[Change, ...] = ()
    more_rows: tuple[dict, ...] = ()
    reassign_exons: bool = False
    marks: tuple[Mark | Span | Ruler, ...] = ()
    ruler: Ruler | None = None

    def __post_init__(self):
        assert self.drawing.strip(), f"{self.name} has no drawing"
        if isinstance(self.expected, dict):
            assert sorted(self.expected) == sorted(CASE_COLUMNS), _column_difference(self.expected, CASE_COLUMNS)
        else:
            assert not self.more_rows, f"{self.name} has more_rows but expects no row"
        for row in self.more_rows:
            assert set(row) <= set(OUTPUT_COLUMN_KINDS), _column_difference(row, OUTPUT_COLUMN_KINDS)


def _column_difference(columns, expected):
    return f"missing: {[c for c in expected if c not in columns]}, unknown: {[c for c in columns if c not in expected]}"


def case_params(cases):
    """
    pytest params (case, change, strand) of each case on each strand, for its change and each equivalent description.
    A case with a known bug is a strict xfail.
    """
    params = []
    for case in cases:
        marks = [pytest.mark.xfail(reason=case.bug, strict=True)] if case.bug else []
        changes = [("vcf", case.change)] + [(f"equivalent{i}", change) for i, change in enumerate(case.equivalent, 1)]
        for label, change in changes:
            for strand, strand_name in [("+", "plus"), ("-", "minus")]:
                params.append(pytest.param(case, change, strand, id=f"{case.name}-{strand_name}-{label}", marks=marks))
    return params


def plain(value):
    """A value of the result table as a plain Python value, with None for a missing value."""
    if isinstance(value, (list, tuple)):
        return type(value)(plain(item) for item in value)
    if isinstance(value, dict):
        return {key: plain(item) for key, item in value.items()}
    if pd.api.types.is_scalar(value) and pd.isna(value):
        return None
    return value.item() if hasattr(value, "item") else value


def write_inputs(transcript, changes, strand, directory):
    """Write the FASTA, GFF3 and VCF file of a run; return their paths and the VCF record of each change."""
    fasta = directory / "genome.fa"
    gff3 = directory / "annotation.gff3"
    vcf = directory / "variants.vcf"
    fasta.write_text(f">{transcript.fasta_contig or transcript.contig}\n{transcript.chromosome(strand)}\n")
    gff3.write_text("##gff-version 3\n" + "\n".join(transcript.gff3(strand)) + "\n")
    records = [change.record(transcript, strand) for change in changes]
    lines = sorted(
        (pos, f"{transcript.contig}\t{pos}\t{change.vcf_id}\t{ref}\t{alt}\t.\t.\t.\n")
        for change, (pos, ref, alt) in zip(changes, records)
    )
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "".join(line for _, line in lines)
    )
    return (str(vcf), str(gff3), str(fasta)), records


def run(case, change, strand, directory, sequences=True):
    """
    Run annotate() on the case, with the given description of its variant and the further variants of the case.
    Return the result table and the VCF record of the given description, None without one.

    :param sequences: the sequences argument of annotate()
    """
    changes = ([] if change is None else [change]) + list(case.more_changes)
    paths, records = write_inputs(case.layout.transcript, changes, strand, directory)
    results = annotate(*paths, reassign_exons=case.reassign_exons, sequences=sequences)
    return results, (records[0] if change is not None else None)


def on_strand(value, strand):
    """The expected value on the strand."""
    return getattr(value, "plus" if strand == "+" else "minus") if isinstance(value, PerStrand) else value


def expected_row(case, strand, more=None):
    """
    The expected value of every output column on the strand, in output order: of the first row, or of the further row
    with the values more, an item of more_rows. SAME_EXONS becomes the transcript_exons of the row.
    """
    values = {**case.layout.ref, **case.expected, **(more or {})}
    row = {column: on_strand(values[column], strand) for column in OUTPUT_COLUMN_KINDS}
    if row["alt_transcript_exons"] is SAME_EXONS:
        row["alt_transcript_exons"] = row["transcript_exons"]
    return row


def check(case, change, strand, directory, sequences=True, columns=OUTPUT_COLUMN_KINDS):
    """
    Run the case and compare the result with its expected values, or check the error that it expects.

    :param sequences: the sequences argument of annotate()
    :param columns: the kind of each column that the result has, in output order. The expected values of the other
        output columns are not checked.
    """
    if isinstance(case.expected, Raises):
        with pytest.raises(case.expected.exception, match=on_strand(case.expected.match, strand)):
            run(case, change, strand, directory, sequences)
        return
    results, record = run(case, change, strand, directory, sequences)

    assert list(results.columns) == list(columns), "the columns are not the expected ones"
    wrong_dtypes = [
        f"{column}: {results[column].dtype}"
        for column, kind in columns.items()
        if results[column].dtype != pd.api.types.pandas_dtype(KIND_DTYPES[kind])
    ]
    assert not wrong_dtypes, f"wrong dtypes: {wrong_dtypes}"
    if isinstance(case.expected, NoRow):
        assert len(results) == 0, f"expected 0 rows, got {len(results)}"
        return

    rows = [expected_row(case, strand, more) for more in ({}, *case.more_rows)]
    if change is not case.change:
        # Another description of the variant: the columns that echo the VCF record take it from the record, unless a
        # further row gives them
        pos, ref, alt = record
        echo = {"ref": ref, "alt": alt, "start": pos - 1, "end": pos - 1 + len(ref)}
        for row, more in zip(rows, ({}, *case.more_rows)):
            row.update({column: value for column, value in echo.items() if column not in more})
    assert len(results) == len(rows), f"expected {len(rows)} rows, got {len(results)}"

    def key(row):
        # variant_id is null for a record without ID
        return row["transcript_id"], row["variant_id"] is not None, row["variant_id"] or ""

    actual_rows = [{column: plain(value) for column, value in results.iloc[i].items()} for i in range(len(results))]
    differences = [
        (f"{key(expected)}: " if len(rows) > 1 else "")
        + f"{column}: expected {expected[column]!r}, got {actual[column]!r}"
        for expected, actual in zip(sorted(rows, key=key), sorted(actual_rows, key=key))
        for column in columns
        if actual[column] != expected[column] or type(actual[column]) is not type(expected[column])
    ]
    assert not differences, "\n".join(differences)


# The drawing of a case. A run of at least ELISION bases that nothing below needs is drawn as "..N..", where N is the
# number of its bases.
ELISION = 10
# The bases drawn on each side of a change
CHANGE_CONTEXT = 6
# The bases drawn on each side of a mark, of each end of a span and of a ruler number
MARK_CONTEXT = 3
# The bases drawn at each end of an exon. The first 3 and the last 3 coding bases are drawn too.
EXON_EDGE = 3
# The bases drawn at each end of a drawn intron, and at the transcript end of a drawn flank
INTRON_EDGE = 6
# The label of a change gives a longer allele as its length
ALLELE = 10


@dataclasses.dataclass
class _Column:
    """A column of the ref and the alt line: a base of the layout, or an inserted base, which has no ref base."""

    region: tuple[str, int]  # ("exon", i), ("intron", i) with i from 0, or ("flank", 0) at 5' and ("flank", 1) at 3'
    ref: str  # the ref base as drawn, or "-" for an inserted base
    alt: str  # the alt base as drawn, or "-" for a deleted base
    pos: int | None = None  # the layout position of the ref base
    tx: int | None = None  # the transcript position of the ref base
    alt_tx: int | None = None  # the transcript position of the alt base
    changed: bool = False
    inserted: bool = False  # a base of an insertion, not an extra alt base of a delins


def _layout_bases(transcript):
    """The region, the transcript position and the base as drawn, of each layout position."""
    introns = transcript.introns if transcript.introns is not None else (INTRON,) * (len(transcript.exons) - 1)
    parts = [("flank", 0, transcript.flanks[0])]
    for i, exon in enumerate(transcript.exons):
        parts += [("exon", i, exon), ("intron", i, introns[i] if i < len(introns) else "")]
    parts.append(("flank", 1, transcript.flanks[1]))
    regions = []
    tx = []
    drawn = []
    length = 0
    for name, i, seq in parts:
        regions += [(name, i)] * len(seq)
        tx += list(range(length, length + len(seq))) if name == "exon" else [None] * len(seq)
        drawn += seq if name == "exon" else seq.lower()
        length += len(seq) if name == "exon" else 0
    return regions, tx, drawn


def _changed_columns(regions, tx, drawn, start, ref, alt):
    """
    The columns of a change from start, with its ref and alt in upper case. The ref and the alt are aligned from the
    left, and the extra alt bases follow the ref bases. These go into the region of the bases on both sides of them, or
    into the exon at an exon edge. An alt base is upper case where it replaces a coding base or lies between coding
    bases.
    """
    sides = [p for p in (start + len(ref) - 1, start + len(ref)) if 0 <= p < len(regions)]
    exon_sides = [p for p in sides if regions[p][0] == "exon"]
    region = regions[exon_sides[0]] if exon_sides else regions[sides[0]]
    upper = bool(exon_sides) and all(drawn[p].isupper() for p in exon_sides if regions[p] == region)
    columns = []
    for j in range(max(len(ref), len(alt))):
        base = alt[j] if j < len(alt) else "-"
        if j < len(ref):
            p = start + j
            coding = regions[p][0] == "exon" and drawn[p].isupper()
            columns.append(_Column(regions[p], drawn[p], base if coding else base.lower(), p, tx[p], changed=True))
        else:
            columns.append(_Column(region, "-", base if upper else base.lower(), changed=True, inserted=not ref))
    return columns


def _touched(change, transcript):
    """The layout positions of the ref of a change, or the two positions on both sides of an insertion."""
    _, ref, _, _ = change.parse()
    start = change.locate(transcript)
    if ref:
        return list(range(start, start + len(ref)))
    return [p for p in (start - 1, start) if 0 <= p < len(transcript.layout)]


def _change_label(change):
    """ "ref>alt" as written in the change, "-" for an empty allele."""
    match = re.fullmatch(r"[ACGTacgt]*\[([ACGTacgt]*)>([ACGTacgt]*)\][ACGTacgt]*", change.text)
    return ">".join(f"{len(allele)} nt" if len(allele) > ALLELE else allele or "-" for allele in match.groups())


def _ruler_label(ruler):
    assert ruler.unit in ("tx", "CDS", "layout") and ruler.line in ("ref", "alt"), ruler
    assert not (ruler.unit == "layout" and ruler.line == "alt"), "layout positions are positions of the ref line"
    return ("alt " if ruler.line == "alt" else "") + ruler.unit


def render(transcript, change=None, marks=(), ruler=None, more_changes=(), equivalent=()):
    """
    The layout block of a drawing: the transcript and the change, 5' to 3' in transcript orientation, as lines without
    indentation. The lines are, in this order:

    - the ruler, if given;
    - the ref line: exons in [ ], | at each exon junction, 5' UTR and 3' UTR in lower case, the coding region in upper
      case, with a space before each codon of the annotated frame and before the 3' UTR. An intron or a flank that the
      change, one of more_changes or an equivalent description touches is drawn in lower case, outside the brackets.
      Other introns and flanks are left out;
    - the alt line, if the change alters the layout. The ref and the alt line are aligned: "-" fills the gap of a
      deletion in the alt line and of an insertion in the ref line, so each column holds one position in both lines.
      An insertion at an exon edge goes into the exon;
    - a line with ^ under the changed bases and "ref>alt", and one more such line for each of more_changes, with its ID;
    - the marks, in their order.

    :param marks: Mark, Span and Ruler lines
    :param ruler: the Ruler line above the ref line
    """
    regions, tx, drawn = _layout_bases(transcript)
    start = len(drawn)
    ref = ""
    alt = ""
    if change is not None:
        _, ref, alt, _ = change.parse()
        start = change.locate(transcript)
    columns = []
    for p in range(len(drawn) + 1):
        if p == start:
            columns += _changed_columns(regions, tx, drawn, start, ref, alt)
        if p < len(drawn) and not start <= p < start + len(ref):
            columns.append(_Column(regions[p], drawn[p], drawn[p], p, tx[p]))

    # The exons are drawn, and each intron or flank that a change touches
    others = [*more_changes, *equivalent]
    shown = {region for region in regions if region[0] == "exon"}
    shown |= {regions[p] for c in [change, *others] if c is not None for p in _touched(c, transcript)}
    columns = [column for column in columns if column.region in shown]
    alt_tx = 0
    for column in columns:
        if column.region[0] == "exon" and column.alt != "-":
            column.alt_tx = alt_tx
            alt_tx += 1
    index = {
        ("layout", "ref"): {column.pos: i for i, column in enumerate(columns) if column.pos is not None},
        ("tx", "ref"): {column.tx: i for i, column in enumerate(columns) if column.tx is not None},
        ("tx", "alt"): {column.alt_tx: i for i, column in enumerate(columns) if column.alt_tx is not None},
    }
    coding = [t for t, base in enumerate("".join(transcript.exons)) if base.isupper()]
    cds_start = coding[0]
    cds_end = coding[-1] + 1
    first_coding = {"ref": cds_start}
    first_coding["alt"] = sum(1 for column in columns[: index["tx", "ref"][cds_start]] if column.alt_tx is not None)

    def find(line, position, unit="tx"):
        """(column, offset) of a position. The transcript length is at offset 1 from the last base."""
        found = index["layout" if unit == "layout" else "tx", line]
        position += first_coding[line] if unit == "CDS" else 0
        if position in found:
            return found[position], 0
        if unit != "layout" and position == len(found):
            return found[position - 1], 1
        raise ValueError(f"the {line} line has no {unit} position {position}")

    # The columns that are drawn in any case
    pinned = set()

    def pin(first, last, context=0):
        pinned.update(range(first - context, last + context + 1))

    changed = [i for i, column in enumerate(columns) if column.changed]
    if changed:
        pin(changed[0], changed[-1], CHANGE_CONTEXT)
    for c in others:
        touched = [index["layout", "ref"][p] for p in _touched(c, transcript)]
        pin(min(touched), max(touched), CHANGE_CONTEXT)
    for region in shown:
        members = [i for i, column in enumerate(columns) if column.region == region]
        edge = EXON_EDGE if region[0] == "exon" else INTRON_EDGE
        if region != ("flank", 0):
            pin(members[0], members[0] + edge - 1)
        if region != ("flank", 1):
            pin(members[-1] - edge + 1, members[-1])
    pin(index["tx", "ref"][cds_start], index["tx", "ref"][min(cds_start + 2, cds_end - 1)])
    pin(index["tx", "ref"][max(cds_end - 3, cds_start)], index["tx", "ref"][cds_end - 1])
    for item in [*marks, *([ruler] if ruler else [])]:
        if isinstance(item, Ruler):
            ends = [find(item.line, position, item.unit)[0] for position in item.positions]
        else:
            assert item.start < item.end, item
            ends = [find(item.line, item.start)[0], find(item.line, item.end - 1)[0]]
        for first, last in [(ends[0], ends[-1])] if isinstance(item, Mark) else zip(ends, ends):
            pin(first, last, MARK_CONTEXT)

    # A space before each codon and before the 3' UTR, unless a bracket stands there. An insertion goes after it.
    groups = {cds_start, *range(cds_start + transcript.frame, cds_end, 3), cds_end}
    opens = [False] * len(columns)
    for i, column in enumerate(columns):
        if column.region[0] == "exon" and column.tx in groups:
            j = i
            while j > 0 and columns[j - 1].inserted and columns[j - 1].region == column.region:
                j -= 1
            opens[j] = j > 0 and columns[j - 1].region == column.region

    # The runs to elide, as first column: column after the run. In the coding region, a run holds whole codons.
    elided = {}
    i = 0
    while i < len(columns):
        end = i
        while end < len(columns) and end not in pinned and columns[end].region == columns[i].region:
            end += 1
        first = i
        last = end
        if any(column.tx is not None and cds_start <= column.tx < cds_end for column in columns[i:end]):
            cuts = [
                j for j in range(i, min(end + 1, len(columns))) if opens[j] and columns[j].region == columns[i].region
            ]
            first, last = (cuts[0], cuts[-1]) if cuts else (i, i)
        if last - first >= ELISION:
            elided[first] = last
        i = max(end, i + 1)

    # The ref and the alt line, and the x of each drawn column in them
    ref_line = "5' "
    alt_line = "5' "
    xs = {}
    previous = None
    i = 0
    while i < len(columns):
        column = columns[i]
        if column.region == previous:
            separator = " " if opens[i] else ""
        else:
            exon_before = previous is not None and previous[0] == "exon"
            exon_now = column.region[0] == "exon"
            separator = "]" * exon_before + "|" * (exon_before and exon_now) + "[" * exon_now
        ref_line += separator
        alt_line += separator
        previous = column.region
        if i in elided:
            marker = f"..{elided[i] - i}.."
            ref_line += marker
            alt_line += marker
            i = elided[i]
            continue
        xs[i] = len(ref_line)
        ref_line += column.ref
        alt_line += column.alt
        i += 1
    closing = "]" * (previous[0] == "exon") + " 3'"
    ref_line += closing
    alt_line += closing

    def x(found):
        column, offset = found
        return xs[column] + offset

    def under(first, last, char, label):
        """char from x first to x last, then the label."""
        return " " * first + char * (last - first + 1) + (f" {label}" if label else "")

    def arrow(first, last, label):
        """An arrow from x first to x last, with the label inside if it fits, else after it."""
        width = last - first + 1
        inner = f" {label} "
        dashes = width - 2 - len(inner)
        if dashes >= 4:
            return " " * first + "<" + "-" * (dashes // 2) + inner + "-" * (dashes - dashes // 2) + ">"
        return " " * first + ("<" + "-" * (width - 2) + ">" if width > 1 else "<>") + f" {label}"

    def numbers(item):
        line = ""
        for position in item.positions:
            at = x(find(item.line, position, item.unit))
            if line and at <= len(line):
                raise ValueError(f"ruler number {position} touches the number before it: drop one of them")
            line += " " * (at - len(line)) + str(position)
        return line

    rows = [(_ruler_label(ruler), numbers(ruler))] if ruler else []
    rows.append(("ref", ref_line))
    if ref != alt:
        rows.append(("alt", alt_line))
    if changed:
        rows.append(("", under(xs[changed[0]], xs[changed[-1]], "^", _change_label(change))))
    for c in more_changes:
        touched = [xs[index["layout", "ref"][p]] for p in _touched(c, transcript)]
        rows.append(("", under(min(touched), max(touched), "^", f"{_change_label(c)} {c.vcf_id}")))
    for item in marks:
        if isinstance(item, Ruler):
            rows.append((_ruler_label(item), numbers(item)))
            continue
        first = x(find(item.line, item.start))
        last = x(find(item.line, item.end - 1))
        rows.append(
            (
                "",
                under(first, last, item.char, item.label) if isinstance(item, Mark) else arrow(first, last, item.label),
            )
        )
    width = max(len(label) for label, _ in rows) + 1
    return "\n".join((label.ljust(width) + body).rstrip() for label, body in rows)


def render_case(case):
    """The layout block of the drawing of a case, which its drawing must hold."""
    return render(case.layout.transcript, case.change, case.marks, case.ruler, case.more_changes, case.equivalent)
