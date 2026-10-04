"""
Place a variant relative to the exon boundaries of a transcript.

A VCF record names one placement of a variant. An indel in a repeat has more placements that give the same alt
sequence, and near an exon boundary they can disagree on whether the exon or the intron changes. A delins whose REF
and ALT differ in length has two, because its bases can be matched from either end. The functions here use all
placements and the splice dinucleotides to find the alt coding region, or to report why it is unknown.

Coordinates are 0-based, half-open and on the plus strand. A boundary is a point between two bases: boundary x lies
between base x - 1 and base x.
"""

from collections import OrderedDict
from dataclasses import dataclass, field

from Bio.Seq import Seq

# Reasons for an unknown alt transcript
SPLICE_SITE_DESTROYED = "splice_site_destroyed"
EXON_BOUNDARY_AMBIGUOUS = "exon_boundary_ambiguous"

_CHUNK_SIZE = 1024
# A ReferenceSequence keeps at most this many chunks, about 16 MB of bases. Without a bound, a whole-genome VCF would
# keep most of the genome in memory.
_MAX_CHUNKS = 16384


class ReferenceSequence:
    """
    Upper case bases of one chromosome of a pyfaidx.Fasta, fetched in chunks. It keeps the most recently used
    ``max_chunks`` chunks.
    """

    def __init__(self, fasta, chromosome, max_chunks=_MAX_CHUNKS):
        self._sequence = fasta[chromosome]
        self._chunks = OrderedDict()
        self._max_chunks = max_chunks

    def _chunk(self, index):
        if index in self._chunks:
            self._chunks.move_to_end(index)
            return self._chunks[index]
        start = index * _CHUNK_SIZE
        chunk = self._chunks[index] = self._sequence[start : start + _CHUNK_SIZE].seq.upper()
        if len(self._chunks) > self._max_chunks:
            self._chunks.popitem(last=False)
        return chunk

    def base(self, position):
        """The base at ``position``, or None outside the chromosome."""
        if position < 0:
            return None
        chunk = self._chunk(position // _CHUNK_SIZE)
        offset = position % _CHUNK_SIZE
        return chunk[offset] if offset < len(chunk) else None

    def bases(self, start, end):
        """The bases in [start, end), cut at the chromosome ends."""
        start = max(start, 0)
        if end <= start:
            return ""
        return "".join(
            self._chunk(index)[max(start - index * _CHUNK_SIZE, 0) : end - index * _CHUNK_SIZE]
            for index in range(start // _CHUNK_SIZE, (end - 1) // _CHUNK_SIZE + 1)
        )


@dataclass(frozen=True)
class Placement:
    """
    One placement of a variant: the reference bases [start, end) become ``alt``.
    An insertion has start == end; its bases go between base start - 1 and base start.
    ``match_left`` says how the changed bases are matched (map_boundary).
    """

    start: int
    end: int
    alt: str
    match_left: bool = True

    @property
    def length_change(self):
        return len(self.alt) - (self.end - self.start)

    def apply(self, bases, offset):
        """The alt version of ``bases``, which start at ``offset`` and cover the placement."""
        return bases[: self.start - offset] + self.alt + bases[self.end - offset :]

    def map_boundary(self, boundary, insert_left):
        """
        The alt position of a boundary. If ``match_left``, changed bases are matched from the left: reference base
        start + i becomes alt base start + i. Otherwise they are matched from the right: reference base end - 1 - i
        becomes alt base start + len(alt) - 1 - i. Both give the same position unless REF and ALT differ in length
        and both are non-empty. An insertion exactly at the boundary goes left of it if ``insert_left``, else right
        of it.
        """
        if self.start == self.end == boundary:
            return boundary + len(self.alt) if insert_left else boundary
        if boundary <= self.start:
            return boundary
        if boundary >= self.end:
            return boundary + self.length_change
        if self.match_left:
            return self.start + min(boundary - self.start, len(self.alt))
        return self.start + max(len(self.alt) - (self.end - boundary), 0)


@dataclass(frozen=True)
class ExonBoundary:
    """
    An exon edge of a transcript. ``exon_on_left`` is True at the exon end and False at the exon start.
    ``dinucleotide`` holds the two reference intron bases next to the boundary, the splice donor or acceptor
    site. It is None at the transcript start and end, which are no splice sites.
    """

    position: int
    exon_on_left: bool
    dinucleotide: str | None

    def alt_position(self, placement):
        """
        The alt position a placement gives this boundary (Placement.map_boundary). An insertion at the boundary goes
        into the exon. A change across the boundary maps it base for base, so the boundary does not move to a splice
        dinucleotide elsewhere in the alt bases.
        """
        return placement.map_boundary(self.position, insert_left=self.exon_on_left)

    def is_valid(self, alt_bases, offset, position):
        """
        Whether the reference splice dinucleotide stays next to the alt ``position``, on the intron side.
        ``alt_bases`` start at ``offset``. Every position is valid at a transcript end.
        """
        if self.dinucleotide is None:
            return True
        index = position - offset
        if self.exon_on_left:
            return alt_bases[index : index + 2] == self.dinucleotide
        return index >= 2 and alt_bases[index - 2 : index] == self.dinucleotide


@dataclass
class TranscriptEffect:
    """
    Effect of a variant on the coding region of one transcript.

    ``unknown_reason`` is SPLICE_SITE_DESTROYED or EXON_BOUNDARY_AMBIGUOUS if the alt transcript is unknown,
    else None. ``alt_coding`` maps each given coding row (start, end) to its alt bases, if the alt is known.
    ``utr5`` and ``utr3`` give the change of the UTR next to the coding region as (ref, alt), in transcript
    orientation: the ref bases right before the coding region (5'UTR) or right after it (3'UTR), and the alt bases
    that replace them. If the variant leaves that UTR unchanged, ref equals alt, and both are mostly empty.
    """

    unknown_reason: str | None = None
    alt_coding: dict = field(default_factory=dict)
    utr5: tuple[str, str] = ("", "")
    utr3: tuple[str, str] = ("", "")


def trim_alleles(start, ref, alt):
    """
    Remove the bases that REF and ALT share at the start, then those they share at the end.

    :param start: 0-based position of the first REF base
    :return: (start, ref, alt) of the changed bases. ref is empty for an insertion, alt for a deletion.
    """
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    ref, alt = ref[prefix:], alt[prefix:]
    suffix = 0
    while suffix < min(len(ref), len(alt)) and ref[-1 - suffix] == alt[-1 - suffix]:
        suffix += 1
    return start + prefix, ref[: len(ref) - suffix], alt[: len(alt) - suffix]


def equivalent_placements(start, ref, alt, reference):
    """
    All placements of a trimmed variant that give the same alt sequence, from left to right by where their length
    change lies. A deletion or insertion shifts through a repeat. A delins whose REF and ALT differ in length has two
    placements: its bases matched from the right, then from the left. Any other variant has one placement.

    :param start, ref, alt: the changed bases, as returned by trim_alleles
    :param reference: ReferenceSequence of the chromosome
    """
    if ref and alt:
        if len(ref) == len(alt):
            return [Placement(start, start + len(ref), alt)]
        # Where a boundary inside REF falls in ALT is unknown, so both ways to match the bases are placements. With
        # matching from the right, the length change lies at the left end, so that placement comes first.
        return [Placement(start, start + len(ref), alt, match_left=False), Placement(start, start + len(ref), alt)]

    if ref:
        length = len(ref)
        first = start
        while first > 0 and reference.base(first - 1) == reference.base(first + length - 1):
            first -= 1
        last = start
        while reference.base(last + length) is not None and reference.base(last) == reference.base(last + length):
            last += 1
        return [Placement(position, position + length, "") for position in range(first, last + 1)]

    inserted, position = alt, start
    while position > 0 and reference.base(position - 1) == inserted[-1]:
        inserted, position = inserted[-1] + inserted[:-1], position - 1
    placements = [Placement(position, position, inserted)]
    while reference.base(position) == inserted[0]:
        inserted, position = inserted[1:] + inserted[0], position + 1
        placements.append(Placement(position, position, inserted))
    return placements


def variant_placements(start, ref, alt, reference):
    """
    The placements of a VCF record. Its ALT is a sequence, not a symbolic allele or a breakend.

    :param start: 0-based position of the first REF base
    :return: list of Placement, empty if REF and ALT are equal
    """
    start, ref, alt = trim_alleles(start, ref.upper(), alt.upper())
    if not ref and not alt:
        return []
    return equivalent_placements(start, ref, alt, reference)


def exon_boundaries(exons, reference):
    """
    The ExonBoundary at both edges of each exon of one transcript.

    :param exons: (start, end) of each exon
    """
    if not exons:
        return []
    transcript_start = min(start for start, _ in exons)
    transcript_end = max(end for _, end in exons)
    boundaries = []
    for start, end in exons:
        acceptor = None if start == transcript_start else reference.bases(start - 2, start)
        donor = None if end == transcript_end else reference.bases(end, end + 2)
        boundaries += [ExonBoundary(start, False, acceptor), ExonBoundary(end, True, donor)]
    return boundaries


def _exon_bases(bases, offset, exons, start, end):
    """
    The bases of the exons in [start, end), joined in genomic order. ``bases`` start at ``offset`` and cover
    [start, end). An exon whose end lies before its start, e.g. in the alt of a deleted exon, adds no bases.
    """
    parts = sorted((max(exon_start, start), min(exon_end, end)) for exon_start, exon_end in exons)
    return "".join(
        bases[part_start - offset : part_end - offset] for part_start, part_end in parts if part_start < part_end
    )


def _touches(placement, start, end, insertion_at_start=False, insertion_at_end=False):
    """
    Whether a placement changes a base in [start, end). An insertion at an edge counts only if the matching
    flag is set, i.e. if its bases go into the interval there.
    """
    if placement.start == placement.end:
        position = placement.start
        return (
            start < position < end
            or (position == start and insertion_at_start)
            or (position == end and insertion_at_end)
        )
    return placement.start < end and placement.end > start


def _start_codon_position(placements, edge, strand, reference, exon_edge):
    """
    The alt position of the coding region's 5' edge. A scanning ribosome starts at the first ATG, so this is the
    5'-most position, over all placements, where an ATG starts in the alt (in transcript orientation). Without such a
    position, it is the position of the placement shifted farthest into the 5'UTR.

    :param exon_edge: alt position of the exon boundary at the edge, or None if the edge lies inside an exon. The
        coding region does not start before it, in transcript orientation.
    """
    plus = strand == "+"
    candidates = {placement.map_boundary(edge, insert_left=plus) for placement in placements}
    if exon_edge is not None:
        candidates = {m for m in candidates if (m >= exon_edge if plus else m <= exon_edge)}
    # The alt bases reach 3 bases past every candidate
    offset = max(min(edge, min(placement.start for placement in placements)) - 3, 0)
    end = max(edge, max(placement.end for placement in placements)) + 3
    alt_bases = placements[0].apply(reference.bases(offset, end), offset)
    if plus:
        with_atg = [m for m in candidates if alt_bases[m - offset : m - offset + 3] == "ATG"]
        if with_atg:
            return min(with_atg)
        position = placements[0].map_boundary(edge, insert_left=True)
        return position if exon_edge is None else max(position, exon_edge)
    with_atg = [m for m in candidates if m - offset >= 3 and alt_bases[m - offset - 3 : m - offset] == "CAT"]
    if with_atg:
        return max(with_atg)
    position = placements[-1].map_boundary(edge, insert_left=False)
    return position if exon_edge is None else min(position, exon_edge)


def _stop_codon_position(placements, edge, strand, exon_edge):
    """
    The alt position of the coding region's 3' edge: the position of the placement shifted farthest into the 3'UTR,
    so an indel that fits into the UTR leaves the coding region unchanged.

    :param exon_edge: alt position of the exon boundary at the edge, or None if the edge lies inside an exon. The
        coding region does not end after it, in transcript orientation.
    """
    if strand == "+":
        position = placements[-1].map_boundary(edge, insert_left=False)
        return position if exon_edge is None else min(position, exon_edge)
    position = placements[0].map_boundary(edge, insert_left=True)
    return position if exon_edge is None else max(position, exon_edge)


def place_in_transcript(placements, coding_rows, exons, reference, strand, coding_region):
    """
    Effect of a variant on the coding region of one transcript.

    The window of the variant is the union of its placements. If the window touches a coding row, or the splice
    dinucleotide at an exon edge of a coding row (not at a transcript end), the variant has an effect.

    Each placement maps each nearby exon boundary to an alt position (ExonBoundary.alt_position), and a position is
    valid if the reference splice dinucleotide stays next to it. A placement is valid if all its positions are valid.
    Positions of different placements are not mixed, because that could make two exons overlap. Without a valid
    placement, the splice site is destroyed. If the valid placements put a boundary at different positions, the
    boundary is ambiguous. Both make the alt transcript unknown. Several placements come from an indel in a repeat, or
    from a delins whose REF and ALT differ in length. A transcript end has no splice dinucleotide, so every position
    is valid there: valid placements that put a transcript end at different alt positions make the boundary
    ambiguous.

    The exon boundaries decide which bases the mRNA holds. A coding row edge at an exon boundary takes the boundary's
    alt position, except at the two edges of the whole coding region, where the exon can hold UTR bases next to the
    coding region. At the stop codon, the edge takes the position of the valid placement shifted farthest into the
    3'UTR, so an indel that fits into the UTR leaves the coding region unchanged. At the start codon, it takes the
    5'-most position, over all valid placements, where an ATG starts in the alt, because a scanning ribosome starts at
    the first ATG. Without such a position, it takes the position of the valid placement shifted farthest into the
    5'UTR. Of the two placements of a delins, the one with its length change at the UTR end counts as shifted into
    the UTR. If such an edge is also an exon edge, its position stays on the exon side of the exon boundary's alt
    position. So for an insertion at that exon boundary, the splice site decides whether its bases enter the mRNA, and
    the coding region edge rule decides whether they are coding. The alt bases that these rules leave outside the
    coding region, and the deleted UTR bases, give the UTR change (TranscriptEffect.utr5 and utr3).

    :param placements: equivalent placements of the variant (equivalent_placements)
    :param coding_rows: (start, end) of the coding rows (CDS plus stop codon, one per exon) near the variant
    :param exons: (start, end) of every exon of the transcript
    :param reference: ReferenceSequence of the chromosome
    :param strand: strand of the transcript, "+" or "-"
    :param coding_region: (start, end) of the transcript's whole coding region: the smallest coding row start and the
        largest coding row end
    :return: TranscriptEffect, or None if the variant touches no coding row and no splice dinucleotide
    """
    if not placements:
        return None

    boundaries = exon_boundaries(exons, reference)
    exon_starts = {b.position for b in boundaries if not b.exon_on_left}
    exon_ends = {b.position for b in boundaries if b.exon_on_left}
    splice_sites = {(b.position, b.exon_on_left) for b in boundaries if b.dinucleotide is not None}
    coding_start, coding_end = coding_region

    # An insertion at an exon edge of a coding row goes into the coding row, but not at an edge of the whole coding
    # region: there its bases go into the UTR. At the start codon they can be coding, but only through another
    # placement inside the coding row.
    targets = []
    for start, end in coding_rows:
        targets.append(
            (start, end, start in exon_starts and start != coding_start, end in exon_ends and end != coding_end)
        )
        if (start, False) in splice_sites:
            targets.append((start - 2, start, False, False))
        if (end, True) in splice_sites:
            targets.append((end, end + 2, False, False))
    if not any(_touches(placement, *target) for placement in placements for target in targets):
        return None

    window_start = min(placement.start for placement in placements)
    window_end = max(placement.end for placement in placements)
    offset = max(window_start - 4, 0)
    alt_bases = placements[0].apply(reference.bases(offset, window_end + 4), offset)

    # A placement is valid if it keeps the splice site at every nearby boundary. Taking each boundary from another
    # placement could make exons overlap.
    nearby = [b for b in boundaries if window_start - 2 <= b.position <= window_end + 2]
    placements = [
        placement
        for placement in placements
        if all(b.is_valid(alt_bases, offset, b.alt_position(placement)) for b in nearby)
    ]
    if not placements:
        return TranscriptEffect(unknown_reason=SPLICE_SITE_DESTROYED)
    positions = {tuple(b.alt_position(placement) for b in nearby) for placement in placements}
    if len(positions) > 1:
        return TranscriptEffect(unknown_reason=EXON_BOUNDARY_AMBIGUOUS)
    resolved = {(b.position, b.exon_on_left): position for b, position in zip(nearby, positions.pop())}

    # An exon edge outside the window is not resolved; every placement maps it to the same position
    leftmost, rightmost = placements[0], placements[-1]
    plus = strand == "+"
    alt_coding = {}
    # Alt position of each exon: the valid placements agree on every boundary, and they all agree outside the window
    alt_exons = [
        (
            ExonBoundary(start, False, None).alt_position(leftmost),
            ExonBoundary(end, True, None).alt_position(leftmost),
        )
        for start, end in exons
    ]
    left = right = ("", "")
    for start, end in coding_rows:
        exon_start, exon_end = resolved.get((start, False)), resolved.get((end, True))
        if start == coding_start:
            alt_start = (
                _start_codon_position(placements, start, strand, reference, exon_start)
                if plus
                else _stop_codon_position(placements, start, strand, exon_start)
            )
        else:
            alt_start = exon_start if exon_start is not None else leftmost.map_boundary(start, insert_left=True)
        if end == coding_end:
            alt_end = (
                _stop_codon_position(placements, end, strand, exon_end)
                if plus
                else _start_codon_position(placements, end, strand, reference, exon_end)
            )
        else:
            alt_end = exon_end if exon_end is not None else rightmost.map_boundary(end, insert_left=False)
        # Guard: exon edges come from one placement and keep their order, but the coding region edge rules can take
        # the two edges of a row from different placements
        if alt_end < alt_start:
            return TranscriptEffect(unknown_reason=SPLICE_SITE_DESTROYED)
        bases_start = min(start, window_start)
        bases = placements[0].apply(reference.bases(bases_start, max(end, window_end)), bases_start)
        alt_coding[start, end] = bases[alt_start - bases_start : alt_end - bases_start]

        # The UTR bases next to the coding region change if the window reaches past its edge, or if the edge rule
        # puts alt bases outside it. Outside [bases_start, max(end, window_end)), the alt equals the reference.
        if start == coding_start:
            ref_bases = _exon_bases(reference.bases(bases_start, start), bases_start, exons, bases_start, start)
            alt_bases = _exon_bases(bases, bases_start, alt_exons, bases_start, alt_start)
            left = (ref_bases, alt_bases)
        if end == coding_end:
            side_end = max(end, window_end)
            alt_side_end = side_end + leftmost.length_change
            ref_bases = _exon_bases(reference.bases(end, side_end), end, exons, end, side_end)
            alt_bases = _exon_bases(bases, bases_start, alt_exons, alt_end, alt_side_end)
            right = (ref_bases, alt_bases)

    # In transcript orientation, the left side is the 5'UTR on the plus strand and the 3'UTR on the minus strand
    if plus:
        return TranscriptEffect(alt_coding=alt_coding, utr5=left, utr3=right)
    utr5 = tuple(str(Seq(side).reverse_complement()) for side in right)
    utr3 = tuple(str(Seq(side).reverse_complement()) for side in left)
    return TranscriptEffect(alt_coding=alt_coding, utr5=utr5, utr3=utr3)
