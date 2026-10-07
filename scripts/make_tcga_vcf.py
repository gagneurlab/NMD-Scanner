# /// script
# requires-python = "==3.12.*"
# dependencies = [
#     "pyfaidx==0.9.0.4",
# ]
# [tool.uv]
# exclude-newer = "2026-10-07T16:00:00Z"
# ///
"""
Build resources/TCGA_benchmark/tcga_dataset.vcf from tcga_dataset.csv of NMDEff, with REF and ALT on the forward strand
of GRCh38.

Run it with:

    uv run scripts/make_tcga_vcf.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa

Steps:
1. Download tcga_dataset.csv of NMDEff (https://github.com/hjkng/nmdeff) at a pinned commit and check its sha256.
   With --csv, read a local copy instead, and check its sha256 the same way.
2. Parse the substitution of each row from HGVSc (c.123G>T, optionally prefixed by the transcript ID). HGVSc has the
   alleles on the transcript strand. Fail on anything that is not a single-base substitution.
3. Read the strand of the transcript from the GENCODE GFF3 (the ID of the table has no version, the ID of the GFF3 has
   one). On a minus-strand transcript, complement REF and ALT. For a transcript missing from the GFF3, use the
   orientation whose REF matches the FASTA.
4. Check that each REF matches the FASTA, and fail otherwise.
5. Write one VCF row per table row, in table order. POS is the start of the table (start equals end for every row).
   ID, QUAL, FILTER and INFO are ".".

The summary on stderr lists the rows, the complemented rows, and the rows whose strand came from the GFF3 or from the
FASTA.
"""

import argparse
import csv
import gzip
import hashlib
import io
import logging
import re
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from pyfaidx import Fasta

logger = logging.getLogger("make_tcga_vcf")

# tcga_dataset.csv was last changed in this commit, which is also the tag v1.0 of NMDEff. Same pin as train_model.py.
TCGA_URL = "https://raw.githubusercontent.com/hjkng/nmdeff/08c92768fcb689236a833db6a2f2d9bcbe919f12/tcga_dataset.csv"
TCGA_SHA256 = "7458af15a699a643de3638114af99a17d7eb95b3f80bce187e45347b9c1dc3ca"
DEFAULT_OUT = Path(__file__).resolve().parent.parent / "resources" / "TCGA_benchmark" / "tcga_dataset.vcf"

HGVSC_SUBSTITUTION = re.compile(r"^(?:[^:]+:)?c\.\d+([ACGT])>([ACGT])$")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True, slots=True)
class TableRow:
    """The columns of one row of tcga_dataset.csv that the VCF needs."""

    chromosome: str
    start: int
    end: int
    transcript_id: str
    hgvsc: str


@dataclass(frozen=True, slots=True)
class VcfRow:
    chrom: str
    pos: int
    ref: str
    alt: str
    complemented: bool
    strand_source: str


def read_table(csv_path: Path | None) -> list[TableRow]:
    if csv_path is None:
        logger.info("Downloading %s", TCGA_URL)
        with urllib.request.urlopen(TCGA_URL, timeout=60) as response:
            data = response.read()
    else:
        data = csv_path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != TCGA_SHA256:
        raise ValueError(f"the table has sha256 {digest}, expected {TCGA_SHA256}")
    return [
        TableRow(
            chromosome=record["chromosome"],
            start=int(record["start"]),
            end=int(record["end"]),
            transcript_id=record["Transcript_ID"],
            hgvsc=record["HGVSc"],
        )
        for record in csv.DictReader(io.StringIO(data.decode("utf-8"), newline=""))
    ]


def read_strands(gff3: Path) -> dict[str, str]:
    """
    Map each transcript ID without version to its strand. A transcript line is any line whose ID attribute starts with
    ENST (exon, CDS and UTR lines have IDs like exon:ENST...).
    """
    strands = {}
    opener = gzip.open if gff3.suffix == ".gz" else open
    with opener(gff3, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            attributes = fields[8]
            if not attributes.startswith("ID=ENST"):
                continue
            transcript = attributes.split(";", 1)[0][len("ID=") :].split(".")[0]
            strands[transcript] = fields[6]
    return strands


def complement(allele: str) -> str:
    return allele.translate(COMPLEMENT)


def build_row(record: TableRow, strands: dict[str, str], fasta: Fasta) -> VcfRow:
    """
    Convert the HGVSc substitution of one table row to genomic alleles. Fail if the substitution is not a single base,
    if the table has start != end, or if neither orientation matches the FASTA.
    """
    hgvsc = record.hgvsc
    match = HGVSC_SUBSTITUTION.match(hgvsc)
    if match is None:
        raise ValueError(f"{hgvsc} ({record.transcript_id}) is not a single-base substitution")
    if record.start != record.end:
        raise ValueError(f"{hgvsc} ({record.transcript_id}) has start {record.start} and end {record.end}")
    chrom = record.chromosome
    pos = record.start
    ref = match.group(1)
    alt = match.group(2)
    genomic_ref = fasta[chrom][pos - 1 : pos].seq.upper()

    strand = strands.get(record.transcript_id.split(".")[0])
    if strand is None:
        strand_source = "fasta"
        complemented = ref != genomic_ref
    elif strand in ("+", "-"):
        strand_source = "gff3"
        complemented = strand == "-"
    else:
        raise ValueError(f"unknown strand {strand} for {record.transcript_id}")
    if complemented:
        ref = complement(ref)
        alt = complement(alt)
    if ref != genomic_ref:
        raise ValueError(f"{chrom}:{pos} {hgvsc} ({record.transcript_id}): REF {ref} but the FASTA has {genomic_ref}")
    return VcfRow(chrom, pos, ref, alt, complemented, strand_source)


def write_vcf(rows: list[VcfRow], path: Path) -> None:
    with open(path, "w", newline="") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            handle.write(f"{row.chrom}\t{row.pos}\t.\t{row.ref}\t{row.alt}\t.\t.\t.\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, type=Path, help="GENCODE GFF3 of GRCh38, optionally gzip-compressed")
    parser.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA with chr names and a .fai index")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output VCF (default: %(default)s)")
    parser.add_argument("--csv", type=Path, help="local copy of tcga_dataset.csv; its sha256 is still checked")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    table = read_table(args.csv)
    strands = read_strands(args.gff3)
    fasta = Fasta(str(args.fasta))
    rows = [build_row(record, strands, fasta) for record in table]
    write_vcf(rows, args.out)
    from_gff3 = sum(row.strand_source == "gff3" for row in rows)
    from_fasta = sum(row.strand_source == "fasta" for row in rows)
    logger.info("rows: %d", len(rows))
    logger.info("complemented rows: %d", sum(row.complemented for row in rows))
    logger.info(
        "strand from the GFF3: %d rows, from the FASTA (transcript not in the GFF3): %d rows", from_gff3, from_fasta
    )
    logger.info("wrote %s", args.out)


if __name__ == "__main__":
    main()
