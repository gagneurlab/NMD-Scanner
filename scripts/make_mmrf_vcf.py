# /// script
# requires-python = "==3.12.*"
# dependencies = [
#     "pyfaidx==0.9.0.4",
# ]
# [tool.uv]
# exclude-newer = "2026-10-07T16:00:00Z"
# ///
"""
Build resources/MMRF_benchmark/MMRF_TARGET_dataset.vcf from MMRF_TARGET_dataset.csv of NMDEff, with REF and ALT on the
forward strand of GRCh38. It imports the steps of make_tcga_vcf.py, which must be in the same folder.

Run it with:

    uv run scripts/make_mmrf_vcf.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa

Steps:
1. Download MMRF_TARGET_dataset.csv of NMDEff (https://github.com/hjkng/nmdeff) at the commit that make_tcga_vcf.py
   pins, and check its sha256. With --csv, read a local copy instead, and check its sha256 the same way.
2. Build each VCF row with build_row of make_tcga_vcf.py: parse HGVSc, take the strand from the GENCODE GFF3 or the
   FASTA, complement REF and ALT on the minus strand, and check REF against the FASTA.
3. Write one VCF row per table row, in table order. POS is the start of the table. The table has no end column.
   ID, QUAL, FILTER and INFO are ".".

The summary on stderr lists the rows, the complemented rows, and the rows whose strand came from the GFF3 or from the
FASTA.
"""

import argparse
import csv
import hashlib
import io
import logging
import urllib.request
from pathlib import Path

from make_tcga_vcf import TableRow, build_row, read_strands, write_vcf
from pyfaidx import Fasta

logger = logging.getLogger("make_mmrf_vcf")

# The same commit as TCGA_URL in make_tcga_vcf.py, which is also the tag v1.0 of NMDEff.
MMRF_URL = (
    "https://raw.githubusercontent.com/hjkng/nmdeff/08c92768fcb689236a833db6a2f2d9bcbe919f12/MMRF_TARGET_dataset.csv"
)
MMRF_SHA256 = "c063914da054292c6c6b867eb58650d5378c74d5c852f54c874eab8e51cb4a19"
DEFAULT_OUT = Path(__file__).resolve().parent.parent / "resources" / "MMRF_benchmark" / "MMRF_TARGET_dataset.vcf"


def read_mmrf_table(csv_path: Path | None) -> list[TableRow]:
    if csv_path is None:
        logger.info("Downloading %s", MMRF_URL)
        with urllib.request.urlopen(MMRF_URL, timeout=60) as response:
            data = response.read()
    else:
        data = csv_path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != MMRF_SHA256:
        raise ValueError(f"the table has sha256 {digest}, expected {MMRF_SHA256}")
    return [
        TableRow(
            chromosome=record["chromosome"],
            start=int(record["start"]),
            # The table has no end column. A single-base substitution ends where it starts.
            end=int(record["start"]),
            transcript_id=record["Transcript_ID"],
            hgvsc=record["HGVSc"],
        )
        for record in csv.DictReader(io.StringIO(data.decode("utf-8"), newline=""))
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, type=Path, help="GENCODE GFF3 of GRCh38, optionally gzip-compressed")
    parser.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA with chr names and a .fai index")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output VCF (default: %(default)s)")
    parser.add_argument("--csv", type=Path, help="local copy of MMRF_TARGET_dataset.csv; its sha256 is still checked")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    table = read_mmrf_table(args.csv)
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
