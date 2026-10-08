# /// script
# requires-python = "==3.12.*"
# dependencies = [
#     "pandas==3.0.6",
#     "pyfaidx==0.9.0.4",
# ]
# [tool.uv]
# exclude-newer = "2026-10-07T16:00:00Z"
# ///
"""
Build resources/TCGA_benchmark/tcga_dataset.vcf and resources/MMRF_benchmark/MMRF_TARGET_dataset.vcf from the tables
of NMDEff (https://github.com/hjkng/nmdeff), with REF and ALT on the forward strand of GRCh38.

Run it with:

    uv run scripts/make_benchmark_vcfs.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa

Each VCF has one row per table row, in table order. ID, QUAL, FILTER and INFO are ".". add_alleles explains how REF
and ALT come from HGVSc. train_model.py and validate_model.py import the shared steps from this file.
"""

import argparse
import gzip
import hashlib
import io
import logging
import re
import urllib.request
from pathlib import Path

import pandas as pd
from pyfaidx import Fasta

logger = logging.getLogger("make_benchmark_vcfs")

RESOURCES = Path(__file__).resolve().parent.parent / "resources"
# Both tables were last changed in this commit, which is also the tag v1.0 of NMDEff
NMDEFF_URL = "https://raw.githubusercontent.com/hjkng/nmdeff/08c92768fcb689236a833db6a2f2d9bcbe919f12/"
NMDEFF_SHA256 = {
    "tcga_dataset.csv": "7458af15a699a643de3638114af99a17d7eb95b3f80bce187e45347b9c1dc3ca",
    "MMRF_TARGET_dataset.csv": "c063914da054292c6c6b867eb58650d5378c74d5c852f54c874eab8e51cb4a19",
}
VCFS = {
    "tcga_dataset.csv": RESOURCES / "TCGA_benchmark" / "tcga_dataset.vcf",
    "MMRF_TARGET_dataset.csv": RESOURCES / "MMRF_benchmark" / "MMRF_TARGET_dataset.vcf",
}

# The columns that join a measured row to the NMD-Scanner row of its transcript and variant
KEY = ["transcript", "chrom", "pos", "ref", "alt"]
HGVSC_SUBSTITUTION = re.compile(r"^(?:[^:]+:)?c\.\d+([ACGT])>([ACGT])$")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def check_sha256(data: bytes, expected: str, name: str) -> bytes:
    digest = hashlib.sha256(data).hexdigest()
    if digest != expected:
        raise ValueError(f"{name} has sha256 {digest}, expected {expected}")
    return data


def download(url: str, expected_sha256: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        return check_sha256(response.read(), expected_sha256, url)


def read_nmdeff(name: str) -> pd.DataFrame:
    """Download a table of NMDEff at the pinned commit and check its sha256."""
    return pd.read_csv(io.BytesIO(download(NMDEFF_URL + name, NMDEFF_SHA256[name])))


def read_strands(gff3: Path) -> dict[str, str]:
    """
    Map each transcript ID without version to its strand. A transcript line is any line whose ID attribute starts with
    ENST (exon, CDS and UTR lines have IDs like exon:ENST...).
    """
    strands = {}
    with gzip.open(gff3, "rt") if gff3.suffix == ".gz" else open(gff3) as handle:
        for line in handle:
            fields = line.split("\t")
            if not line.startswith("#") and fields[8].startswith("ID=ENST"):
                strands[fields[8].split(";", 1)[0].removeprefix("ID=").split(".")[0]] = fields[6]
    return strands


def add_alleles(table: pd.DataFrame, strands: dict[str, str], fasta: Fasta) -> pd.DataFrame:
    """
    Rename the columns chromosome, start and Transcript_ID of a table of Kim et al. 2024 to chrom, pos and transcript,
    and add REF and ALT on the forward strand.

    HGVSc (c.123G>T, optionally prefixed by the transcript ID) has the alleles on the transcript strand. On a
    minus-strand transcript of the GENCODE GFF3, complement them. For a transcript missing from the GFF3 (the ID of
    the table has no version, the ID of the GFF3 has one), use the orientation whose REF matches the FASTA. Fail if
    HGVSc is not a single-base substitution, or if REF does not match the FASTA.
    """
    table = table.rename(columns={"chromosome": "chrom", "start": "pos", "Transcript_ID": "transcript"})
    refs, alts = [], []
    for chrom, pos, transcript, hgvsc in zip(
        table["chrom"], table["pos"], table["transcript"], table["HGVSc"], strict=True
    ):
        match = HGVSC_SUBSTITUTION.match(hgvsc)
        if match is None:
            raise ValueError(f"{hgvsc} ({transcript}) is not a single-base substitution")
        ref, alt = match.groups()
        genomic_ref = fasta[chrom][pos - 1 : pos].seq.upper()
        strand = strands.get(transcript.split(".")[0])
        if strand == "-" or (strand is None and ref != genomic_ref):
            ref, alt = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
        if ref != genomic_ref:
            raise ValueError(f"{chrom}:{pos} {hgvsc} ({transcript}): REF {ref} but the FASTA has {genomic_ref}")
        refs.append(ref)
        alts.append(alt)
    return table.assign(ref=refs, alt=alts)


def write_vcf(variants: pd.DataFrame, path: Path) -> None:
    """Write one VCF row per row of variants, from its columns chrom, pos, ref and alt."""
    with open(path, "w", newline="") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in variants[["chrom", "pos", "ref", "alt"]].itertuples(index=False):
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")


def add_key(features: pd.DataFrame) -> pd.DataFrame:
    """
    Add the columns of KEY to the result of nmd_scanner.annotate. NMD-Scanner's start is 0-based, so start + 1 is pos.
    """
    return features.assign(
        transcript=features["transcript_id"].astype(str).str.split(".").str[0],
        chrom=features["chrom"].astype(str),
        pos=features["start"].astype("int64") + 1,
        ref=features["ref"].astype(str),
        alt=features["alt"].astype(str),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, type=Path, help="GENCODE GFF3 of GRCh38, optionally gzip-compressed")
    parser.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA with chr names and a .fai index")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    strands = read_strands(args.gff3)
    fasta = Fasta(str(args.fasta))
    for name, path in VCFS.items():
        table = add_alleles(read_nmdeff(name), strands, fasta)
        write_vcf(table, path)
        logger.info("Wrote %d rows to %s", len(table), path)


if __name__ == "__main__":
    main()
