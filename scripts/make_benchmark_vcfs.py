"""
Build resources/TCGA_benchmark/tcga_dataset.vcf and resources/MMRF_benchmark/MMRF_TARGET_dataset.vcf from the tables
of NMDEff (https://github.com/hjkng/nmdeff), with REF and ALT on the forward strand of GRCh38.

Run it with:

    uv run --group scripts python scripts/make_benchmark_vcfs.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa

Each VCF has one row per table row, in table order. ID, QUAL, FILTER and INFO are ".". add_alleles explains how REF
and ALT come from HGVSc. train_model.py and validate_model.py import the shared steps from this file.
"""

import argparse
import hashlib
import logging
import urllib.request
from pathlib import Path

import polars as pl
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
HGVSC_SUBSTITUTION = r"^(?:[^:]+:)?c\.\d+(?P<ref>[ACGT])>(?P<alt>[ACGT])$"
GFF3_COLUMNS = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def check_sha256(data: bytes, expected: str, name: str) -> bytes:
    digest = hashlib.sha256(data).hexdigest()
    if digest != expected:
        raise ValueError(f"{name} has sha256 {digest}, expected {expected}")
    return data


def download(url: str, expected_sha256: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        return check_sha256(response.read(), expected_sha256, url)


def read_nmdeff(name: str) -> pl.DataFrame:
    """Download a table of NMDEff at the pinned commit and check its sha256."""
    return pl.read_csv(download(NMDEFF_URL + name, NMDEFF_SHA256[name]), infer_schema_length=None)


def read_transcripts(gff3: Path) -> pl.DataFrame:
    """
    Read the transcript lines of the GFF3, one row per transcript. The columns are transcript, strand and gene.
    transcript is the transcript ID without version. gene is the ID without version of the gene in the Parent of the
    line. A transcript line is any line whose ID attribute starts with ENST (exon, CDS and UTR lines have IDs like
    exon:ENST...).
    """
    transcripts = (
        pl.scan_csv(
            gff3,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            quote_char=None,
            new_columns=GFF3_COLUMNS,
            infer_schema=False,
        )
        .select(
            transcript=pl.col("attributes").str.extract(r"^ID=(ENST[^.;]*)"),
            strand=pl.col("strand"),
            gene=pl.col("attributes").str.extract(r";Parent=([^;.]+)"),
        )
        .filter(pl.col("transcript").is_not_null())
        .unique()
        .collect()
    )
    if (orphans := transcripts.filter(pl.col("gene").is_null())).height > 0:
        raise ValueError(f"the line of {orphans['transcript'][0]} in {gff3} has no Parent")
    return transcripts


def add_alleles(table: pl.DataFrame, transcripts: pl.DataFrame, fasta: Fasta) -> pl.DataFrame:
    """
    Rename the columns chromosome, start and Transcript_ID of a table of Kim et al. 2024 to chrom, pos and transcript,
    and add REF and ALT on the forward strand. transcripts is the table of read_transcripts.

    HGVSc (c.123G>T, optionally prefixed by the transcript ID) has the alleles on the transcript strand. On a
    minus-strand transcript of the GENCODE GFF3, complement them. For a transcript missing from the GFF3 (the ID of
    the table has no version, the ID of the GFF3 has one), use the orientation whose REF matches the FASTA. Fail if
    HGVSc is not a single-base substitution, or if REF does not match the FASTA.
    """
    table = table.rename({"chromosome": "chrom", "start": "pos", "Transcript_ID": "transcript"})
    table = table.with_columns(
        pl.col("HGVSc").str.extract_groups(HGVSC_SUBSTITUTION).struct.unnest(),
        genomic_ref=pl.Series(
            [fasta[chrom][pos - 1 : pos].seq.upper() for chrom, pos in zip(table["chrom"], table["pos"])]
        ),
        strand=pl.col("transcript")
        .str.split(".")
        .list.first()
        .replace_strict(transcripts["transcript"], transcripts["strand"], default=None),
    )
    if (invalid := table.filter(pl.col("ref").is_null())).height > 0:
        row = invalid.row(0, named=True)
        raise ValueError(f"{row['HGVSc']} ({row['transcript']}) is not a single-base substitution")
    strand = pl.col("strand")
    complement = (strand == "-") | (strand.is_null() & (pl.col("ref") != pl.col("genomic_ref")))
    table = table.with_columns(
        pl.when(complement).then(pl.col("ref", "alt").replace(COMPLEMENT)).otherwise(pl.col("ref", "alt"))
    )
    if (mismatches := table.filter(pl.col("ref") != pl.col("genomic_ref"))).height > 0:
        row = mismatches.row(0, named=True)
        raise ValueError(
            f"{row['chrom']}:{row['pos']} {row['HGVSc']} ({row['transcript']}): REF {row['ref']} but the FASTA has"
            f" {row['genomic_ref']}"
        )
    return table.drop("genomic_ref", "strand")


def join_features(
    measured: pl.DataFrame,
    features: pl.DataFrame,
    transcripts: pl.DataFrame,
    cohort: str,
    logger: logging.Logger,
) -> pl.DataFrame:
    """
    Join each measured row to the NMD-Scanner row of its transcript and variant, and keep the rows with
    nmd_efficiency_prediction_status "ok". Log the rows that this drops, and why, to logger. The rows keep the order of
    measured. transcripts is the table of read_transcripts.

    A measured row has no NMD-Scanner row if the GFF3 lacks its transcript, or if its variant touches no coding region
    of the transcript (see "Technical Notes.md"). A REF that does not match the FASTA would also give no row, but
    train_model.py and validate_model.py fail on it before.
    """
    rows = measured.join(features, on=KEY, validate="m:1", maintain_order="left")
    unmatched = measured.join(features, on=KEY, how="anti")["transcript"]
    in_gff3 = unmatched.is_in(transcripts["transcript"].implode())
    if len(unmatched) > 0:
        logger.info(
            "%s: %d of %d rows have no NMD-Scanner row of their transcript and variant: %d on transcripts missing from"
            " the GFF3, %d with a variant that touches no coding region of its transcript (%s)",
            cohort,
            len(unmatched),
            len(measured),
            (~in_gff3).sum(),
            in_gff3.sum(),
            " ".join(unmatched.filter(in_gff3).unique().sort()) or "none",
        )
    status = rows["nmd_efficiency_prediction_status"]
    dropped = status.filter(status != "ok").value_counts(sort=True)
    if len(dropped) > 0:
        logger.info(
            "%s: %d of the %d joined rows have an nmd_efficiency_prediction_status other than ok: %s",
            cohort,
            dropped["count"].sum(),
            len(rows),
            ", ".join(f"{count} {value}" for value, count in dropped.iter_rows()),
        )
    return rows.filter(status == "ok")


def write_vcf(variants: pl.DataFrame, path: Path) -> None:
    """Write one VCF row per row of variants, from its columns chrom, pos, ref and alt."""
    columns = variants.select(
        pl.col("chrom").alias("#CHROM"),
        POS="pos",
        ID=pl.lit("."),
        REF="ref",
        ALT="alt",
        QUAL=pl.lit("."),
        FILTER=pl.lit("."),
        INFO=pl.lit("."),
    )
    with open(path, "w", newline="") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        columns.write_csv(handle, separator="\t", quote_style="never")


def add_key(features: pl.DataFrame) -> pl.DataFrame:
    """
    Add the columns transcript and pos of KEY to an nmd_scanner.annotate result as a polars table of
    nmd_scanner.to_arrow. NMD-Scanner's start is 0-based, so start + 1 is pos.
    """
    return features.with_columns(
        transcript=pl.col("transcript_id").str.split(".").list.first(),
        pos=pl.col("start") + 1,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, type=Path, help="GENCODE GFF3 of GRCh38, optionally gzip-compressed")
    parser.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA with chr names and a .fai index")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    transcripts = read_transcripts(args.gff3)
    fasta = Fasta(str(args.fasta))
    for name, path in VCFS.items():
        table = add_alleles(read_nmdeff(name), transcripts, fasta)
        write_vcf(table, path)
        logger.info("Wrote %d rows to %s", len(table), path)


if __name__ == "__main__":
    main()
