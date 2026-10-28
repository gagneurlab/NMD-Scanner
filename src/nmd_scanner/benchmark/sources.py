"""
The sources of the benchmarks, and the steps that the dataset modules share to turn a study table into a benchmark
table (see nmd_scanner.benchmark).
"""

import hashlib
import io
import logging
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import polars as pl
from pyfaidx import Fasta

logger = logging.getLogger(__name__)

# Both tables of NMDEff (https://github.com/hjkng/nmdeff) were last changed in this commit, which is also its tag v1.0
NMDEFF_URL = "https://raw.githubusercontent.com/hjkng/nmdeff/08c92768fcb689236a833db6a2f2d9bcbe919f12/"
# The columns of every benchmark table, in order. tissue, if a table has it, and source_row follow them.
COLUMNS = ["chrom", "pos", "ref", "alt", "transcript", "gene", "NMD_efficiency"]
HGVSC_SUBSTITUTION = r"^(?:[^:]+:)?c\.\d+(?P<ref>[ACGT])>(?P<alt>[ACGT])$"
GFF3_COLUMNS = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def check_sha256(data: bytes, expected: str, name: str) -> bytes:
    """Return data if its sha256 is expected, and fail otherwise. name names data in the error."""
    digest = hashlib.sha256(data).hexdigest()
    if digest != expected:
        raise ValueError(f"{name} has sha256 {digest}, expected {expected}")
    return data


class RemoteFile:
    """
    A read-only file at a URL that answers HTTP range requests. Each read is one request, so zipfile can read single
    members of a remote zip.
    """

    def __init__(self, url: str, size: int) -> None:
        self.url = url
        self.size = size
        self.position = 0

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.position

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        self.position = offset + {io.SEEK_SET: 0, io.SEEK_CUR: self.position, io.SEEK_END: self.size}[whence]
        return self.position

    def read(self, n: int = -1) -> bytes:
        end = self.size if n < 0 else min(self.position + n, self.size)
        if end <= self.position:
            return b""
        request = urllib.request.Request(self.url, headers={"Range": f"bytes={self.position}-{end - 1}"})
        with urllib.request.urlopen(request, timeout=120) as response:
            # A server that ignores the range answers 200 and sends the whole file
            if response.status != 206:
                raise ValueError(f"{self.url} answered a range request with HTTP {response.status}")
            data = response.read()
        self.position += len(data)
        return data


def fetch(downloads: Path, url: str, sha256: str, member: str | None = None, size: int | None = None) -> bytes:
    """
    Return the file at url, or with member, the member of the zip of size bytes at url. Check its sha256. Fail if
    member comes without size.

    The file is read from the directory downloads if it is there, and otherwise downloaded and saved there. Its name in
    downloads is the last part of url, or of member. A member comes by HTTP range requests, so only its bytes are
    downloaded.
    """
    path = downloads / Path(member or urllib.parse.urlparse(url).path).name
    if path.exists():
        return check_sha256(path.read_bytes(), sha256, str(path))
    if member is None:
        with urllib.request.urlopen(url, timeout=120) as response:
            data = response.read()
    elif size is None:
        raise ValueError(f"the member {member} of {url} needs the size of the zip")
    else:
        # zipfile checks the CRC32 of the member
        with zipfile.ZipFile(RemoteFile(url, size)) as archive:
            data = archive.read(member)
    check_sha256(data, sha256, url)
    downloads.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return data


def read_transcripts(gff3: Path) -> pl.DataFrame:
    """
    Read the transcript lines of the GFF3, one row per transcript. The columns are transcript, strand and gene.
    transcript is the transcript ID without version. gene is the ID without version of the gene in the Parent of the
    line. A transcript line is any line whose ID attribute starts with ENST (exon, CDS and UTR lines have IDs like
    exon:ENST...).
    """
    transcripts_df = (
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
    if transcripts_df.is_empty():
        raise ValueError(f"{gff3} has no transcript line with an ID that starts with ENST")
    return transcripts_df


def hgvsc_alleles(table_df: pl.DataFrame, transcripts_df: pl.DataFrame, fasta: Fasta) -> pl.DataFrame:
    """
    Rename the columns chromosome, start and Transcript_ID of a table of Kim et al. 2024 to chrom, pos and transcript,
    and add ref and alt on the forward strand. transcripts_df is the table of read_transcripts.

    HGVSc (c.123G>T, optionally prefixed by the transcript ID) has the alleles on the transcript strand. On a
    minus-strand transcript of the GFF3, complement them. For a transcript missing from the GFF3, use the orientation
    whose REF matches the FASTA. Fail if HGVSc is not a single-base substitution. The transcript IDs of the table must
    lack the version, like those of read_transcripts.
    """
    table_df = table_df.rename({"chromosome": "chrom", "start": "pos", "Transcript_ID": "transcript"})
    table_df = table_df.with_columns(
        pl.col("HGVSc").str.extract_groups(HGVSC_SUBSTITUTION).struct.unnest(),
        genomic_ref=pl.Series(
            [fasta[chrom][pos - 1 : pos].seq.upper() for chrom, pos in zip(table_df["chrom"], table_df["pos"])],
            dtype=pl.String,
        ),
        strand=pl.col("transcript").replace_strict(
            transcripts_df["transcript"], transcripts_df["strand"], default=None
        ),
    )
    if (invalid_df := table_df.filter(pl.col("ref").is_null())).height > 0:
        row = invalid_df.row(0, named=True)
        raise ValueError(f"{row['HGVSc']} ({row['transcript']}) is not a single-base substitution")
    strand = pl.col("strand")
    complement = (strand == "-") | (strand.is_null() & (pl.col("ref") != pl.col("genomic_ref")))
    return table_df.with_columns(
        pl.when(complement).then(pl.col("ref", "alt").replace(COMPLEMENT)).otherwise(pl.col("ref", "alt"))
    ).drop("genomic_ref", "strand")


def finish_table(table_df: pl.DataFrame, transcripts_df: pl.DataFrame, fasta: Fasta) -> pl.DataFrame:
    """
    Add the gene of each transcript, check each REF against the FASTA, and keep the benchmark columns: COLUMNS,
    tissue if table_df has it, and source_row. gene comes from transcripts_df (see read_transcripts), and it is null
    for a transcript that the GFF3 lacks. Raise a ValueError at the first REF that differs from the FASTA.
    """
    table_df = table_df.with_columns(
        fasta_ref=pl.Series(
            [
                fasta[chrom][pos - 1 : pos - 1 + len(ref)].seq.upper()
                for chrom, pos, ref in table_df.select("chrom", "pos", "ref").iter_rows()
            ],
            dtype=pl.String,
        ),
        gene=pl.col("transcript").replace_strict(transcripts_df["transcript"], transcripts_df["gene"], default=None),
    )
    if (mismatches_df := table_df.filter(pl.col("ref") != pl.col("fasta_ref"))).height > 0:
        row = mismatches_df.row(0, named=True)
        raise ValueError(
            f"{row['chrom']}:{row['pos']} ({row['transcript']}): REF {row['ref']} but the FASTA has {row['fasta_ref']}"
        )
    return table_df.select(*COLUMNS, *(["tissue"] if "tissue" in table_df.columns else []), "source_row")
