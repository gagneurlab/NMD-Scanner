# Import dependencies
import gzip
import itertools
import logging
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import polars as pl
from Bio.Seq import Seq
from pyfaidx import Fasta

from nmd_scanner._polars_bio import pb

# Create the functions used for reading in the files (VCF, GFF3, FASTA)

GFF3_SUFFIXES = (".gff3", ".gff3.gz", ".gff", ".gff.gz")

# The polars-bio errors for a VCF without its header
VCF_HEADER_ERRORS = ("empty input", "invalid record", "missing header")

STOP_CODONS = {"TAA", "TAG", "TGA"}
MITOCHONDRIAL_CHROMOSOMES = {"MT", "M", "chrM", "chrMT"}
# vertebrate mitochondrial code (NCBI translation table 2): AGA and AGG are stop codons, TGA codes for Trp
MITOCHONDRIAL_STOP_CODONS = {"TAA", "TAG", "AGA", "AGG"}

logger = logging.getLogger(__name__)


def read_vcf(vcf_path):
    # TODO: adjust this function to also get structural variants and then also adjust downstream analysis
    #  (especially the inclusion of the Variants into the reference CDS sequence to create the alternative CDS).
    #  Until then, extract_ptc skips symbolic alleles and breakends.

    """
    Reads a VCF file, plain or gzip-compressed, with polars-bio into a pandas DataFrame with the
    columns Chromosome, Start, End, ID, Ref and Alt, one row per record in file order, and a
    RangeIndex. QUAL, FILTER and INFO are not read.

    Start and End are 0-based half-open: Start is POS - 1, and End is Start plus the length of REF,
    also for a symbolic allele with an INFO END. Every other field stays text as written, e.g.
    ``01``, ``007`` and ``NA``. polars-bio returns ``.`` in ID and ALT as an empty string, which
    becomes ``.`` again: a VCF field is never empty.

    The VCF needs its header, at least the ``##fileformat`` and ``#CHROM`` lines. It must be
    left-normalized and single-allelic (one ALT allele per record), e.g. produced by
    ``bcftools norm -m- -f reference.fa``. Multi-allelic records are rejected because the
    downstream variant application assumes exactly one ALT allele per row.

    polars-bio rejects some files that VCF 4.3 allows, and so does read_vcf: a file with an empty
    last line, a file that starts with a UTF-8 byte order mark, and a record with POS 0.

    :raises FileNotFoundError: if ``vcf_path`` does not exist
    :raises OSError: if ``vcf_path`` cannot be opened, e.g. because it is a directory
    :raises ValueError: if polars-bio cannot read the file, e.g. because it has no header, or if
        the VCF has multi-allelic records
    """
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    # open raises the OSError that fits, e.g. IsADirectoryError; polars-bio's error does not say which
    with open(vcf_path, "rb"):
        pass
    try:
        # info_fields=[]: no INFO field is parsed, nothing reads them
        df = (
            pb.scan_vcf(os.fspath(vcf_path), info_fields=[], use_zero_based=True)
            .select(
                pl.col("chrom").alias("Chromosome"),
                pl.col("start").cast(pl.Int64).alias("Start"),
                # not the end polars-bio gives, which follows INFO END for a symbolic allele
                (pl.col("start").cast(pl.Int64) + pl.col("ref").str.len_chars().cast(pl.Int64)).alias("End"),
                pl.col("id").replace("", ".").alias("ID"),
                pl.col("ref").alias("Ref"),
                pl.col("alt").replace("", ".").alias("Alt"),
            )
            .collect()
            .to_pandas()
        )
    # polars < 2 wraps a polars-bio read error in ValueError or ComputeError. polars 2 passes on the
    # bare Exception that DataFusion raises, so no narrower type exists to catch.
    except Exception as error:
        message = f"Cannot read {os.fspath(vcf_path)!r} as VCF ({error})."
        # polars-bio reads the header when it registers the file as a table
        if "Failed to register table" in str(error) and any(text in str(error) for text in VCF_HEADER_ERRORS):
            message += " The VCF needs its header, at least the ##fileformat and #CHROM lines."
        raise ValueError(message) from error

    # Reject multi-allelic records: the VCF spec allows several ALT alleles per record, but the
    # rest of the pipeline assumes one. polars-bio joins the ALT alleles of a record with "|", which
    # its documentation does not say. A symbolic allele (<...>) and the mate of a breakend ([...[ or
    # ]...]) may contain "|" themselves, so their text is removed first.
    alleles = df["Alt"].str.replace(r"<[^>]*>", "<>", regex=True).str.replace(r"[\[\]][^\[\]]*[\[\]]", "[]", regex=True)
    multiallelic = alleles.str.contains("|", regex=False)
    if multiallelic.any():
        n_multiallelic = int(multiallelic.sum())
        raise ValueError(
            f"{n_multiallelic} multi-allelic record(s) found in {vcf_path}. "
            "NMD-Scanner requires a left-normalized, single-allelic VCF. "
            "Split and normalize first, e.g. `bcftools norm -m- -f reference.fa`."
        )

    return df


def detect_annotation_format(path):
    """
    Detects the format of the annotation file ``path`` (str or path-like) from its filename suffix,
    gzip-compressed or not.

    :return: "gff3" for a ``.gff3`` or ``.gff`` file
    :raises ValueError: for any other suffix
    """
    path = os.fspath(path)
    if path.lower().endswith(GFF3_SUFFIXES):
        return "gff3"
    raise ValueError(
        f"Cannot detect annotation format from filename: {path!r}. "
        f"Expected one of the suffixes {', '.join(GFF3_SUFFIXES)}."
    )


def read_annotation(path, fasta, reassign_exons=False):
    """
    Reads a GFF3 gene annotation file into a pandas DataFrame with the exon rows and the coding
    regions that ``rules.extract_ptc`` takes (see ``read_gff3``).

    The coding regions are the CDS rows: they include the stop codon, one row per transcript and
    exon. On them, the column has_start_codon says whether the coding region of the transcript
    starts with an annotated start codon, and has_stop_codon whether it ends in an annotated stop
    codon. On the exon rows, both are NA.

    The format is checked from the filename suffix (``.gff3`` or ``.gff``, gzip-compressed or not).

    :param path: Path to the GFF3 file
    :param fasta: Reference genome (pyfaidx.Fasta object), which shows whether a CDS ends in a stop
        codon, and for an Ensembl GFF3 whether it starts with one
    :param reassign_exons: Recompute the exon numbers with ``compute_exon_numbers``
    :return: DataFrame, as ``read_gff3`` returns it
    :raises ValueError: for an unknown suffix
    """
    detect_annotation_format(path)
    annotation = read_gff3(path, fasta)
    if reassign_exons:
        logger.info("Recomputing exon numbers.")
        annotation = compute_exon_numbers(annotation)
    return annotation


def read_gff3(gff3_path, fasta):
    """
    Reads a GFF3 file into a pandas DataFrame with the exon rows and the coding regions that
    ``rules.extract_ptc`` takes.

    Two GFF3 flavors are supported. The flavor is detected from the attributes of the rows that
    ``_read_gff3_rows`` reads, so gene_type and transcript_type (GENCODE) or biotype (Ensembl) must
    be on the transcript rows, not only on the gene rows:

    - GENCODE: every feature row carries ``gene_id``/``transcript_id`` directly, and they are kept
      (e.g. hg19 lift37 has ``ENST00000477874.1_2`` there but ``ENST00000477874.1`` in ``ID``).
      Only chrY PAR transcripts reuse the chrX value on that attribute (the "_PAR_Y" suffix only
      shows up in ``ID``/``Parent``), so for those the ``ID``/``Parent`` hierarchy gives
      ``gene_id``/``transcript_id``.
    - Ensembl: ``gene_id``/``transcript_id`` are never given directly on exon/CDS rows and are
      resolved from the ``ID``/``Parent`` hierarchy (``gene:``/``transcript:`` prefixed). Exon rows
      carry their number as ``rank``; a CDS row takes the number of the exon it lies in. Without
      ``rank``, ``compute_exon_numbers`` computes the numbers.

    The coding regions are the CDS rows. A GFF3 CDS includes the start and the stop codon, so the
    result has no start_codon or stop_codon rows. On the CDS rows, the column has_start_codon says
    whether the coding region of the transcript starts with an annotated start codon, and
    has_stop_codon whether it ends in an annotated stop codon. On the exon rows, both are NA. The two
    flavors get the coding regions and the two columns as follows:

    - GENCODE: has_start_codon is True if the transcript has start_codon rows. A start codon can be
      a non-ATG codon such as CTG. A cds_start_NF transcript has none. has_stop_codon is True if the
      transcript has stop_codon rows. A few cds_end_NF transcripts without stop_codon rows have a
      CDS that ends in a complete stop codon; ``_trim_cds_end_nf_stop_codons`` removes it from the
      CDS.
    - Ensembl: the GFF3 has no start_codon or stop_codon rows. has_start_codon comes from the first
      3 CDS bases in ``fasta`` (``_has_start_codon_from_sequence``), and has_stop_codon from the
      last 3 (``_has_stop_codon_from_sequence``).

    The file is read with polars-bio (see ``_read_gff3_rows``), plain or gzip-compressed. Only the
    attributes in ``GFF3_ATTRIBUTES`` are read.

    :param gff3_path: Path to the GFF3 file
    :param fasta: Reference genome (pyfaidx.Fasta object)
    :return: DataFrame with the exon and CDS rows in file order and a RangeIndex, the columns
        Chromosome, Source, Feature, Start, End, Score, Strand, Frame, has_start_codon,
        has_stop_codon, gene_id and transcript_id, and the other attributes of ``GFF3_ATTRIBUTES``
        that the file has, except ID and Parent. exon_number (nullable integer) is there for an
        Ensembl GFF3, and for a GENCODE GFF3 that has the attribute.
    :raises ValueError: if the file has no ID or Parent attribute, e.g. because it is a GTF, if its
        flavor is neither GENCODE nor Ensembl, or if an exon or CDS row breaks a rule of
        ``_check_strand_and_phase`` or ``_check_cds_in_exons``
    """
    if not os.path.exists(gff3_path):
        raise FileNotFoundError(f"GFF3 file not found: {gff3_path}")
    # open raises the OSError that fits, e.g. IsADirectoryError; polars-bio's error does not say which
    with open(gff3_path, "rb"):
        pass

    df = _read_gff3_rows(gff3_path)
    columns = set(df.columns)
    if not {"ID", "Parent"} & columns:
        # polars-bio reads a GTF as GFF3 rows without attributes
        raise ValueError(
            f"Cannot read {os.fspath(gff3_path)!r} as GFF3: no row has an ID or Parent attribute, "
            "e.g. because it is a GTF."
        )

    if {"gene_type", "transcript_type"} <= columns:
        df = _normalize_gencode_gff3(df)
        _check_strand_and_phase(df, gff3_path)
        _check_cds_in_exons(df, gff3_path)
        df = _trim_cds_end_nf_stop_codons(df, fasta)
        # the CDS includes the stop codon; the stop_codon rows only say whether there is one
        is_stop = df["Feature"] == "stop_codon"
        if not is_stop.any() and (df["Feature"] == "CDS").any():
            logger.warning(
                "No stop_codon rows found next to the CDS rows: every transcript is treated as having no annotated "
                "stop codon (no 3'UTR length, no stop codon distance, every in-frame stop is premature). "
                "A GENCODE GFF3 marks the stop codons with stop_codon rows."
            )
        # the CDS starts with the start codon; the start_codon rows only say whether there is one
        is_start = df["Feature"] == "start_codon"
        if not is_start.any() and (df["Feature"] == "CDS").any():
            logger.warning(
                "No start_codon rows found next to the CDS rows: every transcript is treated as having no annotated "
                "start codon (no start loss, no start codon position, no PTC distance to the start codon). "
                "A GENCODE GFF3 marks the start codons with start_codon rows."
            )
        stop_transcripts = df.loc[is_stop, "transcript_id"]
        start_transcripts = df.loc[is_start, "transcript_id"]
        df = _set_cds_flag(df[~(is_stop | is_start)], "has_stop_codon", stop_transcripts)
        df = _set_cds_flag(df, "has_start_codon", start_transcripts)
    elif "biotype" in columns:
        df = _normalize_ensembl_gff3(df)
        _check_strand_and_phase(df, gff3_path)
        _check_cds_in_exons(df, gff3_path)
        df = _has_stop_codon_from_sequence(df, fasta)
        df = _has_start_codon_from_sequence(df, fasta)
    else:
        raise ValueError(
            "Unrecognized GFF3 flavor: expected GENCODE-style attributes "
            "(gene_type/transcript_type) or Ensembl-style attributes (biotype)."
        )

    # ID/Parent are GFF3 structural columns; _resolve_ids_via_hierarchy has used them by this
    # point. Drop them so a stray "ID" column doesn't collide with the unrelated "ID" column
    # read_vcf uses for the VCF record ID once CDS rows are joined to it.
    df = df.drop(columns=[c for c in ("ID", "Parent") if c in df.columns])

    return df.reset_index(drop=True)


def _check_strand_and_phase(df, gff3_path):
    """
    Raises a ValueError if an exon or CDS row has a strand other than + or -, or if a CDS row has a
    phase (Frame) other than 0, 1 or 2. The analysis reads the strand of each transcript and the
    phase of its 5'-most CDS row, and GFF3 requires a phase on each CDS row. The error gives the
    count of such rows and names the transcript of the first one, in file order.

    :param df: Rows of the GFF3 (DataFrame) with Chromosome, Feature, Start, End, Strand, Frame and
        transcript_id, as the normalization of its flavor gives them
    :param gff3_path: Path to the GFF3 file, for the error message
    :raises ValueError: for such a row
    """
    path = os.fspath(gff3_path)
    rows = df[df["Feature"].isin(["exon", "CDS"])]
    unstranded = rows[~rows["Strand"].isin(["+", "-"])]
    if not unstranded.empty:
        first = unstranded.iloc[0]
        raise ValueError(
            f"Cannot use {path!r}: {len(unstranded)} exon or CDS row(s) have a strand other than + or -. The first "
            f"is {_where(first)}, with strand {first['Strand']!r}. NMD-Scanner needs strand + or - on each exon and "
            "CDS row."
        )
    cds = rows[rows["Feature"] == "CDS"]
    unphased = cds[~cds["Frame"].isin(["0", "1", "2"])]
    if not unphased.empty:
        first = unphased.iloc[0]
        raise ValueError(
            f"Cannot use {path!r}: {len(unphased)} CDS row(s) have a phase other than 0, 1 or 2. The first is "
            f"{_where(first)}, with phase {first['Frame']!r}. GFF3 requires a phase on each CDS row."
        )


def _check_cds_in_exons(df, gff3_path):
    """
    Raises a ValueError if a CDS row does not lie inside an exon row of its transcript, i.e. if no exon row of the
    transcript, on the same chromosome, starts at or before the CDS row and ends at or after it. The analysis reads
    the CDS bases from the exons of the transcript. CDS rows that share bases are allowed, e.g. the two CDS rows at a
    ribosomal slippage site. The check skips a transcript without exon rows: the output columns have a null case for it.
    The error gives the count of such rows and names the first one, in file order.

    :param df: Rows of the GFF3 (DataFrame) with Chromosome, Feature, Start, End and transcript_id, as the
        normalization of its flavor gives them
    :param gff3_path: Path to the GFF3 file, for the error message
    :raises ValueError: for such a row
    """
    key = ["Chromosome", "transcript_id"]
    exons = df.loc[df["Feature"] == "exon", [*key, "Start", "End"]].sort_values([*key, "Start"])
    # The largest End of the exon rows of the transcript that start at or before each exon row
    exons["reach"] = exons.groupby(key, sort=False)["End"].cummax()
    cds = df.loc[df["Feature"] == "CDS", [*key, "Start", "End"]]
    cds = cds[cds["transcript_id"].isin(exons["transcript_id"])]
    # Each CDS row takes the reach of the last exon row of its transcript that starts at or before it
    reach = pd.merge_asof(
        cds.assign(row=cds.index).sort_values("Start"),
        exons[[*key, "Start", "reach"]].sort_values("Start", kind="stable"),
        on="Start",
        by=key,
    ).set_index("row")["reach"]
    outside = cds.index[~(reach.reindex(cds.index) >= cds["End"])]
    if len(outside):
        first = df.loc[outside[0]]
        raise ValueError(
            f"Cannot use {os.fspath(gff3_path)!r}: {len(outside)} CDS row(s) do not lie inside an exon row of their "
            f"transcript. The first is {_where(first)}. NMD-Scanner needs each CDS row inside an exon row of its "
            "transcript."
        )


def _where(row):
    """The feature, transcript and location (1-based, as in the GFF3) of a row, for an error message."""
    return (
        f"the {row['Feature']} row of transcript {row['transcript_id']} at "
        f"{row['Chromosome']}:{row['Start'] + 1}-{row['End']}"
    )


# The GFF3 attributes that read_gff3 and its helpers use. polars-bio reads no other attribute.
GFF3_ATTRIBUTES = (
    "ID",
    "Parent",
    "gene_id",
    "transcript_id",
    "gene_type",
    "transcript_type",
    "biotype",
    "tag",
    "exon_number",
    "rank",
    "ensembl_end_phase",
)
# The feature types of the rows that read_gff3 uses. Ensembl has many transcript types (mRNA,
# lnc_RNA, ...), so any row whose Parent is a gene is read too.
GFF3_FEATURES = ("transcript", "exon", "CDS", "start_codon", "stop_codon")


def _read_gff3_rows(gff3_path):
    """
    Reads the rows of a GFF3 file that ``read_gff3`` uses, with polars-bio: the rows of a type in
    ``GFF3_FEATURES`` and the rows whose Parent is a gene.

    :param gff3_path: Path to the GFF3 file, plain or gzip-compressed
    :return: DataFrame with the columns Chromosome, Source, Feature (category), Start, End, Score
        (float64, NaN if missing), Strand and Frame (text, "." if missing), and the attributes in
        ``GFF3_ATTRIBUTES`` that have a value in at least one row, as text (NaN where a row has
        none). Start and End are 0-based half-open int64. polars-bio reads Score as float32, so it
        keeps about 7 significant digits. In attribute values, polars-bio percent-decodes only %3B,
        %3D, %26, %2C and %09, with upper-case hex digits; it leaves every other escape, e.g. %25 or
        %3b, as it is. A decoded %2C is a comma, so ``tag=a%2Ccds_end_NF`` reads as the two tags a and
        cds_end_NF.
    :raises ValueError: if polars-bio cannot read the file, if it skips a data line of the file, or if
        a row has a start after its end
    """
    try:
        with ThreadPoolExecutor(max_workers=1) as executor:
            # the data lines are counted while polars-bio reads the file
            data_lines = executor.submit(_count_gff3_data_lines, gff3_path)
            # All rows are collected and filtered afterwards: polars-bio filters this predicate in
            # Python, which takes longer, and the row count shows the lines that polars-bio skipped.
            rows = (
                pb.scan_gff(os.fspath(gff3_path), attr_fields=list(GFF3_ATTRIBUTES), use_zero_based=True)
                # polars-bio's GFF frame takes the columns of select as one list
                .select(
                    [
                        pl.col("chrom").alias("Chromosome"),
                        pl.col("source").alias("Source"),
                        pl.col("type").cast(pl.Categorical).alias("Feature"),
                        pl.col("start").cast(pl.Int64).alias("Start"),
                        pl.col("end").cast(pl.Int64).alias("End"),
                        pl.col("score").cast(pl.Float64).alias("Score"),
                        pl.col("strand").alias("Strand"),
                        pl.col("phase").cast(pl.String).fill_null(".").alias("Frame"),
                        *GFF3_ATTRIBUTES,
                    ]
                )
                .collect()
            )
    # polars < 2 wraps a polars-bio read error in ValueError or ComputeError. polars 2 passes on the
    # bare Exception that DataFusion raises, so no narrower type exists to catch.
    except Exception as error:
        raise ValueError(f"Cannot read {os.fspath(gff3_path)!r} as GFF3 ({error}).") from error
    # polars-bio skips a line that it cannot parse, and logs it at DEBUG only
    if rows.height < data_lines.result():
        raise ValueError(
            f"Cannot read {os.fspath(gff3_path)!r} as GFF3: polars-bio read {rows.height} of its "
            f"{data_lines.result()} data lines. It skips a line that it cannot parse, e.g. a line whose Start or "
            "Score is not a number or whose columns are not separated by tabs."
        )
    rows = rows.filter(pl.col("Feature").is_in(GFF3_FEATURES) | pl.col("Parent").str.starts_with("gene:"))
    # polars-bio gives Start as a 32-bit integer without sign, so a start of 0 wraps to the largest value
    reversed_rows = rows.filter(pl.col("Start") >= pl.col("End"))
    if reversed_rows.height:
        first = reversed_rows.row(0, named=True)
        raise ValueError(
            f"Cannot read {os.fspath(gff3_path)!r} as GFF3: {reversed_rows.height} row(s) have a start after their "
            f"end, e.g. the {first['Feature']} row at {first['Chromosome']}:{first['Start'] + 1}-{first['End']}. "
            "polars-bio reads a start of 0 as 4294967296."
        )
    # polars-bio gives an attribute that no row has as a column without values
    absent = [name for name in GFF3_ATTRIBUTES if rows[name].null_count() == rows.height]
    return rows.drop(absent).to_pandas()


def _count_gff3_data_lines(gff3_path):
    """
    Counts the data lines of a GFF3 file, plain or gzip-compressed: the lines before a ##FASTA
    directive that are not blank and do not start with "#". polars-bio skips blank lines too.
    """
    count = 0
    with open(gff3_path, "rb") as file:
        gzipped = file.read(2) == b"\x1f\x8b"
        file.seek(0)
        for line in gzip.GzipFile(fileobj=file) if gzipped else file:
            if line.startswith(b"##FASTA"):
                break
            if not line.startswith(b"#") and not line.isspace():
                count += 1
    return count


def _has_stop_codon_from_sequence(df, fasta):
    """
    Sets has_stop_codon on the CDS rows of an Ensembl GFF3 from the last 3 CDS bases in ``fasta``.

    A transcript has a stop codon if the last 3 bases of its CDS are a stop codon (see
    ``_last_codons``), unless its last coding exon ends mid-codon: the exon attribute
    ``ensembl_end_phase`` is 1 or 2 there. Ensembl annotates no stop codon for a transcript tagged
    cds_end_NF either, but the GFF3 has no cds_end_NF tag. So a cds_end_NF transcript whose CDS ends
    in stop codon bases gets has_stop_codon True here: 13 transcripts in Ensembl 108, none on chr22.

    :param df: Exon and CDS rows of an Ensembl GFF3 (DataFrame) with transcript_id and exon_number
    :param fasta: Reference genome (pyfaidx.Fasta object)
    :return: DataFrame with the column has_stop_codon (see ``_set_cds_flag``)
    """
    df = df.reset_index(drop=True)
    codons = _last_codons(df[df["Feature"] == "CDS"], fasta)

    mid_codon: set[tuple] = set()
    if "ensembl_end_phase" in df.columns:
        exons = df[(df["Feature"] == "exon") & pd.to_numeric(df["ensembl_end_phase"], errors="coerce").isin([1, 2])]
        mid_codon = set(zip(exons["transcript_id"], exons["exon_number"]))
    last_exon = list(zip(codons["transcript_id"], df.loc[codons["last_cds_index"], "exon_number"]))
    ends_mid_codon = pd.Series([key in mid_codon for key in last_exon], index=codons.index, dtype=bool)

    stops = codons[codons["is_stop"] & ~ends_mid_codon]
    logger.info(
        "Stop codons from the FASTA: %d transcripts end in a stop codon; %d more end in stop codon bases, "
        "but their last coding exon ends mid-codon (ensembl_end_phase 1 or 2).",
        len(stops),
        int((codons["is_stop"] & ends_mid_codon).sum()),
    )
    logger.warning(
        "Ensembl GFF3 has no cds_end_NF tag: a cds_end_NF transcript whose CDS ends in stop codon bases "
        "gets a stop codon (13 transcripts in Ensembl 108)."
    )
    return _set_cds_flag(df, "has_stop_codon", stops["transcript_id"])


def _has_start_codon_from_sequence(df, fasta):
    """
    Sets has_start_codon on the CDS rows of an Ensembl GFF3 from the first 3 CDS bases in ``fasta``.

    An Ensembl GFF3 has no start_codon rows. So a transcript has a start codon here if its CDS starts
    with ATG and its 5'-most CDS row has phase 0 (see ``_first_codons``). A non-ATG start codon,
    e.g. CTG, is not found. A transcript without the 5' end of its CDS (cds_start_NF) whose CDS
    starts with ATG in phase 0 gets a start codon.

    :param df: Exon and CDS rows of an Ensembl GFF3 (DataFrame) with transcript_id
    :param fasta: Reference genome (pyfaidx.Fasta object)
    :return: DataFrame with the column has_start_codon (see ``_set_cds_flag``)
    """
    codons = _first_codons(df[df["Feature"] == "CDS"], fasta)
    starts = codons[(codons["codon"] == "ATG") & (codons["phase"] == 0)]
    logger.info(
        "Start codons from the FASTA: %d of %d transcripts start with ATG in phase 0.", len(starts), len(codons)
    )
    logger.warning(
        "Ensembl GFF3 has no start_codon rows: a CDS has a start codon if it starts with ATG in phase 0. "
        "A non-ATG start codon, e.g. CTG, is not found."
    )
    return _set_cds_flag(df, "has_start_codon", starts["transcript_id"])


def _set_cds_flag(df, column, transcripts):
    """
    Returns a copy of ``df`` with the column ``column`` (nullable boolean): on the CDS rows, whether
    their transcript is in ``transcripts``; NA on the other rows.
    """
    flag = df["transcript_id"].isin(transcripts).astype("boolean")
    return df.assign(**{column: flag.where(df["Feature"] == "CDS")})


def _trim_cds_end_nf_stop_codons(df, fasta):
    """
    Removes the last codon from the CDS of a GENCODE GFF3 for transcripts tagged cds_end_NF,
    without stop_codon rows, whose CDS ends in a complete stop codon (see ``_last_codons``).

    cds_end_NF means that the end of the CDS is not confirmed, so these 3 bases are no annotated
    stop codon. The GENCODE GTF of the same release has them as UTR.

    This concerns 10 transcripts in GENCODE 42, none on chr22. GENCODE GFF3 repeats the transcript
    tags on every row, so the tags of the CDS rows are used.

    :param df: Exon, CDS and stop_codon rows of a GENCODE GFF3 (DataFrame) with transcript_id and tag
    :param fasta: Reference genome (pyfaidx.Fasta object)
    :return: DataFrame with the CDS rows shortened
    """
    df = df.reset_index(drop=True)
    if "tag" not in df.columns:
        return df
    cds = df[df["Feature"] == "CDS"]
    cds_end_nf = cds["tag"].fillna("").str.split(",").map(lambda tags: "cds_end_NF" in tags)
    stop_transcripts = df.loc[df["Feature"] == "stop_codon", "transcript_id"]
    codons = _last_codons(cds[cds_end_nf & ~cds["transcript_id"].isin(stop_transcripts)], fasta)
    stops = codons[codons["is_stop"] & codons["in_frame"]]

    df = df.copy()
    for cds_index, start, end in stops.explode("pieces")["pieces"]:
        # the stop codon piece is at the 3' end of its CDS row
        if df.at[cds_index, "Strand"] == "+":
            df.at[cds_index, "End"] = start
        else:
            df.at[cds_index, "Start"] = end
    logger.info(
        "Stop codons: %d cds_end_NF transcripts without stop_codon rows end in a stop codon; "
        "it is removed from the CDS.",
        len(stops),
    )
    # a CDS row that held only stop codon bases is gone
    return df[df["Start"] < df["End"]]


def _first_codons(cds, fasta):
    """
    Reads the first codon of each transcript's CDS: its first 3 bases, read in the direction of
    transcription. Transcripts with fewer than 3 CDS bases or on a chromosome that is not in
    ``fasta`` are left out.

    :param cds: CDS rows (DataFrame) with Chromosome, Start, End, Strand, Frame and transcript_id
    :param fasta: Reference genome (pyfaidx.Fasta object)
    :return: DataFrame with one row per transcript and the columns transcript_id, codon (upper case)
        and phase (the Frame of the CDS row at the 5' end, NaN if it is not a number)
    """
    plus = cds["Strand"] == "+"
    # each transcript's CDS rows from its 5' end: by Start ascending on the plus strand, by End
    # descending on the minus strand
    cds = cds.assign(_from_5_prime=cds["Start"].where(plus, -cds["End"])).sort_values(
        ["transcript_id", "_from_5_prime"]
    )
    rows = zip(
        cds["transcript_id"],
        cds["Chromosome"],
        cds["Strand"],
        cds["Start"],
        cds["End"],
        pd.to_numeric(cds["Frame"], errors="coerce"),
    )

    codons = []
    for transcript_id, transcript_rows in itertools.groupby(rows, key=lambda row: row[0]):
        transcript_rows = list(transcript_rows)
        _, chrom, strand, _, _, phase = transcript_rows[0]
        if chrom not in fasta:
            continue
        # the first 3 bases, from the 5' end
        codon = ""
        missing = 3
        for _, _, _, start, end, _ in transcript_rows:
            take = min(missing, end - start)
            if strand == "+":
                codon += str(fasta[chrom][start : start + take]).upper()
            else:
                codon += str(Seq(str(fasta[chrom][end - take : end]).upper()).reverse_complement())
            missing -= take
            if missing == 0:
                break
        if not missing:
            codons.append((transcript_id, codon, phase))
    return pd.DataFrame(codons, columns=["transcript_id", "codon", "phase"])


def _last_codons(cds, fasta):
    """
    Reads the last codon of each transcript's CDS: its last 3 bases, read in the direction of
    transcription. Transcripts with fewer than 3 CDS bases or on a chromosome that is not in
    ``fasta`` are left out.

    :param cds: CDS rows (DataFrame) with Chromosome, Start, End, Strand, Frame and transcript_id
    :param fasta: Reference genome (pyfaidx.Fasta object)
    :return: DataFrame with one row per transcript and the columns

        - transcript_id
        - last_cds_index: index of the CDS row at the 3' end
        - pieces: the codon bases as (cds_index, Start, End) tuples, from 5' to 3'. A codon split
          across an intron has 2 pieces.
        - is_stop: the codon is TAA, TAG or TGA; on the mitochondrial chromosome TAA, TAG, AGA or AGG
        - in_frame: the CDS row at the 3' end, after its phase (Frame), holds a multiple of 3 bases
    """
    plus = cds["Strand"] == "+"
    # each transcript's CDS rows from its 3' end: by End descending on the plus strand, by Start ascending on the minus strand
    cds = cds.assign(_from_3_prime=(-cds["End"]).where(plus, cds["Start"])).sort_values(
        ["transcript_id", "_from_3_prime"]
    )
    rows = zip(
        cds.index,
        cds["transcript_id"],
        cds["Chromosome"],
        cds["Strand"],
        cds["Start"],
        cds["End"],
        pd.to_numeric(cds["Frame"], errors="coerce"),
    )

    codons = []
    skipped = 0
    for transcript_id, transcript_rows in itertools.groupby(rows, key=lambda row: row[1]):
        transcript_rows = list(transcript_rows)
        last_cds_index, _, chrom, strand, last_start, last_end, phase = transcript_rows[0]

        # the last 3 bases, from the 3' end: (cds_index, Start, End)
        parts = []
        missing = 3
        for cds_index, _, _, _, start, end, _ in transcript_rows:
            take = min(missing, end - start)
            parts.append((cds_index, end - take, end) if strand == "+" else (cds_index, start, start + take))
            missing -= take
            if missing == 0:
                break
        if missing:
            continue
        if chrom not in fasta:
            skipped += 1
            continue

        codon = "".join(str(fasta[chrom][s:e]).upper() for _, s, e in sorted(parts, key=lambda part: part[1]))
        if strand == "-":
            codon = str(Seq(codon).reverse_complement())
        codons.append(
            (
                transcript_id,
                last_cds_index,
                parts[::-1],
                codon in (MITOCHONDRIAL_STOP_CODONS if chrom in MITOCHONDRIAL_CHROMOSOMES else STOP_CODONS),
                not pd.isna(phase) and (last_end - last_start - int(phase)) % 3 == 0,
            )
        )

    if skipped:
        logger.warning("No stop codon check for %d transcripts: their chromosome is not in the FASTA.", skipped)
    return pd.DataFrame(codons, columns=["transcript_id", "last_cds_index", "pieces", "is_stop", "in_frame"])


def _resolve_ids_via_hierarchy(df, transcript_mask, keep_attributes=False):
    """
    Resolve gene_id/transcript_id for every row of a GFF3 dataframe from the ID/Parent
    hierarchy rather than from the flat gene_id/transcript_id attributes.

    This matters because some GFF3 flavors do not give exon/CDS rows a gene_id/transcript_id
    attribute at all (Ensembl), and others reuse the same gene_id/transcript_id attribute value
    for genomically distinct copies of a feature (GENCODE's chrY PAR transcripts: the "_PAR_Y"
    suffix that disambiguates them from their chrX counterpart only shows up in ID/Parent, not
    in the gene_id/transcript_id attributes).

    With ``keep_attributes``, the gene_id/transcript_id attributes stay where a row has them, and
    the hierarchy only fills in the rows without one and the "_PAR_Y" rows. The ID can differ from
    the attribute (GENCODE lift37: ``ENST00000477874.1_2`` vs ``ENST00000477874.1``).

    :param df: Attribute-parsed GFF3 dataframe; must have ``ID`` and ``Parent`` columns.
    :param transcript_mask: Boolean mask selecting the rows that define a transcript (their
        ``ID`` is the transcript identifier, their ``Parent`` is the gene identifier).
    :return: A copy of ``df`` with (overwritten) ``gene_id`` and ``transcript_id`` columns, set
        for the transcript rows themselves and for their children (rows whose ``Parent`` points
        at a transcript row's ``ID``); unrelated rows (e.g. ``gene`` rows) are left as NaN.
    """

    def strip_prefix(series):
        stripped = series.astype(str).str.replace(r"^(?:gene|transcript):", "", regex=True)
        return stripped.where(series.notna(), other=pd.NA)

    df = df.copy()
    ids = strip_prefix(df["ID"])
    parents = strip_prefix(df["Parent"])

    tx_rows = transcript_mask & ids.notna() & parents.notna()
    tx_to_gene = pd.Series(parents[tx_rows].to_numpy(), index=ids[tx_rows].to_numpy())

    transcript_id = pd.Series(pd.NA, index=df.index, dtype="object")
    gene_id = pd.Series(pd.NA, index=df.index, dtype="object")

    transcript_id[tx_rows] = ids[tx_rows]
    gene_id[tx_rows] = parents[tx_rows]

    is_child = parents.isin(tx_to_gene.index)
    transcript_id[is_child] = parents[is_child]
    gene_id[is_child] = parents[is_child].map(tx_to_gene)

    if keep_attributes:
        for name, resolved in (("transcript_id", transcript_id), ("gene_id", gene_id)):
            if name in df.columns:
                use_hierarchy = df[name].isna() | resolved.astype(str).str.endswith("_PAR_Y")
                resolved = resolved.where(use_hierarchy, df[name])
            df[name] = resolved
        return df

    df["transcript_id"] = transcript_id
    df["gene_id"] = gene_id
    return df


def _normalize_gencode_gff3(df):
    """
    Returns the exon, CDS, start_codon and stop_codon rows of a GENCODE-flavored GFF3 dataframe
    (see ``read_gff3``), with gene_id/transcript_id from the attributes
    (from ID/Parent where missing and for chrY PAR copies).
    """
    transcript_mask = df["Feature"] == "transcript"
    df = _resolve_ids_via_hierarchy(df, transcript_mask, keep_attributes=True)

    child_mask = df["Feature"].isin(["exon", "CDS", "start_codon", "stop_codon"])
    _check_ids_resolved(df, child_mask)

    child_df = df.loc[child_mask].reset_index(drop=True)
    # the attribute is text; an Ensembl GFF3 gives the numbers as nullable integer too
    if "exon_number" in child_df.columns:
        child_df["exon_number"] = child_df["exon_number"].astype("Int64")
    return child_df


def _normalize_ensembl_gff3(df):
    """
    Returns the exon and CDS rows of an Ensembl-flavored GFF3 dataframe (see ``read_gff3``),
    with gene_id/transcript_id resolved from ID/Parent and exon numbers.
    """
    transcript_mask = df["ID"].notna() & df["Parent"].astype(str).str.startswith("gene:")
    df = _resolve_ids_via_hierarchy(df, transcript_mask)

    child_mask = df["Feature"].isin(["exon", "CDS"])
    _check_ids_resolved(df, child_mask)

    child_df = df.loc[child_mask].reset_index(drop=True)

    # Exon rows carry their number as the "rank" attribute, CDS rows carry nothing. Without
    # "rank", compute the numbers from genomic order.
    if "rank" in child_df.columns:
        rank = pd.to_numeric(child_df["rank"], errors="coerce")
        if rank[child_df["Feature"] == "exon"].notna().all():
            child_df["exon_number"] = rank.where(child_df["Feature"] == "exon").astype("Int64")
            return _assign_exon_numbers_to_cds(child_df)
    return compute_exon_numbers(child_df)


def _check_ids_resolved(df, mask):
    unresolved = int(df.loc[mask, "transcript_id"].isna().sum())
    if unresolved:
        raise ValueError(f"Could not resolve transcript_id via the ID/Parent hierarchy for {unresolved} row(s).")


def read_fasta(fasta_path):
    """
    Reads a genome FASTA file using pyfaidx.Fasta and returns a pyfaidx.Fasta object.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    return Fasta(fasta_path)


def compute_exon_numbers(annotation):
    """
    Compute exon numbers for the Features exon, CDS and stop_codon in an annotation DataFrame.
    Exon numbers are assigned based on genomic order per transcript and strand.
    CDS and stop_codon features inherit the exon number of the exon they overlap.

    The exon numbers follow the direction of transcription: exon 1 has the smallest Start on
    the + strand and the largest Start on the - strand.

    :param annotation: DataFrame with Feature, transcript_id, Start, End and Strand, e.g. as
        ``read_annotation`` returns it
    :return: A copy of ``annotation`` with a RangeIndex and the computed numbers in the column
        exon_number (nullable integer). Rows that are not exon, CDS or stop_codon rows keep their
        exon_number; a CDS or stop_codon row without an overlapping exon keeps its exon_number too.
    """
    # a unique index for _assign_exon_numbers_to_cds
    annotation_df = annotation.reset_index(drop=True)

    # The computed numbers are ints, so hold the column as nullable integer. A caller may give it as
    # text, e.g. "1", with missing values on features without one.
    if "exon_number" in annotation_df.columns:
        annotation_df["exon_number"] = annotation_df["exon_number"].astype("Int64")
    else:
        annotation_df["exon_number"] = pd.Series(pd.NA, index=annotation_df.index, dtype="Int64")

    # Step 1: Compute exon numbers for exon features, by genomic order per transcript and strand
    exons = annotation_df[(annotation_df["Feature"] == "exon") & annotation_df["transcript_id"].notna()]
    sign = np.where(exons["Strand"] == "+", 1, -1)
    order = exons.assign(_key=exons["Start"].to_numpy() * sign).sort_values(["transcript_id", "_key"], kind="stable")
    annotation_df.loc[order.index, "exon_number"] = (order.groupby("transcript_id").cumcount() + 1).astype("Int64")

    # Step 2: Assign exon numbers to CDS and stop_codon features
    return _assign_exon_numbers_to_cds(annotation_df)


def _assign_exon_numbers_to_cds(df):
    """
    Gives every CDS and stop_codon row of ``df`` the exon_number of the exon of its transcript that
    it overlaps most (the first one on a tie). Rows without an overlapping exon keep their
    exon_number. The exon rows must have an exon_number already.

    :param df: DataFrame with Feature, transcript_id, Start, End and an Int64 exon_number column;
        the index must be unique.
    :return: A copy of ``df`` with the exon_number of the CDS and stop_codon rows set
    """
    df = df.copy()
    exons = df.loc[df["Feature"] == "exon", ["transcript_id", "Start", "End", "exon_number"]]
    cds = df.loc[df["Feature"].isin(["CDS", "stop_codon"]), ["transcript_id", "Start", "End"]]
    cds = cds.rename_axis("cds_index").reset_index()
    pairs = cds.merge(exons, on="transcript_id", suffixes=("", "_exon"))
    pairs = pairs[(pairs["Start_exon"] <= pairs["End"]) & (pairs["End_exon"] >= pairs["Start"])]
    overlap = np.minimum(pairs["End"], pairs["End_exon"]) - np.maximum(pairs["Start"], pairs["Start_exon"])
    best = pairs.loc[overlap.groupby(pairs["cds_index"]).idxmax()]
    df.loc[best["cds_index"].to_numpy(), "exon_number"] = best["exon_number"].to_numpy()
    return df
