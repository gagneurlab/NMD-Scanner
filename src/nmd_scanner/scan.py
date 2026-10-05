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

# A GTF is rejected by its file name
GTF_SUFFIXES = (".gtf", ".gtf.gz")
GTF_NOT_SUPPORTED = "GTF input is no longer supported. Use the GFF3 of the same GENCODE or Ensembl release."
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
    #  (especially the inclusion of the Variants into the reference CDS sequence to create the alternative CDS)

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
    except (ValueError, pl.exceptions.ComputeError) as error:
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
    :raises ValueError: for a ``.gtf`` file, since GTF input is no longer supported, and for any
        other suffix
    """
    path = os.fspath(path)
    lowered = path.lower()
    if lowered.endswith(GTF_SUFFIXES):
        raise ValueError(f"Cannot read {path!r}: {GTF_NOT_SUPPORTED}")
    if lowered.endswith(GFF3_SUFFIXES):
        return "gff3"
    raise ValueError(f"Cannot detect annotation format from filename: {path!r}. Expected one of {GFF3_SUFFIXES}.")


def read_annotation(path, fasta=None, reassign_exons=False):
    """
    Reads a GFF3 gene annotation file into a pandas DataFrame with the exon rows and the coding
    regions that ``rules.extract_ptc`` takes (see ``read_gff3``).

    The coding regions are the CDS rows: they include the stop codon, one row per transcript and
    exon. On them, the column has_stop_codon says whether the coding region of the transcript ends
    in an annotated stop codon. On the exon rows, has_stop_codon is NA.

    The format is checked from the filename suffix (``.gff3`` or ``.gff``, gzip-compressed or not).
    GTF input is no longer supported.

    :param path: Path to the GFF3 file
    :param fasta: Reference genome (pyfaidx.Fasta object), which shows whether a CDS ends in a stop
        codon. Required: without it, read_annotation raises a ValueError.
    :param reassign_exons: Recompute the exon numbers with ``compute_exon_numbers``
    :return: DataFrame, as ``read_gff3`` returns it
    :raises ValueError: for a GTF file name, for an unknown suffix, or without ``fasta``
    """
    detect_annotation_format(path)
    if fasta is None:
        raise ValueError(
            "Reading a GFF3 needs the reference genome FASTA, which shows whether a CDS ends in a stop codon."
        )

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
      (as in the GTF; e.g. hg19 lift37 has ``ENST00000477874.1_2`` there but ``ENST00000477874.1``
      in ``ID``). Only chrY PAR transcripts reuse the chrX value on that attribute (the "_PAR_Y"
      suffix only shows up in ``ID``/``Parent``), so for those the ``ID``/``Parent`` hierarchy gives
      ``gene_id``/``transcript_id``, which is what the GTF has too.
    - Ensembl: ``gene_id``/``transcript_id`` are never given directly on exon/CDS rows and are
      resolved from the ``ID``/``Parent`` hierarchy (``gene:``/``transcript:`` prefixed). Exon rows
      carry their number as ``rank``; a CDS row takes the number of the exon it lies in. Without
      ``rank``, ``compute_exon_numbers`` computes the numbers.

    The coding regions are the CDS rows. A GFF3 CDS includes the stop codon, so the result has no
    stop_codon rows. On the CDS rows, the column has_stop_codon says whether the coding region of
    the transcript ends in an annotated stop codon. On the exon rows, it is NA. Both flavors get the
    coding regions and has_stop_codon of the GTF of the same release:

    - GENCODE: has_stop_codon is True if the transcript has stop_codon rows. A few cds_end_NF
      transcripts have 3 more CDS bases than in the GTF; ``_trim_cds_end_nf_stop_codons`` removes
      them.
    - Ensembl: the GFF3 has no stop_codon rows, so has_stop_codon comes from the last 3 CDS bases
      in ``fasta`` (``_has_stop_codon_from_sequence``).

    The file is read with polars-bio (see ``_read_gff3_rows``), plain or gzip-compressed. Only the
    attributes in ``GFF3_ATTRIBUTES`` are read.

    :param gff3_path: Path to the GFF3 file
    :param fasta: Reference genome (pyfaidx.Fasta object)
    :return: DataFrame with the exon and CDS rows in file order and a RangeIndex, the columns
        Chromosome, Source, Feature, Start, End, Score, Strand, Frame, has_stop_codon, gene_id and
        transcript_id, and the other attributes of ``GFF3_ATTRIBUTES`` that the file has, except ID
        and Parent. exon_number (nullable integer) is there for an Ensembl GFF3, and for a GENCODE
        GFF3 that has the attribute.
    :raises ValueError: if the file has no ID or Parent attribute, e.g. because it is a GTF, or if its
        flavor is neither GENCODE nor Ensembl
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
            f"Cannot read {os.fspath(gff3_path)!r} as GFF3: no row has an ID or Parent attribute. "
            f"If it is a GTF: {GTF_NOT_SUPPORTED}"
        )

    if {"gene_type", "transcript_type"} <= columns:
        df = _trim_cds_end_nf_stop_codons(_normalize_gencode_gff3(df), fasta)
        # the CDS includes the stop codon; the stop_codon rows only say whether there is one
        is_stop = df["Feature"] == "stop_codon"
        if not is_stop.any() and (df["Feature"] == "CDS").any():
            logger.warning(
                "No stop_codon rows found next to the CDS rows: every transcript is treated as having no annotated "
                "stop codon (no 3'UTR length, no stop codon distance, every in-frame stop is premature). "
                "A GENCODE GFF3 marks the stop codons with stop_codon rows."
            )
        df = _set_has_stop_codon(df[~is_stop], df.loc[is_stop, "transcript_id"])
    elif "biotype" in columns:
        df = _has_stop_codon_from_sequence(_normalize_ensembl_gff3(df), fasta)
    else:
        raise ValueError(
            "Unrecognized GFF3 flavor: expected GENCODE-style attributes "
            "(gene_type/transcript_type) or Ensembl-style attributes (biotype)."
        )

    # ID/Parent are GFF3 structural columns with no GTF equivalent; they served their purpose
    # in _resolve_ids_via_hierarchy. Drop them so a stray "ID" column doesn't collide with the
    # unrelated "ID" column read_vcf uses for the VCF record ID once CDS rows are joined to it.
    df = df.drop(columns=[c for c in ("ID", "Parent") if c in df.columns])

    return df.reset_index(drop=True)


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
GFF3_FEATURES = ("transcript", "exon", "CDS", "stop_codon")


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
    except (ValueError, pl.exceptions.ComputeError) as error:
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
    Sets has_stop_codon on the CDS rows of an Ensembl GFF3 as the Ensembl GTF has it.

    The Ensembl GTF has stop_codon rows for a transcript if the last 3 bases of its CDS are a stop
    codon (see ``_last_codons``), unless its last coding exon ends mid-codon or it is tagged
    cds_end_NF. The GFF3 gives the first condition by the exon attribute ``ensembl_end_phase`` (1 or
    2), but has no cds_end_NF tag. So a cds_end_NF transcript whose CDS ends in stop codon bases gets
    has_stop_codon True here, but has no stop codon in the GTF: 13 transcripts in Ensembl 108, none
    on chr22.

    :param df: Exon and CDS rows of an Ensembl GFF3 (DataFrame) with transcript_id and exon_number
    :param fasta: Reference genome (pyfaidx.Fasta object)
    :return: DataFrame with the column has_stop_codon (see ``_set_has_stop_codon``)
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
        "gets a stop codon, unlike in the Ensembl GTF (13 transcripts in Ensembl 108)."
    )
    return _set_has_stop_codon(df, stops["transcript_id"])


def _set_has_stop_codon(df, stop_transcripts):
    """
    Returns a copy of ``df`` with the column has_stop_codon (nullable boolean): on the CDS rows,
    whether their transcript is in ``stop_transcripts``; NA on the other rows.
    """
    has_stop = df["transcript_id"].isin(stop_transcripts).astype("boolean")
    return df.assign(has_stop_codon=has_stop.where(df["Feature"] == "CDS"))


def _trim_cds_end_nf_stop_codons(df, fasta):
    """
    Removes the last codon from the CDS of a GENCODE GFF3 where the GENCODE GTF has it as UTR.

    This concerns transcripts tagged cds_end_NF, without stop_codon rows, whose CDS ends in a
    complete stop codon (see ``_last_codons``): 10 transcripts in GENCODE 42, none on chr22. All
    other CDS and stop_codon rows are the same as in the GTF. GENCODE GFF3 repeats the transcript
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
        "it is removed from the CDS, as in the GENCODE GTF.",
        len(stops),
    )
    # a CDS row that held only stop codon bases is gone
    return df[df["Start"] < df["End"]]


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
        parts, missing = [], 3
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
    the hierarchy only fills in the rows without one and the "_PAR_Y" rows. The attribute is the
    GTF's value; the ID can differ from it (GENCODE lift37: ``ENST00000477874.1_2`` vs
    ``ENST00000477874.1``).

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
    Returns the exon, CDS and stop_codon rows of a GENCODE-flavored GFF3 dataframe (see
    ``read_gff3``), with gene_id/transcript_id from the attributes
    (from ID/Parent where missing and for chrY PAR copies).
    """
    transcript_mask = df["Feature"] == "transcript"
    df = _resolve_ids_via_hierarchy(df, transcript_mask, keep_attributes=True)

    child_mask = df["Feature"].isin(["exon", "CDS", "stop_codon"])
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

    On + Strand: Smallest exon number is the Start, Largest exon number is the End.
    On - Strand: Smallest exon number is the Start, Largest exon number is the End.
    (was different for hg19: the smallest exon number was the end, that is why we need to adjust it here.)

    :param annotation: DataFrame with Feature, transcript_id, Start, End and Strand, e.g. as
        ``read_annotation`` returns it
    :return: A copy of ``annotation`` with a RangeIndex and the computed numbers in the column
        exon_number (nullable integer). Rows that are not exon, CDS or stop_codon rows keep their
        exon_number; a CDS or stop_codon row without an overlapping exon keeps its exon_number too.
    """
    # a unique index for _assign_exon_numbers_to_cds
    gtf_df = annotation.reset_index(drop=True)

    # The computed numbers are ints, so hold the column as nullable integer. A caller may give it as
    # text, e.g. "1", with missing values on features without one.
    if "exon_number" in gtf_df.columns:
        gtf_df["exon_number"] = gtf_df["exon_number"].astype("Int64")
    else:
        gtf_df["exon_number"] = pd.Series(pd.NA, index=gtf_df.index, dtype="Int64")

    # Step 1: Compute exon numbers for exon features, by genomic order per transcript and strand
    exons = gtf_df[(gtf_df["Feature"] == "exon") & gtf_df["transcript_id"].notna()]
    sign = np.where(exons["Strand"] == "+", 1, -1)
    order = exons.assign(_key=exons["Start"].to_numpy() * sign).sort_values(["transcript_id", "_key"], kind="stable")
    gtf_df.loc[order.index, "exon_number"] = (order.groupby("transcript_id").cumcount() + 1).astype("Int64")

    # Step 2: Assign exon numbers to CDS and stop_codon features
    return _assign_exon_numbers_to_cds(gtf_df)


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
