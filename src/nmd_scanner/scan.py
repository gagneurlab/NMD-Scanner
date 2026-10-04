# Import dependencies
import itertools
import logging
import os

import numpy as np
import pandas as pd
import pyranges as pr
from Bio.Seq import Seq
from pyfaidx import Fasta

# Create the functions used for reading in the files (VCF, GTF, FASTA)

GTF_SUFFIXES = (".gtf", ".gtf.gz")
GFF3_SUFFIXES = (".gff3", ".gff3.gz", ".gff", ".gff.gz")

STOP_CODONS = {"TAA", "TAG", "TGA"}
MITOCHONDRIAL_CHROMOSOMES = {"MT", "M", "chrM", "chrMT"}
# vertebrate mitochondrial code (NCBI translation table 2): AGA and AGG are stop codons, TGA codes for Trp
MITOCHONDRIAL_STOP_CODONS = {"TAA", "TAG", "AGA", "AGG"}

logger = logging.getLogger(__name__)


def read_vcf(vcf_path):
    # TODO: adjust this function to also get structural variants and then also adjust downstream analysis
    #  (especially the inclusion of the Variants into the reference CDS sequence to create the alternative CDS)

    """
    Read a single VCF file into a PyRanges object with adjusted coordinates.

    The VCF must be left-normalized and single-allelic (one ALT allele per
    record), e.g. produced by ``bcftools norm -m- -f reference.fa``.
    Multi-allelic records (comma-separated ALT) are rejected because the
    downstream variant application assumes exactly one ALT allele per row.
    """
    df = pd.read_csv(
        vcf_path,
        comment="#",
        sep="\t",
        header=None,
        names=["Chromosome", "Start", "ID", "Ref", "Alt", "Qual", "Filter", "Info"],
    )

    # Reject multi-allelic records: the VCF spec allows comma-separated ALT,
    # but the rest of the pipeline assumes one ALT allele per row.
    multiallelic = df["Alt"].astype(str).str.contains(",")
    if multiallelic.any():
        n_multiallelic = int(multiallelic.sum())
        raise ValueError(
            f"{n_multiallelic} multi-allelic record(s) found in {vcf_path}. "
            "NMD-Scanner requires a left-normalized, single-allelic VCF. "
            "Split and normalize first, e.g. `bcftools norm -m- -f reference.fa`."
        )

    # Adjust coordinates to 0-based
    df["Start"] = df["Start"] - 1
    df["End"] = df["Start"] + df["Ref"].str.len()

    # Keep only relevant columns
    gr = pr.PyRanges(df[["Chromosome", "Start", "End", "ID", "Ref", "Alt", "Qual", "Filter", "Info"]])
    return gr


def read_gtf(gtf_path):
    """
    Reads a GTF file into a PyRanges object, with its rows as they are in the file. A GTF CDS excludes
    the stop codon, which has its own stop_codon rows. ``merge_stop_codons_into_cds`` builds the coding
    regions that ``rules.extract_ptc`` takes. ``read_annotation`` reads a GTF straight into them.
    """
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")
    try:
        return pr.read_gtf(gtf_path)
    except ValueError as error:
        raise ValueError(
            f"Cannot read {os.fspath(gtf_path)!r} as GTF ({error}). "
            "If it is a GFF3, pass it with --annotation (or read_annotation), under a .gff3 file name."
        ) from error


def merge_stop_codons_into_cds(df, transcript_col="transcript_id"):
    """
    Returns the coding regions of a GTF, in the form ``rules.extract_ptc`` takes: one CDS row per
    transcript and exon, the union of its CDS and stop_codon rows.

    A GTF CDS excludes the stop codon, which has its own stop_codon rows. A stop codon split across an
    intron has two of them, and the second one can lie in an exon without CDS. A transcript without
    stop_codon rows, e.g. one tagged cds_end_NF, keeps its CDS as it is: its coding region ends without
    a stop codon.

    The rows are merged per transcript and exon_number. To recompute the exon numbers with
    ``compute_exon_numbers``, do so before the merge.

    :param df: GTF rows (DataFrame) with Feature, Start, End, exon_number and ``transcript_col``. Only
        the CDS and stop_codon rows are used.
    :param transcript_col: The name of the column that indicates the transcript ID
    :return: DataFrame with the CDS rows, each extended by the stop codon bases of its exon, plus one
        CDS row for each exon that holds only stop codon bases. exon_number is int. The column
        has_stop_codon says whether the coding region of the transcript ends in an annotated stop codon,
        i.e. whether the transcript has stop_codon rows.
    :raises ValueError: if stop codon bases do not touch or overlap the CDS of their exon.
    """
    df = df[df["Feature"].isin(["CDS", "stop_codon"])].copy()
    df["exon_number"] = df["exon_number"].astype(int)
    keys = [transcript_col, "exon_number"]

    is_stop = df["Feature"] == "stop_codon"
    stop_rows = df[is_stop]
    if stop_rows.empty and not df.empty:
        logger.warning(
            "No stop_codon rows found next to the CDS rows: every transcript is treated as having no annotated "
            "stop codon (no 3'UTR length, no stop codon distance, every in-frame stop is premature). "
            "A GTF CDS excludes the stop codon, so the merge needs the stop_codon rows of the GTF."
        )
    stops = stop_rows.groupby(keys, observed=True).agg(stop_start=("Start", "min"), stop_end=("End", "max"))

    cds = df[~is_stop].merge(stops, left_on=keys, right_index=True, how="left")
    has_stop = cds["stop_start"].notna()
    # coordinates are half-open, so touching intervals share one coordinate
    gap = has_stop & ((cds["stop_start"] > cds["End"]) | (cds["stop_end"] < cds["Start"]))
    if gap.any():
        raise ValueError(
            "stop_codon rows do not touch the CDS of their exon in transcripts: "
            + ", ".join(sorted(cds.loc[gap, transcript_col].astype(str).unique()))
        )
    start_dtype, end_dtype = cds["Start"].dtype, cds["End"].dtype
    cds["Start"] = cds[["Start", "stop_start"]].min(axis=1).astype(start_dtype)
    cds["End"] = cds[["End", "stop_end"]].max(axis=1).astype(end_dtype)
    cds = cds.drop(columns=["stop_start", "stop_end"])

    # exons with stop codon bases but no CDS row, e.g. the second part of a split stop codon
    cds_keys = pd.MultiIndex.from_frame(cds[keys])
    stop_only = stops[~stops.index.isin(cds_keys)]
    extra = stop_rows.drop_duplicates(keys).set_index(keys).loc[stop_only.index].reset_index()
    extra["Start"] = stop_only["stop_start"].to_numpy().astype(start_dtype)
    extra["End"] = stop_only["stop_end"].to_numpy().astype(end_dtype)
    extra["Feature"] = "CDS"

    logger.info(
        "Stop codons from stop_codon rows: %d transcripts, %d of them with stop codon bases in an exon without CDS.",
        stops.index.get_level_values(transcript_col).nunique(),
        extra[transcript_col].nunique(),
    )
    coding = pd.concat([cds, extra[cds.columns]], ignore_index=True)
    coding["has_stop_codon"] = coding[transcript_col].isin(stop_rows[transcript_col])
    return coding


def detect_annotation_format(path):
    """
    Detects whether ``path`` (str or path-like) is a GTF or a GFF3 file from its filename suffix
    (gzip-compressed or not).
    """
    path = os.fspath(path)
    lowered = path.lower()
    if lowered.endswith(GTF_SUFFIXES):
        return "gtf"
    if lowered.endswith(GFF3_SUFFIXES):
        return "gff3"
    raise ValueError(
        f"Cannot detect annotation format from filename: {path!r}. Expected one of {GTF_SUFFIXES + GFF3_SUFFIXES}."
    )


def read_annotation(path, fasta=None, fmt=None, reassign_exons=False):
    """
    Reads a gene annotation file, GTF or GFF3, into a PyRanges object with the coding regions that
    ``rules.extract_ptc`` takes.

    The coding regions are the CDS rows: they include the stop codon, one row per transcript and
    exon. On them, the column has_stop_codon says whether the coding region of the transcript ends
    in an annotated stop codon. On the other rows, has_stop_codon is NA. A GTF CDS excludes the stop
    codon, so ``merge_stop_codons_into_cds`` merges the stop_codon rows of a GTF into its CDS rows.
    The other rows of a GTF stay as they are, and its exon_number becomes Int64. A GFF3 CDS
    includes the stop codon; see ``read_gff3``.

    The format is auto-detected from the filename suffix (``.gtf``/``.gff3``/``.gff``,
    gzip-compressed or not), unless ``fmt`` is given.

    :param path: Path to the GTF or GFF3 file
    :param fasta: Reference genome (pyfaidx.Fasta object). Required for GFF3, which takes from it
        whether a CDS ends in a stop codon; unused for GTF.
    :param fmt: "gtf" or "gff3" to skip the detection from the filename suffix
    :param reassign_exons: Recompute the exon numbers with ``compute_exon_numbers``. For a GTF, this
        runs before the merge, which keys on exon_number.
    """
    if fmt is None:
        fmt = detect_annotation_format(path)
    elif fmt not in ("gtf", "gff3"):
        raise ValueError(f"Unknown annotation format {fmt!r}, expected 'gtf' or 'gff3'.")
    if fmt == "gff3" and fasta is None:
        raise ValueError(
            "Reading a GFF3 needs the reference genome FASTA, which shows whether a CDS ends in a stop codon."
        )

    annotation = read_gtf(path) if fmt == "gtf" else read_gff3(path, fasta)
    if reassign_exons:
        logger.info("Recomputing exon numbers.")
        annotation = compute_exon_numbers(annotation)
    if fmt == "gff3":
        return annotation

    df = annotation.df
    # GTF rows hold exon_number as str, the coding regions as int: cast both before the concat
    other = df[~df["Feature"].isin(["CDS", "stop_codon"])].astype({"exon_number": "Int64"})
    other["has_stop_codon"] = pd.Series(pd.NA, index=other.index, dtype="boolean")
    coding = merge_stop_codons_into_cds(df).astype({"exon_number": "Int64", "has_stop_codon": "boolean"})
    return pr.PyRanges(pd.concat([other, coding], ignore_index=True))


def read_gff3(gff3_path, fasta):
    """
    Reads a GFF3 file into a PyRanges object with the exon rows and the coding regions, in the
    column layout of ``read_gtf``.

    Two GFF3 flavors are supported, auto-detected from the attributes present:

    - GENCODE: every feature row carries ``gene_id``/``transcript_id`` directly, and they are kept
      (as in the GTF; e.g. hg19 lift37 has ``ENST00000477874.1_2`` there but ``ENST00000477874.1``
      in ``ID``). Only chrY PAR transcripts reuse the chrX value on that attribute (the "_PAR_Y"
      suffix only shows up in ``ID``/``Parent``), so for those the ``ID``/``Parent`` hierarchy gives
      ``gene_id``/``transcript_id``, which is what the GTF has too.
    - Ensembl: ``gene_id``/``transcript_id`` are never given directly on exon/CDS rows and are
      resolved from the ``ID``/``Parent`` hierarchy (``gene:``/``transcript:`` prefixed). Exon rows
      carry their number as ``rank``; a CDS row takes the number of the exon it lies in. Without
      ``rank``, ``compute_exon_numbers`` computes the numbers.

    The coding regions are the CDS rows, in the form ``read_annotation`` returns for a GTF too. A
    GFF3 CDS includes the stop codon, so it needs no merge with stop_codon rows, and the result has
    no stop_codon rows. On the CDS rows, the column has_stop_codon says whether the coding region of
    the transcript ends in an annotated stop codon. On the exon rows, it is NA. Both flavors get the
    coding regions and has_stop_codon of the GTF of the same release:

    - GENCODE: has_stop_codon is True if the transcript has stop_codon rows. A few cds_end_NF
      transcripts have 3 more CDS bases than in the GTF; ``_trim_cds_end_nf_stop_codons`` removes
      them.
    - Ensembl: the GFF3 has no stop_codon rows, so has_stop_codon comes from the last 3 CDS bases
      in ``fasta`` (``_has_stop_codon_from_sequence``).

    :param gff3_path: Path to the GFF3 file
    :param fasta: Reference genome (pyfaidx.Fasta object)
    """
    if not os.path.exists(gff3_path):
        raise FileNotFoundError(f"GFF3 file not found: {gff3_path}")

    try:
        df = pr.read_gff3(gff3_path).df
    except ValueError as error:
        raise ValueError(
            f"Cannot read {os.fspath(gff3_path)!r} as GFF3 ({error}). If it is a GTF, give it a .gtf file name."
        ) from error
    columns = set(df.columns)

    if {"gene_type", "transcript_type"} <= columns:
        df = _trim_cds_end_nf_stop_codons(_normalize_gencode_gff3(df), fasta)
        # the CDS includes the stop codon; the stop_codon rows only say whether there is one
        is_stop = df["Feature"] == "stop_codon"
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

    return pr.PyRanges(df)


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
        "gets a stop codon, unlike in the Ensembl GTF (13 transcripts in Ensembl 108). Use the GTF to avoid this."
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

    return df.loc[child_mask].reset_index(drop=True)


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
    return compute_exon_numbers(pr.PyRanges(child_df)).df


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


def compute_exon_numbers(gtf):
    """
    Compute exon numbers for the Features exon, CDS and stop_codon in a GTF PyRanges object.
    Exon numbers are assigned based on genomic order per transcript and strand.
    CDS and stop_codon features inherit the exon number of the exon they overlap.

    On + Strand: Smallest exon number is the Start, Largest exon number is the End.
    On - Strand: Smallest exon number is the Start, Largest exon number is the End.
    (was different for hg19: the smallest exon number was the end, that is why we need to adjust it here.)

    :param gtf: PyRanges object of the GTF
    :return: PyRanges object with new column 'exon_number_computed'
    """
    gtf_df = gtf.df.copy()

    # A GTF read from file has exon_number as str (pandas 3) with missing values on features
    # without one. The computed numbers are ints, so hold the column as nullable integer.
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
    return pr.PyRanges(_assign_exon_numbers_to_cds(gtf_df))


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
