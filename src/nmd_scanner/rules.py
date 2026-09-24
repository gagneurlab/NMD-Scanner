# Import dependencies

import logging

import numpy as np
import pandas as pd
from Bio.Seq import Seq

from nmd_scanner import catch_sequence
from nmd_scanner._polars_bio import pb
from nmd_scanner.schema import PTC_COLUMN_KINDS, apply_schema, empty_table
from nmd_scanner.variant_placement import ReferenceSequence, place_in_transcript, variant_placements

logger = logging.getLogger(__name__)

# An ALT allele that names a structural variant instead of giving its sequence: a symbolic allele such as <DEL>,
# <DUP:TANDEM> or <*>, a breakend such as G]chr2:100] or [chr2:100[G, or a single breakend such as G. or .G
SYMBOLIC_ALT_PATTERN = r"^<[^>]*>$|[\[\]]|^\.[A-Za-z]+$|^[A-Za-z]+\.$"


# Main extract PTC script:
def extract_ptc(cds_df, vcf, fasta, exons_df):
    """
    Main function for extracting reference coding sequence, alternative coding sequence by incorporating the variant, analyzing for premature termination codons (PTCs),
    start and stop loss, and getting the transcript information.

    :param cds_df: Coding regions of the annotation (DataFrame): CDS rows that include the stop codon, one row per
                   transcript and exon, with exon_number and the column has_stop_codon. has_stop_codon says whether
                   the coding region of the transcript ends in an annotated stop codon. ``scan.read_annotation``
                   returns the coding regions as its CDS rows.
    :param vcf: Variants (DataFrame) with Chromosome, Start, End, ID, Ref and Alt, as ``scan.read_vcf`` returns them.
                A record with a symbolic ALT allele or a breakend is skipped (see ``drop_symbolic_alleles``).
    :param fasta: Reference genome sequence (pyfaidx.Fasta object)
    :param exons_df: Exon rows of the annotation (DataFrame)
    :return: analyze_transcript_df: Annotated dataframe with ref and alt CDS information, PTC analysis, start & stop loss analysis and transcript information.
             It has the columns and dtypes of PTC_COLUMN_KINDS (see nmd_scanner.schema). It has one row per
             variant and transcript where the variant touches the coding region or the splice dinucleotide at one
             of its exon edges, also through an equivalent placement of an indel (see
             ``variant_placement.place_in_transcript``). If the alt transcript is unknown, unknown_reason names why,
             and the alt columns are null. It has zero rows if no variant touches a coding region, or every variant
             is skipped or has a reference mismatch.
    :raises ValueError: if cds_df has no has_stop_codon column, i.e. it does not hold the coding regions.
    """

    if "has_stop_codon" not in cds_df.columns:
        raise ValueError(
            "cds_df has no has_stop_codon column. extract_ptc takes the coding regions: CDS rows that include the "
            "stop codon, with has_stop_codon, as scan.read_annotation returns them."
        )

    cds_df_adj = cds_df.copy()
    # annotation attributes are read as text; exon numbers are int in the output tuples (e.g. ref_cds_info)
    cds_df_adj["exon_number"] = cds_df_adj["exon_number"].astype(int)

    # The variant application below cannot apply a structural variant
    vcf = drop_symbolic_alleles(vcf)

    # Join the variants with the coding regions and the splice dinucleotides at their exon edges
    references = {}

    def reference(chromosome):
        if chromosome not in references:
            references[chromosome] = ReferenceSequence(fasta, chromosome)
        return references[chromosome]

    variants, placements = place_variants(vcf, set(cds_df_adj["Chromosome"]), reference)
    intersection_cds_vcf = join_variant_windows(cds_df_adj, variants)
    logger.info("Joining variants with cds entries: done.")

    # Nothing to analyze: the steps below need at least one row
    if intersection_cds_vcf.empty:
        logger.info("No variant overlapped a CDS; there are no results to compute.")
        return empty_table(PTC_COLUMN_KINDS)

    ##########################################################################################
    # TODO: fix minus strand variants (only for TCGA and MMRF VCF!)
    # Fix REF and ALT for minus-strand CDSs
    # mask_minus_strand = intersection_cds_vcf["Strand"] == "-"
    # intersection_cds_vcf.loc[mask_minus_strand, "Ref"] = intersection_cds_vcf.loc[mask_minus_strand, "Ref"].apply(
    #   lambda seq: str(Seq(seq).reverse_complement()))
    # intersection_cds_vcf.loc[mask_minus_strand, "Alt"] = intersection_cds_vcf.loc[mask_minus_strand, "Alt"].apply(
    #   lambda seq: str(Seq(seq).reverse_complement()))
    ##########################################################################################

    # Filter out Variants with a reference mismatch
    mismatched_rows = intersection_cds_vcf[~intersection_cds_vcf["Ref_matches"]].drop_duplicates(
        ["transcript_id", "variant_row"]
    )
    if not mismatched_rows.empty:
        logger.warning("Skipping %d variant-transcript pairs due to reference mismatches.", len(mismatched_rows))
        logger.warning(
            "Reference-mismatched variants:\n%s",
            mismatched_rows[["transcript_id", "Chromosome", "Start_variant", "End_variant", "Ref", "Alt"]].to_string(
                index=False
            ),
        )
    intersection_cds_vcf = intersection_cds_vcf[intersection_cds_vcf["Ref_matches"]].copy()

    if intersection_cds_vcf.empty:
        logger.info("No variant left after the reference check; there are no results to compute.")
        return empty_table(PTC_COLUMN_KINDS)

    # Apply each variant to the coding rows of each transcript, or find why its alt transcript is unknown
    intersection_cds_vcf = apply_variants(intersection_cds_vcf, placements, cds_df_adj, exons_df, reference)
    logger.info("Creating alt CDS sequence: done.")

    if intersection_cds_vcf.empty:
        logger.info("No variant changed a CDS or its splice sites; there are no results to compute.")
        return empty_table(PTC_COLUMN_KINDS)

    # Limit to relevant transcript (to save time)
    relevant_transcripts = intersection_cds_vcf["transcript_id"].unique()
    cds_df_adj = cds_df_adj[cds_df_adj["transcript_id"].isin(relevant_transcripts)].copy()

    # Fetch reference sequence for all CDS entries per (relevant) transcripts
    cds_df_adj = catch_sequence.add_exon_cds_sequence(cds_df_adj, fasta)  # for faster access

    # get full reference CDS per transcript by stiching exon CDS regions, plus alternative CDS with CDS exon information for both ref and alt
    results_df = create_reference_cds(intersection_cds_vcf, cds_df_adj)
    logger.info("Create reference CDS: done.")

    # Get transcript sequence for relevant transcripts (speed up process) + length and transcript exon information (Tuple: exon number & exon length)
    exons_df = exons_df[exons_df["transcript_id"].isin(relevant_transcripts)].copy()
    # annotation attributes are read as text; align with cds_df_adj so exon numbers are int
    # everywhere they end up together in a tuple (e.g. transcript_exon_info, *_stop_codon_exons).
    exons_df["exon_number"] = exons_df["exon_number"].astype(int)
    exon_seqs = get_transcript_sequence(exons_df, fasta)
    logger.info("Get transcript sequence: done.")

    # Locate the coding region (CDS plus stop codon) in each transcript sequence
    exons_by_transcript = dict(list(exons_df.groupby("transcript_id")))
    cds_ranges = {
        transcript_id: cds_range_in_transcript(exons_by_transcript[transcript_id], cds_group)
        for transcript_id, cds_group in cds_df_adj.groupby("transcript_id")
        if transcript_id in exons_by_transcript
    }

    # Validate that the CDS is present inside the transcript sequence, to make sure the transcript sequence was computed correctly
    exon_seqs_indexed = exon_seqs.set_index("transcript_id")

    def check_cds_in_transcript(row):
        transcript_id = row["transcript_id"]

        # Skip if transcript_id not found
        if transcript_id not in exon_seqs_indexed.index:
            return False

        transcript_seq = exon_seqs_indexed.loc[transcript_id, "transcript_sequence"]
        ref_cds_seq = row["ref_cds_seq"]

        # Check if CDS is a substring of the transcript
        return ref_cds_seq in transcript_seq

    results_df["cds_in_transcript"] = results_df.apply(check_cds_in_transcript, axis=1)

    # TODO: Analyze reference and alternative CDS for start / stop codons
    analysis_df = analyze_sequence(results_df)
    loss_df = start_stop_loss(analysis_df)
    logger.info("Analyzing sequence: done.")

    # Annotate transcript information (transcript start, end, sequence, length, exon info) in case of start or stop loss
    # transcript sequences are in: exon_seqs_subset
    logger.info("Annotating transcript information in case of start/stop loss.")
    transcript_starts = exon_seqs.set_index("transcript_id")["start"].to_dict()
    loss_df["transcript_start"] = loss_df["transcript_id"].map(transcript_starts)
    transcript_ends = exon_seqs.set_index("transcript_id")["end"].to_dict()
    loss_df["transcript_end"] = loss_df["transcript_id"].map(transcript_ends)
    transcript_sequences = exon_seqs.set_index("transcript_id")[
        "transcript_sequence"
    ].to_dict()  # create map of transcript-id to transcript sequence
    loss_df["transcript_seq"] = loss_df["transcript_id"].map(transcript_sequences)
    transcript_lengths = exon_seqs.set_index("transcript_id")["transcript_length"].to_dict()
    loss_df["transcript_length"] = loss_df["transcript_id"].map(transcript_lengths)
    # Object dtype keeps the positions as int next to None.
    loss_df[["cds_start_in_transcript", "cds_end_in_transcript"]] = pd.DataFrame(
        [cds_ranges.get(transcript_id) or (None, None) for transcript_id in loss_df["transcript_id"]],
        index=loss_df.index,
        columns=["cds_start_in_transcript", "cds_end_in_transcript"],
        dtype=object,
    )

    # Add exon information to dataframe
    transcript_exon_info = exon_seqs.set_index("transcript_id")["transcript_exon_info"].to_dict()
    loss_df["transcript_exon_info"] = loss_df["transcript_id"].map(transcript_exon_info)

    # Splice alternative CDS into reference transcript sequence to create alternative transcript sequence and measure new length
    loss_df["alt_transcript_seq"] = loss_df.apply(
        lambda row: (
            splice_alt_cds_into_transcript(row, row["transcript_seq"])
            if pd.notnull(row["transcript_seq"]) and pd.notnull(row["alt_cds_seq"])
            else None
        ),
        axis=1,
    )
    loss_df["alt_transcript_length"] = pd.Series(
        [len(seq) if isinstance(seq, str) else None for seq in loss_df["alt_transcript_seq"]],
        index=loss_df.index,
        dtype=object,
    )
    # The alt CDS starts where the ref CDS starts, shifted by the length change of the 5'UTR
    loss_df["alt_cds_start_in_transcript"] = pd.Series(
        [
            cds_start - len(utr5_ref) + len(utr5_alt) if isinstance(seq, str) else None
            for seq, cds_start, (utr5_ref, utr5_alt) in zip(
                loss_df["alt_transcript_seq"], loss_df["cds_start_in_transcript"], loss_df["utr5_change"]
            )
        ],
        index=loss_df.index,
        dtype=object,
    )
    loss_df = loss_df.drop(columns=["utr5_change", "utr3_change"])

    # Classify the first in-frame stop codon of the alternative transcript, and analyze the transcript sequence
    # (e.g., frame, length, stop codon position, etc.) in case of start or stop loss
    analyze_transcript_df = analyze_transcript(loss_df)

    return apply_schema(analyze_transcript_df, PTC_COLUMN_KINDS)


# Functions used for extracting PTC:


# polars-bio joins the intervals as 32-bit integers with sign
MAX_JOIN_COORDINATE = 2**31 - 1


def drop_symbolic_alleles(vcf):
    """
    Returns the variants without the records whose ALT is a symbolic allele or a breakend (see SYMBOLIC_ALT_PATTERN),
    and logs a warning with their count. Such an ALT names a structural variant instead of giving its sequence. The
    variant application would insert its text into the sequence, e.g. "<DEL>" into the alt CDS.

    :param vcf: Variants (DataFrame) with the column Alt
    :return: The rows of ``vcf`` whose ALT is a sequence, in their order and with their index
    """

    symbolic = vcf["Alt"].astype(str).str.contains(SYMBOLIC_ALT_PATTERN, regex=True)
    if symbolic.any():
        logger.warning(
            "Skipping %d variant(s) with a symbolic ALT allele or a breakend, e.g. <DEL> or G]chr2:100]. "
            "NMD-Scanner cannot apply structural variants yet.",
            int(symbolic.sum()),
        )
    return vcf[~symbolic]


def join_variants_to_cds(cds_df, vcf):
    """
    Joins every CDS row to the variants that overlap it.

    A CDS row and a variant overlap if they are on the same Chromosome and their 0-based half-open
    intervals share at least one base: a variant that ends at the Start of a CDS row, or starts at
    its End, does not overlap it. Strand is ignored. polars-bio computes the overlaps.

    :param cds_df: CDS rows (DataFrame) with Chromosome, Start and End
    :param vcf: Variants (DataFrame) with Chromosome, Start and End
    :return: DataFrame with one row per overlapping CDS row and variant, and a RangeIndex. It has the
        columns of cds_df, then the columns of vcf except Chromosome, also if no variant overlaps a
        CDS row. A column of vcf that cds_df has too gets the suffix "_variant", e.g. Start_variant
        and End_variant. The columns keep their dtypes, except a text column of cds_df with fewer
        distinct values than half its rows, which becomes category to save memory. The rows come in
        the order of cds_df; the variants of one CDS row by Start, then by End descending, then in
        the order of vcf.
    :raises ValueError: if an End of cds_df or vcf is above 2**31 - 1, the largest coordinate that
        polars-bio joins
    """

    for name, df in (("cds_df", cds_df), ("vcf", vcf)):
        if (df["End"] > MAX_JOIN_COORDINATE).any():
            raise ValueError(
                f"Cannot join {name}: it has an End above {MAX_JOIN_COORDINATE}, the largest coordinate that polars-bio joins."
            )

    # a category column holds each distinct text once, which makes the copies of the rows below smaller
    repeated_text = [
        name
        for name, column in cds_df.items()
        if pd.api.types.is_string_dtype(column) and column.nunique() < len(cds_df) / 2
    ]
    cds_df = cds_df.astype(dict.fromkeys(repeated_text, "category"))

    def intervals(df):
        frame = pd.DataFrame(
            {
                "chrom": df["Chromosome"].astype(str).to_numpy(),
                "start": df["Start"].to_numpy(dtype="int64"),
                "end": df["End"].to_numpy(dtype="int64"),
                "row": np.arange(len(df)),
            }
        )
        # polars-bio reads this to treat the intervals as 0-based half-open
        frame.attrs["coordinate_system_zero_based"] = True
        return frame

    order = ["row_cds", "start_variant", "end_variant", "row_variant"]
    # A DataFrame, not a LazyFrame: the collect of a LazyFrame calls logging.info(), which configures
    # the root logger, and shows a tqdm bar. The sort copies only the 4 columns selected before it.
    pairs = (
        pb.overlap(intervals(cds_df), intervals(vcf), suffixes=("_cds", "_variant"), output_type="polars.DataFrame")
        .select(order)
        # polars-bio returns the pairs in no fixed order
        .sort(order, descending=[False, False, True, False])
    )
    row_cds, row_variant = pairs["row_cds"].to_numpy(), pairs["row_variant"].to_numpy()
    # the copies of the rows below take the most memory, so the pairs go first
    del pairs
    cds_rows = cds_df.iloc[row_cds].reset_index(drop=True)
    variant_rows = vcf.drop(columns="Chromosome").iloc[row_variant].reset_index(drop=True)
    return cds_rows.join(variant_rows, rsuffix="_variant")


def place_variants(vcf, chromosomes, reference):
    """
    Finds the equivalent placements of each variant on the given chromosomes, and checks its whole REF against the
    reference genome.

    :param vcf: Variants (DataFrame) with Chromosome, Start, End, Ref and Alt, as ``scan.read_vcf`` returns them
    :param chromosomes: Chromosomes to keep, e.g. those with a coding region
    :param reference: Function that returns the ReferenceSequence of a chromosome
    :return: Tuple (variants, placements). variants holds the kept rows plus variant_row (a key into placements),
             Window_Start and Window_End (the union of the placements) and Ref_matches. Variants whose REF equals
             their ALT are dropped. placements maps variant_row to the list of Placement.
    """

    variants = vcf[vcf["Chromosome"].isin(chromosomes)].reset_index(drop=True)
    variants["variant_row"] = variants.index
    placements = {}
    windows = []
    ref_matches = []
    for row, chromosome, start, end, ref, alt in zip(
        variants.index, variants["Chromosome"], variants["Start"], variants["End"], variants["Ref"], variants["Alt"]
    ):
        chromosome_reference = reference(chromosome)
        placements[row] = variant_placements(int(start), str(ref), str(alt), chromosome_reference)
        windows.append(
            (min(p.start for p in placements[row]), max(p.end for p in placements[row])) if placements[row] else (0, 0)
        )
        ref_matches.append(chromosome_reference.bases(int(start), int(end)) == str(ref).upper())

    variants["Window_Start"] = pd.Series([start for start, _ in windows], index=variants.index, dtype="int64")
    variants["Window_End"] = pd.Series([end for _, end in windows], index=variants.index, dtype="int64")
    variants["Ref_matches"] = pd.Series(ref_matches, index=variants.index, dtype=bool)
    has_placements = pd.Series([bool(placements[row]) for row in variants.index], index=variants.index, dtype=bool)
    return variants[has_placements], placements


def join_variant_windows(cds_df, variants):
    """
    Joins the coding rows with the variant windows that come within 3 bases of them. This finds every coding row
    that a variant can change, also through a splice dinucleotide or an insertion next to it.

    :param cds_df: Coding rows (DataFrame): CDS rows that include the stop codon
    :param variants: Variants with Window_Start and Window_End (place_variants)
    :return: DataFrame with one row per coding row and variant, as ``join_variants_to_cds`` returns it: the coding
             row columns, then the variant columns. Start_variant and End_variant hold the VCF interval of the
             variant.
    """

    coding = cds_df.assign(
        Coding_Start=cds_df["Start"],
        Coding_End=cds_df["End"],
        Start=(cds_df["Start"] - 2).clip(lower=0),
        End=cds_df["End"] + 2,
    )
    windows = variants.assign(
        VCF_Start=variants["Start"],
        VCF_End=variants["End"],
        Start=(variants["Window_Start"] - 1).clip(lower=0),
        End=variants["Window_End"] + 1,
    )
    joined = join_variants_to_cds(coding, windows)
    joined["Start"], joined["End"] = joined.pop("Coding_Start"), joined.pop("Coding_End")
    joined["Start_variant"], joined["End_variant"] = joined.pop("VCF_Start"), joined.pop("VCF_End")
    return joined


def apply_variants(intersection_cds_vcf, placements, cds_df, exons_df, reference):
    """
    Applies each variant to the coding rows of each transcript it joined (see place_in_transcript).

    :param intersection_cds_vcf: Coding rows joined with variants (join_variant_windows)
    :param placements: Placements per variant_row (place_variants)
    :param cds_df: All coding rows; they give the edges of each coding region
    :param exons_df: Exon rows of the annotation (DataFrame with transcript_id, Start, End)
    :param reference: Function that returns the ReferenceSequence of a chromosome
    :return: The rows of the variant-transcript pairs that touch a coding region or its splice dinucleotides, plus
             Exon_Alt_CDS_seq (alt bases of the coding row; None if the alt transcript is unknown),
             UTR5_Ref, UTR5_Alt, UTR3_Ref and UTR3_Alt (the UTR change of the pair, see TranscriptEffect) and
             unknown_reason (None if the alt transcript is known)
    """

    exons_by_transcript = {
        transcript_id: list(zip(group["Start"].astype(int), group["End"].astype(int)))
        for transcript_id, group in exons_df[
            exons_df["transcript_id"].isin(intersection_cds_vcf["transcript_id"])
        ].groupby("transcript_id", observed=True)
    }

    # Start and end of each whole coding region
    coding = cds_df[cds_df["transcript_id"].isin(intersection_cds_vcf["transcript_id"])]
    coding_regions = {
        transcript_id: (int(group["Start"].min()), int(group["End"].max()))
        for transcript_id, group in coding.groupby("transcript_id", observed=True)
    }

    # The loop collects the results per row position in lists. Setting cells of the DataFrame in the loop is slow.
    df = intersection_cds_vcf.reset_index(drop=True)
    starts = df["Start"].to_numpy(dtype="int64")
    ends = df["End"].to_numpy(dtype="int64")
    chromosomes = df["Chromosome"].astype(str).to_numpy()
    strands = df["Strand"].astype(str).to_numpy()
    keep = np.zeros(len(df), dtype=bool)
    alt_coding = [None] * len(df)
    utr5 = [("", "")] * len(df)
    utr3 = [("", "")] * len(df)
    unknown_reasons = [None] * len(df)
    for (transcript_id, variant_row), positions in df.groupby(
        ["transcript_id", "variant_row"], observed=True
    ).indices.items():
        coding_rows = [(int(starts[i]), int(ends[i])) for i in positions]
        effect = place_in_transcript(
            placements[variant_row],
            coding_rows,
            exons_by_transcript.get(transcript_id, []),
            reference(chromosomes[positions[0]]),
            strands[positions[0]],
            coding_regions[transcript_id],
        )
        if effect is None:
            continue
        for i, coding_row in zip(positions, coding_rows):
            keep[i] = True
            if effect.unknown_reason is not None:
                unknown_reasons[i] = effect.unknown_reason
            else:
                alt_coding[i] = effect.alt_coding[coding_row]
                utr5[i] = effect.utr5
                utr3[i] = effect.utr3

    df["Exon_Alt_CDS_seq"] = pd.Series(alt_coding, index=df.index, dtype=object)
    for name, changes in [("UTR5", utr5), ("UTR3", utr3)]:
        df[f"{name}_Ref"] = pd.Series([ref for ref, _ in changes], index=df.index, dtype=object)
        df[f"{name}_Alt"] = pd.Series([alt for _, alt in changes], index=df.index, dtype=object)
    df["unknown_reason"] = pd.Series(unknown_reasons, index=df.index, dtype=object)
    df = df[keep].copy()
    reasons = df.drop_duplicates(["transcript_id", "variant_row"])["unknown_reason"].value_counts()
    if not reasons.empty:
        logger.warning(
            "%d variant-transcript pairs get no prediction because their alt transcript is unknown: %s",
            reasons.sum(),
            ", ".join(f"{count} {reason}" for reason, count in reasons.items()),
        )
    return df


def create_reference_cds(intersection_cds_vcf, cds_df_test):
    """
    Constructs the whole CDS sequence (multiple exons) for transcripts affected by a variant, both in their reference
    and alternative form.
    :param intersection_cds_vcf: DataFrame containing variant-CDS intersection and corresponding alternative CDS sequences
                                 includes: transcript_id, exon_number, Exon_Alt_CDS_seq, and optionally
                                 the UTR change columns and unknown_reason (see apply_variants)
    :param cds_df_test: Reference exon-level CDS data for all transcripts with exon_number
                        includes: transcript_id, exon_number, Start, End, Strand, Exon_CDS_seq, has_stop_codon
    :return: DataFrame with one row per variant-transcript pair, containing full reference and alternative CDS + lengths,
             exon-wise CDS information as tuple (exon number, exon-wise CDS length), has_stop_codon,
             utr5_change and utr3_change (tuples (ref, alt), see TranscriptEffect) and unknown_reason. A pair with
             unknown_reason has None in the alt columns.
    """

    results = []

    # Only transcripts with a variant. transcript_id can be category (see join_variants_to_cds).
    # observed=True: with the pandas 2 default (False), unused categories would be groups too.
    for transcript_id, var_df in intersection_cds_vcf.groupby("transcript_id", observed=True):
        # 1. Get reference exons
        ref_exons = cds_df_test[cds_df_test["transcript_id"] == transcript_id].copy()
        ref_exons = ref_exons.sort_values("Start")

        # Get reference CDS sequence start und stop position for finding position in transcript sequence
        ref_cds_start = ref_exons["Start"].min()
        ref_cds_stop = ref_exons["End"].max()

        # Join reference exon sequences to form full CDS sequence
        ref_seq = "".join(ref_exons["Exon_CDS_seq"].tolist())

        # Collect exon numbers and lengths (for tracking exon contribution later on)
        ref_cds_info = sorted(
            [(row["exon_number"], len(row["Exon_CDS_seq"])) for _, row in ref_exons.iterrows()], key=lambda x: x[0]
        )

        # Get strand info (all should be the same within transcript)
        strand = ref_exons["Strand"].iloc[0]
        ref_seq_final = str(Seq(ref_seq).reverse_complement()) if strand == "-" else ref_seq
        # whether the coding region ends in an annotated stop codon (same within transcript)
        has_stop_codon = bool(ref_exons["has_stop_codon"].iloc[0])

        for variant, cds_df in var_df.groupby(
            ["Chromosome", "Start_variant", "End_variant", "Ref", "Alt"], observed=True
        ):
            # Variant-identifying fields come straight from the group key;
            # ID and gene_id are constant within the group, so read them once.
            chromosome, start_variant, end_variant, ref_allele, alt_allele = variant
            variant_id = cds_df["ID"].iloc[0]
            gene_id = cds_df["gene_id"].iloc[0]

            # An unknown alt transcript gives a row without alt CDS
            unknown_reason = cds_df["unknown_reason"].iloc[0] if "unknown_reason" in cds_df else None
            if isinstance(unknown_reason, str):
                results.append(
                    {
                        "transcript_id": transcript_id,
                        "variant_id": variant_id,
                        "ref_cds_start": ref_cds_start,
                        "ref_cds_stop": ref_cds_stop,
                        "ref_cds_seq": ref_seq_final,
                        "ref_cds_len": len(ref_seq_final),
                        "alt_cds_start": None,
                        "alt_cds_stop": None,
                        "alt_cds_seq": None,
                        "alt_cds_len": None,
                        "chromosome": chromosome,
                        "gene_id": gene_id,
                        "strand": strand,
                        "has_stop_codon": has_stop_codon,
                        "ref": ref_allele,
                        "alt": alt_allele,
                        "start_variant": start_variant,
                        "end_variant": end_variant,
                        "ref_cds_info": ref_cds_info,
                        "alt_cds_info": None,
                        "utr5_change": ("", ""),
                        "utr3_change": ("", ""),
                        "unknown_reason": unknown_reason,
                    }
                )
                continue

            # Sort variant exons
            cds_df = cds_df.sort_values("Start")

            # Copy ref exons for modification
            alt_exons = ref_exons.copy()

            # Replace affected exon sequences with variant versions
            for _, var_row in cds_df.iterrows():
                exon_nr = var_row["exon_number"]
                alt_exons.loc[alt_exons["exon_number"] == exon_nr, "Exon_CDS_seq"] = var_row["Exon_Alt_CDS_seq"]

            # Join and sort alt CDS
            alt_exons = alt_exons.sort_values("Start")

            alt_cds_info = sorted(
                [(row["exon_number"], len(row["Exon_CDS_seq"])) for _, row in alt_exons.iterrows()], key=lambda x: x[0]
            )

            # Get alternative CDS sequence start and stop position for finding position in transcript sequence
            alt_cds_start = alt_exons["Start"].min()
            alt_cds_stop = alt_exons["End"].max()

            alt_seq = "".join(alt_exons["Exon_CDS_seq"].tolist())

            # Apply reverse complement if on minus strand
            alt_seq_final = str(Seq(alt_seq).reverse_complement()) if strand == "-" else alt_seq

            # Append to results
            results.append(
                {
                    "transcript_id": transcript_id,
                    "variant_id": variant_id,
                    "ref_cds_start": ref_cds_start,
                    "ref_cds_stop": ref_cds_stop,
                    "ref_cds_seq": ref_seq_final,
                    "ref_cds_len": len(ref_seq_final),
                    "alt_cds_start": alt_cds_start,
                    "alt_cds_stop": alt_cds_stop,
                    "alt_cds_seq": alt_seq_final,
                    "alt_cds_len": len(alt_seq_final),
                    "chromosome": chromosome,
                    "gene_id": gene_id,
                    "strand": strand,
                    "has_stop_codon": has_stop_codon,
                    "ref": ref_allele,
                    "alt": alt_allele,
                    "start_variant": start_variant,
                    "end_variant": end_variant,
                    "ref_cds_info": ref_cds_info,
                    "alt_cds_info": alt_cds_info,
                    "utr5_change": utr_change(cds_df, "UTR5"),
                    "utr3_change": utr_change(cds_df, "UTR3"),
                    "unknown_reason": None,
                }
            )

    results_df = pd.DataFrame(results)
    # Object dtype keeps the alt positions and lengths as int next to the None of an unknown alt transcript
    for column in ["alt_cds_start", "alt_cds_stop", "alt_cds_len"]:
        if column in results_df:
            results_df[column] = pd.Series([result[column] for result in results], index=results_df.index, dtype=object)
    return results_df


def utr_change(cds_df, utr):
    """
    The UTR change (ref, alt) of a variant-transcript pair, from the columns that apply_variants adds; both empty
    without them.

    :param cds_df: The coding rows of the pair
    :param utr: "UTR5" or "UTR3"
    """

    if f"{utr}_Ref" not in cds_df:
        return "", ""
    return cds_df[f"{utr}_Ref"].iloc[0], cds_df[f"{utr}_Alt"].iloc[0]


def get_transcript_sequence(exons_df, fasta):
    """
    Construct full transcript sequences by concatenating the exon sequences from the FASTA genome reference, grouped by transcript.
    Get transcript length and transcript information as well.
    :param exons_df: DataFrame with the exon rows of the annotation.
                     Must include: transcript_id, strand, chromosome, start, end, exon_number
    :param fasta: Fasta file, reference genome object
    :return: DataFrame with one row per transcript with full transcript sequence, start, end, strand, transcript sequence length, and
             per exon sequence length information for that transcript
    """

    exon_data = []

    # Process each transcript individually
    for transcript_id, group in exons_df.groupby("transcript_id"):
        strand = group.iloc[0]["Strand"]

        if strand not in ["+", "-"]:
            logger.warning("Unknown strand for %s", transcript_id)
            continue

        # Sort by exon start coordinate (strand not considered here yet)
        group_sorted = group.sort_values(by="Start").copy()

        seq_parts = []  # to accumulate exon sequences
        starts = []  # for overall transcript start
        ends = []  # for overall transcript end
        exon_info = []  # for tracking exon_number and length

        # fetch exon sequence and metadata
        for _, row in group_sorted.iterrows():
            chrom = row["Chromosome"]
            start = int(row["Start"])
            end = int(row["End"])
            exon_number = row["exon_number"]

            starts.append(start)
            ends.append(end)

            # Fetch exon sequence from fasta reference genome
            exon_seq = fasta[chrom][start:end]  # .seq
            exon_seq_str = str(exon_seq).upper()
            seq_parts.append(exon_seq_str)

            exon_info.append((exon_number, len(exon_seq_str)))

        # join exon sequences into a full transcript sequence
        joined_seq = "".join(seq_parts)

        # Apply reverse complement for minus strand transcripts
        if strand == "-":
            joined_seq = str(Seq(joined_seq).reverse_complement())
            exon_info = exon_info[::-1]

        exon_data.append(
            {
                "Chromosome": chrom,
                "transcript_id": transcript_id,
                "start": min(starts),
                "end": max(ends),
                "strand": strand,
                "transcript_sequence": joined_seq,
                "transcript_length": len(joined_seq),
                "transcript_exon_info": exon_info,
            }
        )

    exon_seqs = pd.DataFrame(exon_data)
    return exon_seqs


def cds_range_in_transcript(exons, cds):
    """
    Locate the CDS in the transcript sequence, from the exon and CDS coordinates of one transcript.
    Transcript coordinates are 0-based positions in the transcript sequence as built by get_transcript_sequence,
    i.e. read 5' to 3' on both strands. Position 0 is the first base of the 5' exon.

    The start is the transcript position of the 5' CDS base. The end lies one past the 3' CDS base (half-open),
    i.e. after the stop codon if the CDS rows include it. The end is the start plus the summed length of the CDS rows.
    This equals the mapped 3' CDS base plus one, since every CDS row lies inside an exon, also the parts of a stop
    codon split across exons.

    :param exons: Exon rows of one transcript (DataFrame with Start, End, Strand; 0-based half-open genomic coordinates)
    :param cds: CDS rows of the same transcript (DataFrame with Start, End)
    :return: Tuple (start, end) in transcript coordinates, or None if the 5' CDS base lies outside the exons
    """

    strand = exons["Strand"].iloc[0]
    if strand not in ["+", "-"]:
        return None

    # Exons in transcript order, and the genomic position of the 5' CDS base
    exons = exons.sort_values("Start", ascending=(strand == "+"))
    cds_5prime = cds["Start"].min() if strand == "+" else cds["End"].max() - 1

    offset = 0
    for exon_start, exon_end in zip(exons["Start"], exons["End"]):
        if exon_start <= cds_5prime < exon_end:
            start = offset + (cds_5prime - exon_start if strand == "+" else exon_end - 1 - cds_5prime)
            return int(start), int(start + (cds["End"] - cds["Start"]).sum())
        offset += exon_end - exon_start

    return None


def get_exon(cds_pos, exon_info):
    """
    Map a CDS-relative position to the corresponding exon number using exon_info,
    which is a list of (exon_number, exon_length) tuples in CDS order.
    """
    pos_counter = 0
    for exon_number, exon_length in exon_info:
        if cds_pos < pos_counter + exon_length:
            return exon_number
        pos_counter += exon_length
    return exon_info[-1][0]  # fallback


def analyze_sequence(results_df):
    """
    Analyzes reference and alternative CDS for start and stop codons, their positions, and potential premature termination codons (PTCs)

    :param results_df: DataFrame containing CDS sequences and exon information for both reference and alternative sequences, per variant,
                       and has_stop_codon (whether the coding region ends in an annotated stop codon)
    :return: DataFrame with added annotation columns for reference and alternative sequence separately:
             such as start codon position / exon, last codon and its validity as stop codon, first in-frame stop codon + position,
             number and information of all available stop codons, premature stop codon flag.
             The last codon is a valid stop only if it is an annotated stop codon. An in-frame stop codon is premature if it
             lies upstream of the annotated stop codon; without an annotated stop codon, every in-frame stop codon is premature.
    """

    valid_stop_codons = {"TAA", "TAG", "TGA"}
    start_codon = "ATG"

    df = results_df.copy()

    # Initialize result columns for both reference and alternative sequence
    for label in ["ref", "alt"]:
        df[f"{label}_start_codon_pos"] = None
        df[f"{label}_start_codon_exon"] = None  # exon number
        df[f"{label}_last_codon"] = None
        df[f"{label}_valid_stop"] = None
        df[f"{label}_first_stop_codon"] = None
        df[f"{label}_first_stop_pos"] = None
        df[f"{label}_num_stop_codons"] = None
        df[f"{label}_all_stop_codons"] = None
        df[f"{label}_stop_codon_exons"] = None  # exon number
        df[f"{label}_is_premature"] = None

    # Row-wise codon scanning
    for idx, row in df.iterrows():
        has_stop_codon = bool(row["has_stop_codon"])
        for label in ["ref", "alt"]:
            seq = row[f"{label}_cds_seq"]

            exon_info = row[f"{label}_cds_info"]  # for exon number

            # Skip invalid or too-short sequences
            if not isinstance(seq, str) or len(seq) < 3:
                continue

            start_pos = None
            stop_codons = []
            stop_exons = []  # for exon number

            # Scan in codons (step=3)
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i : i + 3]
                if codon == start_codon and start_pos is None:  # first start codon position
                    start_pos = i
                if codon in valid_stop_codons:  # record all stop codons with their positions and exons
                    stop_codons.append((i, codon))
                    stop_exons.append(get_exon(i, exon_info))  # for exon number

            last_codon = seq[-3:]
            # without an annotated stop codon, the last codon is a sense codon or an incomplete one
            is_valid_stop = has_stop_codon and last_codon in valid_stop_codons
            first_stop_pos = stop_codons[0][0] if stop_codons else None
            first_stop = stop_codons[0][1] if stop_codons else None
            # the annotated stop codon is the last codon. Without one, the real stop codon lies downstream of the
            # coding region, so every in-frame stop codon is premature.
            is_premature = first_stop_pos is not None and (not has_stop_codon or first_stop_pos < len(seq) - 3)
            start_exon = get_exon(start_pos, exon_info) if start_pos is not None else None  # for exon number

            # Store results
            df.at[idx, f"{label}_start_codon_pos"] = start_pos
            df.at[idx, f"{label}_start_codon_exon"] = start_exon  # exon number
            df.at[idx, f"{label}_last_codon"] = last_codon
            df.at[idx, f"{label}_valid_stop"] = is_valid_stop
            df.at[idx, f"{label}_first_stop_codon"] = first_stop
            df.at[idx, f"{label}_first_stop_pos"] = first_stop_pos
            df.at[idx, f"{label}_num_stop_codons"] = len(stop_codons)
            df.at[idx, f"{label}_all_stop_codons"] = stop_codons
            df.at[idx, f"{label}_stop_codon_exons"] = stop_exons  # exon number
            df.at[idx, f"{label}_is_premature"] = is_premature

    return df


def start_stop_loss(df):
    """
    Annotates whether a variant caused a start or stop codon loss
    :param df: DataFrame with start & stop codon analysis columns
    :return: Original DataFrames with added columns for "start_loss" and "stop_loss"
    """

    df = df.copy()

    # Start codon loss: the reference CDS starts with a start codon, the alternative CDS does not.
    # A CDS without a leading ATG (e.g. cds_start_NF) has no start codon to lose.
    df["start_loss"] = (df["ref_start_codon_pos"] == 0) & (df["alt_start_codon_pos"] != 0)

    # Stop codon loss: the annotated stop codon no longer encodes a stop in the alternative sequence. A swap to another
    # stop codon, e.g. TAA>TAG, is no loss. Without an annotated stop codon, ref_valid_stop is False: there is no stop
    # codon to lose.
    df["stop_loss"] = (df["ref_valid_stop"] == True) & (df["alt_valid_stop"] != True)

    # Without an alt CDS (unknown alt transcript), start and stop loss are unknown too
    if "unknown_reason" in df:
        unknown = df["unknown_reason"].notna()
        for column in ["start_loss", "stop_loss"]:
            df[column] = df[column].astype(object).where(~unknown, None)

    return df


def splice_alt_cds_into_transcript(row, transcript_seq):
    """
    Splice the alternative CDS sequence into the full transcript sequence to create the alternative transcript.
    A variant can change the UTR next to the coding region, too: utr5_change and utr3_change replace the ref UTR
    bases right before and right after it.
    :param row: A pd.Series row containing "ref_cds_seq" (Reference CDS), "alt_cds_seq" (Alternative / Variant-modified CDS),
                "cds_start_in_transcript" and "cds_end_in_transcript" (from cds_range_in_transcript), and optionally
                utr5_change and utr3_change (tuples (ref, alt) in transcript orientation, see create_reference_cds)
    :param transcript_seq: Full transcript sequence
    :return: Modified (alternative) transcript sequence with the alternative CDS spliced in the correct position,
             or None if the CDS position is unknown, or the transcript does not hold the ref CDS or the ref UTR bases
             there
    """

    ref_cds_seq = row["ref_cds_seq"].upper()
    alt_cds_seq = row["alt_cds_seq"].upper()
    ref_start_idx = row["cds_start_in_transcript"]
    ref_end_idx = row["cds_end_in_transcript"]

    if ref_start_idx is None or transcript_seq[ref_start_idx:ref_end_idx] != ref_cds_seq:
        return None  # Cannot find ref CDS, alignment problem

    utr5_ref, utr5_alt = row.get("utr5_change", ("", ""))
    utr3_ref, utr3_alt = row.get("utr3_change", ("", ""))
    utr5_start = ref_start_idx - len(utr5_ref)
    utr3_end = ref_end_idx + len(utr3_ref)
    if utr5_start < 0 or transcript_seq[utr5_start:ref_start_idx] != utr5_ref:
        return None  # the transcript does not hold the ref 5'UTR bases
    if transcript_seq[ref_end_idx:utr3_end] != utr3_ref:
        return None  # the transcript does not hold the ref 3'UTR bases

    # Replace the reference CDS with the variant-modified / alternative one
    new_transcript_seq = (
        transcript_seq[:utr5_start]
        + utr5_alt
        + alt_cds_seq
        + utr3_alt
        + transcript_seq[utr3_end:]
    )  # fmt: skip

    return new_transcript_seq


def in_frame_codons(seq, start, codons):
    """
    Yield (position, codon) for each codon of the sequence that is one of codons, reading in frame from position start
    to the end of the sequence.
    """

    for i in range(start, len(seq) - 2, 3):
        if seq[i : i + 3] in codons:
            yield i, seq[i : i + 3]


def first_stop_codon(seq, start):
    """Return the position of the first stop codon of the sequence, reading in frame from position start, or None."""

    return next((i for i, _ in in_frame_codons(seq, start, {"TAA", "TAG", "TGA"})), None)


def _common_prefix_length(a, b):
    """Return the length of the longest common prefix of the strings a and b."""

    low, high = 0, min(len(a), len(b))
    while low < high:
        middle = (low + high + 1) // 2
        if a[:middle] == b[:middle]:
            low = middle
        else:
            high = middle - 1
    return low


def annotated_stop_in_alt(ref_seq, alt_seq, stop):
    """
    Return the positions of the annotated stop codon in the alternative sequence, as a tuple of two. A variant upstream
    of the stop codon shifts it by its net indel length. A variant that starts in the stop codon or downstream of it
    does not shift it.

    The variant replaces the part between the longest common prefix and the longest common suffix of both sequences.
    So the result depends on the alternative sequence only, not on how the VCF describes the variant. An indel in a
    repeat has more than one such placement. Taking the prefix first places it 3'-most, as HGVS does; this gives the
    second position. Taking the suffix first places it 5'-most. If that placement ends at or before the stop codon,
    the first position is the shifted stop codon. Otherwise both positions are the same.
    The two differ if the repeat reaches into the stop codon. E.g. TCC inserted right before the stop codon TAA can
    also be placed after its T, inside it. Then the first position is 3 nt downstream of the annotated one, and the
    second is the annotated one. A TAA inserted right before the stop codon TAA can also be placed right after it.
    A variant can replace the first base of the stop codon together with bases upstream of it. If the replacement
    reaches into the stop codon in the alternative sequence, the stop codon maps to the same offset in the replacement
    as in the reference. Otherwise it maps to the position anchored on the unchanged sequence downstream, i.e. shifted
    by the net indel length, so it does not depend on the placement.
    This counts the annotated stop codon as lost only if every representation of the variant removes it. A delins is
    not split into an SNV plus an indel. E.g. ATAG>T at the last sense codon and the stop codon in GTA TAG TAG CAT
    gives GTT TAG CAT: a synonymous change plus the deletion of a stop codon in a stop codon repeat. The protein is
    unchanged and the 3'UTR is 3 nt shorter, but the result is a stop loss.

    :param ref_seq: Reference sequence, e.g. the transcript
    :param alt_seq: Alternative sequence, the reference sequence with the variant applied
    :param stop: Position of the first base of the annotated stop codon in ref_seq
    :return: Tuple of the two positions of the annotated stop codon in alt_seq: from the 5'-most placement if it ends
             at or before the stop codon, and from the 3'-most placement
    """

    net_length = len(alt_seq) - len(ref_seq)
    prefix = _common_prefix_length(ref_seq, alt_seq)
    suffix = _common_prefix_length(ref_seq[::-1], alt_seq[::-1])

    # 3'-most placement: the variant replaces ref_seq[prefix:ref_end] by alt_seq[prefix:alt_end]
    suffix_3prime = min(suffix, min(len(ref_seq), len(alt_seq)) - prefix)
    ref_end, alt_end = len(ref_seq) - suffix_3prime, len(alt_seq) - suffix_3prime
    if ref_end <= stop or alt_end <= stop:
        position_3prime = stop + net_length
    else:
        position_3prime = stop

    # 5'-most placement: the variant ends where the longest common suffix starts
    if len(ref_seq) - suffix <= stop:
        return stop + net_length, position_3prime
    return position_3prime, position_3prime


def ends_at_annotated_stop(row):
    """
    Return whether the reference transcript, read in frame from the annotated start codon, reads its first stop codon
    at the annotated stop codon.

    Not so if the annotated stop codon is out of the reading frame, as in a cds_start_NF CDS read from its first base,
    or if the transcript has no stop codon there. Not so either if an in-frame stop codon lies upstream of it: a
    selenocysteine TGA or a misannotation. The annotation marks selenocysteine codons, but the pipeline does not read
    them, so it cannot tell the two apart.

    :param row: A pd.Series row containing transcript_seq, cds_start_in_transcript and cds_end_in_transcript (from
                cds_range_in_transcript)
    """

    first_stop = first_stop_codon(row["transcript_seq"], row["cds_start_in_transcript"])
    return first_stop == row["cds_end_in_transcript"] - 3


def annotated_stop_distance(row, first_stop):
    """
    Return how far a stop codon of the alternative transcript lies upstream of the annotated stop codon, in nt.
    Positive means upstream, i.e. a PTC. 0 means the stop codon is the annotated one, and negative means it lies
    downstream.

    The stop codon is the annotated one if it lies at one of the two positions from annotated_stop_in_alt. Otherwise
    the distance is to the first of them. After an in-frame indel right before the stop codon, the first position is
    the shifted stop codon, also if the indel can be placed inside it.

    :param row: A pd.Series row containing transcript_seq, alt_transcript_seq and cds_end_in_transcript (from
                cds_range_in_transcript)
    :param first_stop: Position of the stop codon in the alternative transcript, or None
    :return: The distance in nt, or None if first_stop is None
    """

    if first_stop is None:
        return None
    positions = annotated_stop_in_alt(
        row["transcript_seq"], row["alt_transcript_seq"], row["cds_end_in_transcript"] - 3
    )
    return 0 if first_stop in positions else positions[0] - first_stop


def classify_first_stop(row, first_stop):
    """
    Classify the first in-frame stop codon of the alternative transcript, read from the annotated start codon.

    With an annotated stop codon, compare the first stop with it in alternative transcript coordinates (see
    annotated_stop_distance). A stop upstream of it is premature. A stop at its position is neither premature nor a
    stop loss. A stop downstream of it, or no stop up to the transcript end (nonstop), is a stop loss.
    Without an annotated stop codon (has_stop_codon False), there is no position to compare with. A stop inside the
    alternative CDS is premature, and a stop past its end is neither.

    :param row: A pd.Series row containing has_stop_codon, transcript_seq, alt_transcript_seq, alt_cds_seq,
                alt_cds_start_in_transcript, and cds_end_in_transcript (from cds_range_in_transcript)
    :param first_stop: Position of the first in-frame stop codon in the alternative transcript, or None
    :return: Tuple (alt_is_premature, stop_loss)
    """

    if not row["has_stop_codon"]:
        alt_cds_end = row["alt_cds_start_in_transcript"] + len(row["alt_cds_seq"])
        return first_stop is not None and first_stop + 3 <= alt_cds_end, False

    distance = annotated_stop_distance(row, first_stop)
    return distance is not None and distance > 0, distance is None or distance < 0


def analyze_transcript(results_df):
    """
    Analyze the alternative transcript sequence: classify its first in-frame stop codon, and in cases of start or stop
    codon loss, scan for new in-frame start or stop codons.

    The classification reads the alternative transcript in frame from the annotated start codon, at
    alt_cds_start_in_transcript, through the CDS into the 3'UTR (see classify_first_stop). It replaces alt_is_premature and stop_loss from analyze_sequence and
    start_stop_loss, which see only the CDS: after a frameshift, the last codon of the alternative CDS is out of frame.
    The comparison with the annotated stop codon needs a reference transcript that, read the same way, stops there
    (see ends_at_annotated_stop). Rows where it does not, and rows without an alternative transcript, keep those flags.

    :param results_df: DataFrame containing transcript sequence data and annotations, including start_loss and stop_loss flags,
                       alt_cds_start_in_transcript, and the columns that classify_first_stop reads
    :return: pandas DataFrame with additional columns for rescued start / stop codon information
    """

    valid_stop_codons = {"TAA", "TAG", "TGA"}
    start_codon = "ATG"

    df = results_df.copy()

    # Add new columns to store results
    df["transcript_start_codon_pos"] = None
    df["transcript_start_codon_exon"] = None  # for exon number
    df["transcript_last_codon"] = None
    df["transcript_valid_stop"] = None
    df["transcript_first_stop_codon"] = None
    df["transcript_first_stop_pos"] = None
    df["transcript_num_stop_codons"] = None
    df["transcript_all_stop_codons"] = None
    df["transcript_stop_codon_exons"] = None  # for exon number

    for idx, row in df.iterrows():
        seq = row["alt_transcript_seq"]
        cds_start = row["alt_cds_start_in_transcript"]

        exon_info = row["transcript_exon_info"]  # for exon number

        # Skip rows with invalid or too-short sequences, or without a CDS position in the transcript
        if not isinstance(seq, str) or len(seq) < 3 or pd.isna(cds_start):
            continue

        # Read codons in frame from the annotated start codon to the end of the transcript
        scan_start = cds_start
        stop_codons_in_frame = list(in_frame_codons(seq, scan_start, valid_stop_codons))

        # The comparison needs a reference transcript that reads its first stop codon at the annotated one. Otherwise,
        # e.g. for a selenocysteine TGA or an annotated stop codon out of frame, the row keeps the flags from the CDS.
        stop_loss = row["stop_loss"]
        if not row["has_stop_codon"] or ends_at_annotated_stop(row):
            first_stop = stop_codons_in_frame[0][0] if stop_codons_in_frame else None
            is_premature, stop_loss = classify_first_stop(row, first_stop)
            df.at[idx, "alt_is_premature"] = is_premature
            df.at[idx, "stop_loss"] = stop_loss

        start_pos = None
        start_exon = None  # for exon number
        stop_codons = []

        # only analyze rows flagged with start or stop codon loss: skip the others and fill with None values
        # Skips only if we have both start_loss = FALSE and stop_loss = FALSE. If one is true, then don't skip.
        if not (row["start_loss"] or stop_loss):
            continue

        # START LOSS rescue search
        if row["start_loss"]:
            # Walk through sequence starting at CDS start with +1 positions until start codon is found
            for i in range(cds_start, len(seq) - 2):
                codon = seq[i : i + 3]
                if codon == start_codon:
                    start_pos = i

                    start_exon = get_exon(start_pos, exon_info)  # for exon number

                    # From new start codon, scan codons in frame
                    stop_codons = list(in_frame_codons(seq, i, valid_stop_codons))

                    break

        # STOP LOSS readthrough: the in-frame stop codons from the annotated start codon on. Only without a start
        # loss: after a start loss, the scan above already reads on to the transcript end (3'UTR).
        else:
            start_pos = next((i for i, _ in in_frame_codons(seq, scan_start, {start_codon})), None)
            start_exon = get_exon(start_pos, exon_info) if start_pos is not None else None  # for exon number
            stop_codons = stop_codons_in_frame

        stop_exons = [get_exon(i, exon_info) for i, _ in stop_codons]  # for exon number

        # Set last codon in transcript (use last full codon)
        last_codon = seq[-3:] if len(seq) >= 3 else None
        is_valid_stop = last_codon in valid_stop_codons
        first_stop_pos = stop_codons[0][0] if stop_codons else None
        first_stop = stop_codons[0][1] if stop_codons else None

        # Store results
        df.at[idx, "transcript_start_codon_pos"] = start_pos
        df.at[idx, "transcript_start_codon_exon"] = start_exon  # for exon number
        df.at[idx, "transcript_last_codon"] = last_codon
        df.at[idx, "transcript_valid_stop"] = is_valid_stop
        df.at[idx, "transcript_first_stop_codon"] = first_stop
        df.at[idx, "transcript_first_stop_pos"] = first_stop_pos
        df.at[idx, "transcript_num_stop_codons"] = len(stop_codons)
        df.at[idx, "transcript_all_stop_codons"] = stop_codons
        df.at[idx, "transcript_stop_codon_exons"] = stop_exons  # for exon number

    return df
