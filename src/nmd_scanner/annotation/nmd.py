import logging
import pandas as pd
import pyranges as pr
import numpy as np
from tqdm import tqdm
from ..io import extract_fasta_sequence
from .models import CDSAnnotation, NMDPrediction, TranscriptFeatures
from . import utils
from .regressors import nmdeff_predict

logger = logging.getLogger(__name__)

__all__ = ["annotate_cds", "evaluate_nmd_escape"]


def annotate_cds(cds_df, vcf, fasta, exons_df) -> list[CDSAnnotation]:

    """
    Main function for extracting reference coding sequence, alternative coding sequence by incorporating the variant, 
    analyzing for premature termination codons (PTCs), start and stop loss, and getting the transcript information.

    :param cds_df: CDS entries from the GTF file (Dataframe)
    :param vcf: Parsed VCF variant entries (PyRanges object)
    :param fasta: Reference genome sequence (pyfaidx.Fasta object)
    :param exons_df: All exonic entries from the GTF file (DataFrame)
    :return: List of CDSAnnotation objects with core CDS and codon analysis
    """

    # Adjust the last 3 CDS positions to include stop codons
    cds_df_adj = utils.adjust_last_cds_for_stop_codon(cds_df)
    logger.info("Adjusted last CDS for stop codon")

    # adjust exon_number as int datatype
    cds_df_adj["exon_number"] = cds_df_adj["exon_number"].astype(int)

    # Intersect variants with CDS regions
    intersection_cds_vcf = pr.PyRanges(cds_df_adj).join(vcf, how=None, suffix="_variant").df
    logger.info("Joined variants with CDS entries")

    ##########################################################################################
    # TODO: fix minus strand variants (only for TCGA and MMRF VCF!)
    # Fix REF and ALT for minus-strand CDSs
    #mask_minus_strand = intersection_cds_vcf["Strand"] == "-"
    #intersection_cds_vcf.loc[mask_minus_strand, "Ref"] = intersection_cds_vcf.loc[mask_minus_strand, "Ref"].apply(
    #   lambda seq: str(Seq(seq).reverse_complement()))
    #intersection_cds_vcf.loc[mask_minus_strand, "Alt"] = intersection_cds_vcf.loc[mask_minus_strand, "Alt"].apply(
    #   lambda seq: str(Seq(seq).reverse_complement()))
    ##########################################################################################

    # intersection = intersection_test.copy()
    # df3["Exon_CDS_seq"] = df3.apply(lambda row: fasta[row["Chromosome"]][row["Start"]:row["End"]].seq.upper(), axis=1)
    #intersection_cds_vcf["Exon_CDS_seq"] = [
    #    fasta[chrom][start:end].seq.upper()
    #    for chrom, start, end in zip(intersection_cds_vcf["Chromosome"], intersection_cds_vcf["Start"], intersection_cds_vcf["End"])
    #]

    # Fetch reference CDS sequence for each variant region
    logger.info("Creating exon CDS sequence")
    intersection_cds_vcf['Exon_CDS_seq'] = intersection_cds_vcf.apply(lambda row: extract_fasta_sequence(
        row["Chromosome"], row["Start"], row["End"], fasta), axis=1)
    logger.info("Created exon CDS sequence")

    # Apply variant to CDS and compute alternative CDS sequence and lengths
    logger.info("Applying variants to CDS regions")
    tqdm.pandas(desc="Processing variants")
    with tqdm(total=len(intersection_cds_vcf), desc="Applying variants to CDS") as pbar:
        results = []
        for idx, row in intersection_cds_vcf.iterrows():
            results.append(utils.apply_variant_edge_aware_with_lengths(row))
            pbar.update(1)
        intersection_cds_vcf[["Exon_CDS_length", "Exon_Alt_CDS_seq", "Exon_Alt_CDS_length"]] = pd.DataFrame(results, index=intersection_cds_vcf.index)
    logger.info("Created exon CDS and alt CDS sequence")

    ##### New ######
    # Filter out Variants with a reference mismatch
    mismatched_rows = intersection_cds_vcf[intersection_cds_vcf["Exon_Alt_CDS_seq"].isna()]
    logger.warning(f"Skipping {len(mismatched_rows)} variants due to reference mismatches")
    if not mismatched_rows.empty:
        logger.debug(f"Mismatched variants:\n{mismatched_rows[['transcript_id', 'Chromosome', 'Start_variant', 'End_variant', 'Ref', 'Alt']].to_string(index=False)}")
    intersection_cds_vcf = intersection_cds_vcf[intersection_cds_vcf["Exon_Alt_CDS_seq"].notna()].copy()
    ################

    # Save exon-variant merge result
    #output_path = os.path.join(output, "1_variant_exon_output.tsv")
    #intersection_cds_vcf.to_csv(output_path, sep="\t", index=False)
    #print(f"Creating {output_path}: done.")

    # Limit to relevant transcript (to save time)
    relevant_transcripts = intersection_cds_vcf["transcript_id"].unique()
    cds_df_adj = cds_df_adj[cds_df_adj["transcript_id"].isin(relevant_transcripts)].copy()

    # cds_df_adj["Exon_CDS_seq"] = cds_df_adj.apply(lambda row: fasta[row["Chromosome"]][row["Start"]:row["End"]].seq.upper(), axis=1)
    #cds_df_adj["Exon_CDS_seq"] = [
    #    fasta[chrom][start:end].seq.upper()
    #    for chrom, start, end in zip(cds_df_adj["Chromosome"], cds_df_adj["Start"], cds_df_adj["End"])
    #]

    # Fetch reference sequence for all CDS entries per (relevant) transcripts
    cds_df_adj['Exon_CDS_seq'] = cds_df_adj.apply(lambda row: extract_fasta_sequence(
        row["Chromosome"], row["Start"], row["End"], fasta), axis=1)

    # get full reference CDS per transcript by stiching exon CDS regions, plus alternative CDS with CDS exon information for both ref and alt
    results_df = utils.create_reference_cds(intersection_cds_vcf, cds_df_adj)
    logger.info("Created reference CDS")

    # make intermediate output file of results_df
    #output_path = os.path.join(output, "3_create_reference_CDS.tsv")
    #results_df.to_csv(output_path, sep="\t", index=False)
    #print(f"Creating {output_path}: done.")

    # Get transcript sequence for relevant transcripts (speed up process) + length and transcript exon information (Tuple: exon number & exon length)
    exons_df = exons_df[exons_df["transcript_id"].isin(relevant_transcripts)].copy()
    exon_seqs = utils.get_transcript_sequence(exons_df, fasta)
    logger.info("Retrieved transcript sequences")

    # make output file of exon_seqs
    #output_path = os.path.join(output, "4_transcript_sequences.tsv")
    #exon_seqs.to_csv(output_path, sep="\t", index=False)
    #print(f"Creating {output_path}: done.")

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

    logger.info("Validating CDS in transcript sequences")
    with tqdm(total=len(results_df), desc="Validating CDS") as pbar:
        cds_checks = []
        for idx, row in results_df.iterrows():
            cds_checks.append(check_cds_in_transcript(row))
            pbar.update(1)
        results_df["cds_in_transcript"] = cds_checks

    # TODO: Analyze reference and alternative CDS for start / stop codons
    analysis_df = utils.analyze_sequence(results_df)
    loss_df = utils.start_stop_loss(analysis_df) # instead of loss_analysis_df (test for start stop loss)
    logger.info("Analyzed sequences")


    # Annotate transcript information (transcript start, end, sequence, length, exon info) in case of start or stop loss
    # transcript sequences are in: exon_seqs_subset
    logger.info("Annotating transcript information for start/stop loss")
    transcript_starts = exon_seqs.set_index("transcript_id")["start"].to_dict()
    loss_df["transcript_start"] = loss_df["transcript_id"].map(transcript_starts)
    transcript_ends = exon_seqs.set_index("transcript_id")["end"].to_dict()
    loss_df["transcript_end"] = loss_df["transcript_id"].map(transcript_ends)
    transcript_sequences = exon_seqs.set_index("transcript_id")["transcript_sequence"].to_dict()  # create map of transcript-id to transcript sequence
    loss_df["transcript_seq"] = loss_df["transcript_id"].map(transcript_sequences)
    transcript_lengths = exon_seqs.set_index("transcript_id")["transcript_length"].to_dict()
    loss_df["transcript_length"] = loss_df["transcript_id"].map(transcript_lengths)

    # In case of start or stop loss:
    # Splice alternative CDS into reference transcript sequence to create alternative transcript sequence and measure new length
    logger.info("Splicing alternative CDS into transcripts")
    with tqdm(total=len(loss_df), desc="Splicing alt CDS") as pbar:
        alt_seqs = []
        for idx, row in loss_df.iterrows():
            alt_seqs.append(
                utils.splice_alt_cds_into_transcript(row, row["transcript_seq"])
                if pd.notnull(row["transcript_seq"]) else None
            )
            pbar.update(1)
        loss_df["alt_transcript_seq"] = alt_seqs
    loss_df["alt_transcript_length"] = loss_df["alt_transcript_seq"].apply(
        lambda x: len(x) if pd.notnull(x) else None
    )


    # Add exon information to dataframe
    transcript_exon_info = exon_seqs.set_index("transcript_id")["transcript_exon_info"].to_dict()
    loss_df["transcript_exon_info"] = loss_df["transcript_id"].map(transcript_exon_info)


    # Analyze transcript sequence (e.g., frame, length, stop codon position, etc.) in case of start or stop loss
    analyze_transcript_df = utils.analyze_transcript(loss_df)

    # Convert DataFrame to list of CDSAnnotation objects
    cds_annotations = []
    for _, row in analyze_transcript_df.iterrows():
        cds_annotation = CDSAnnotation(
            transcript_id=row.get('transcript_id') or row.get('TRANSCRIPT_ID'),
            variant_id=row.get('variant_id') or row.get('VARIANT_ID', ''),
            ref_cds_start=row.get('ref_cds_start') or row.get('REF_CDS_START'),
            ref_cds_stop=row.get('ref_cds_stop') or row.get('REF_CDS_STOP'),
            ref_cds_seq=row.get('ref_cds_seq') or row.get('REF_CDS_SEQ'),
            ref_cds_len=row.get('ref_cds_len') or row.get('REF_CDS_LEN'),
            alt_cds_start=row.get('alt_cds_start') or row.get('ALT_CDS_START'),
            alt_cds_stop=row.get('alt_cds_stop') or row.get('ALT_CDS_STOP'),
            alt_cds_seq=row.get('alt_cds_seq') or row.get('ALT_CDS_SEQ'),
            alt_cds_len=row.get('alt_cds_len') or row.get('ALT_CDS_LEN'),
            chromosome=row.get('Chromosome') or row.get('chromosome'),
            gene_id=row.get('gene_id') or row.get('GENE_ID'),
            strand=row.get('Strand') or row.get('strand'),
            ref=row.get('Ref') or row.get('ref'),
            alt=row.get('Alt') or row.get('alt'),
            start_variant=row.get('Start_variant') or row.get('start_variant'),
            end_variant=row.get('End_variant') or row.get('end_variant'),
            ref_cds_info=row.get('ref_cds_info'),
            alt_cds_info=row.get('alt_cds_info'),
            cds_in_transcript=row.get('cds_in_transcript'),
            ref_start_codon_pos=row.get('ref_start_codon_pos'),
            ref_start_codon_exon=row.get('ref_start_codon_exon'),
            ref_last_codon=row.get('ref_last_codon'),
            ref_valid_stop=row.get('ref_valid_stop'),
            ref_first_stop_codon=row.get('ref_first_stop_codon'),
            ref_first_stop_pos=row.get('ref_first_stop_pos'),
            ref_num_stop_codons=row.get('ref_num_stop_codons'),
            ref_all_stop_codons=row.get('ref_all_stop_codons'),
            ref_stop_codon_exons=row.get('ref_stop_codon_exons'),
            ref_is_premature=row.get('ref_is_premature'),
            alt_start_codon_pos=row.get('alt_start_codon_pos'),
            alt_start_codon_exon=row.get('alt_start_codon_exon'),
            alt_last_codon=row.get('alt_last_codon'),
            alt_valid_stop=row.get('alt_valid_stop'),
            alt_first_stop_codon=row.get('alt_first_stop_codon'),
            alt_first_stop_pos=row.get('alt_first_stop_pos'),
            alt_num_stop_codons=row.get('alt_num_stop_codons'),
            alt_all_stop_codons=row.get('alt_all_stop_codons'),
            alt_stop_codon_exons=row.get('alt_stop_codon_exons'),
            alt_is_premature=row.get('alt_is_premature'),
            start_loss=row.get('start_loss'),
            stop_loss=row.get('stop_loss'),
            transcript_start=row.get('transcript_start'),
            transcript_end=row.get('transcript_end'),
            transcript_seq=row.get('transcript_seq'),
            transcript_length=row.get('transcript_length'),
            alt_transcript_seq=row.get('alt_transcript_seq'),
            alt_transcript_length=row.get('alt_transcript_length'),
            transcript_exon_info=row.get('transcript_exon_info'),
            transcript_start_codon_pos=row.get('transcript_start_codon_pos'),
            transcript_start_codon_exon=row.get('transcript_start_codon_exon'),
            transcript_last_codon=row.get('transcript_last_codon'),
            transcript_valid_stop=row.get('transcript_valid_stop'),
            transcript_first_stop_codon=row.get('transcript_first_stop_codon'),
            transcript_first_stop_pos=row.get('transcript_first_stop_pos'),
            transcript_num_stop_codons=row.get('transcript_num_stop_codons'),
            transcript_all_stop_codons=row.get('transcript_all_stop_codons'),
            transcript_stop_codon_exons=row.get('transcript_stop_codon_exons')
        )
        cds_annotations.append(cds_annotation)
    
    return cds_annotations


def evaluate_nmd_escape(cds: CDSAnnotation, features: TranscriptFeatures) -> NMDPrediction:
    """
    Evaluate whether a premature stop codon is likely to escape nonsense-mediated decay (NMD) 
    based on established biological rules.
    
    This function applies five NMD escape rules:
    1. Last exon rule: The PTC is in the last exon
    2. 50nt penultimate rule: The PTC is within 50 nucleotides upstream of the last exon junction
    3. Long exon rule: The PTC is in an exon with >407 nucleotides
    4. Start proximal rule: The PTC is within 150 nucleotides of the start codon
    5. Single exon rule: The transcript consists only of a single exon
    
    A PTC is considered to escape NMD if it satisfies any of the above rules.
    
    :param cds: CDSAnnotation object with CDS and codon information
    :param features: TranscriptFeatures object with computed transcript-level features
    :return: NMDPrediction object with rule evaluation results
    """
    # Only relevant for premature stop codons
    if not cds.alt_is_premature:
        return NMDPrediction(
            nmd_last_exon_rule=False,
            nmd_50nt_penultimate_rule=False,
            nmd_long_exon_rule=False,
            nmd_start_proximal_rule=False,
            nmd_single_exon_rule=False,
            nmd_escape=False,
            nmd_efficiency=None
        )
    
    # Extract relevant data from CDSAnnotation object
    stop_pos = cds.alt_first_stop_pos
    exon_info = cds.transcript_exon_info or []
    start_pos = cds.alt_start_codon_pos
    
    # Extract data from TranscriptFeatures
    total_exon_count = features.total_exon_count
    downstream_exon_count = features.downstream_exon_count
    ptc_exon_length = features.ptc_exon_length
    
    # Preprocess exon info
    sorted_exons = sorted(exon_info, key=lambda x: x[0])
    exon_offsets = {}
    offset = 0
    for exon_num, length in sorted_exons:
        exon_offsets[exon_num] = (offset, offset + length)
        offset += length
    
    # Single exon rule
    rule_single_exon = total_exon_count == 1
    
    # Last exon rule
    rule_last_exon = downstream_exon_count == 0 if downstream_exon_count is not None else False
    
    # 50nt from penultimate exon end
    if len(sorted_exons) >= 2:
        penultimate_exon_num, penultimate_len = sorted_exons[-2]
        pen_start, pen_end = exon_offsets.get(penultimate_exon_num, (None, None))
        rule_50nt_penultimate = pen_end is not None and (stop_pos >= pen_end - 50) and (stop_pos < pen_end)
    else:
        rule_50nt_penultimate = False
    
    # Long exon rule (with exon longer than >407nt)
    rule_long_exon = ptc_exon_length is not None and ptc_exon_length > 407
    
    # Start-proximal rule (closer than 150nt from the start codon)
    rule_start_proximal = start_pos is not None and stop_pos is not None and (stop_pos - start_pos) < 150 and (stop_pos - start_pos) >= 0
    
    # NMD escape if any rule is true
    escape = rule_last_exon or rule_50nt_penultimate or rule_long_exon or rule_start_proximal or rule_single_exon

    nmd_efficiency = nmdeff_predict(
        np.array([int(cds.start_loss or 0)]),
        np.array([int(cds.stop_loss or 0)]),
        np.array([total_exon_count or 0]),
        np.array([int(features.ptc_less_than_150nt_to_start or 0)]),
        np.array([int(rule_long_exon)]),
        np.array([int(rule_start_proximal)]),
        np.array([int(rule_single_exon)]),
        np.array([int(escape)]),
        np.array([downstream_exon_count or 0]),
        np.array([int(rule_last_exon)]),
        np.array([features.ptc_to_start_codon or 0]),
        np.array([features.stop_codon_distance or 0]),
        np.array([ptc_exon_length or 0]),
        np.array([features.ptc_to_intron or 0]),
        np.array([features.upstream_exon_count or 0]),
        np.array([int(rule_50nt_penultimate)]),
        np.array([features.utr5_length or 0]),
        np.array([features.utr3_length or 0]),
        np.array([cds.transcript_length or 0])
    )[0]  # Extract single prediction value
    
    return NMDPrediction(
        nmd_last_exon_rule=rule_last_exon,
        nmd_50nt_penultimate_rule=rule_50nt_penultimate,
        nmd_long_exon_rule=rule_long_exon,
        nmd_start_proximal_rule=rule_start_proximal,
        nmd_single_exon_rule=rule_single_exon,
        nmd_escape=escape,
        nmd_efficiency= nmd_efficiency
    )
