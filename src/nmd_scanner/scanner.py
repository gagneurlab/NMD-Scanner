# Import dependencies
import logging
import os
import pandas as pd
from dataclasses import fields
from tqdm import tqdm
from nmd_scanner.io import read_vcf, read_gtf, read_fasta
from .annotation import calculate_transcript_features, evaluate_nmd_escape, annotate_cds, NMDResult, TranscriptFeatures, NMDPrediction

logger = logging.getLogger(__name__)


def annotate_nmd_pandas(vcf_path, gtf_path, fasta_path, output=None, reassign_exons=False, canonical_only=False, output_format=None) -> pd.DataFrame:
    """
    Main function for NMD scanner that returns a DataFrame instead of a list of NMDResult objects.
    This is an alternative to the annotate_nmd function and can be used if you prefer working with DataFrames directly.
    """
    results = annotate_nmd(vcf_path, gtf_path, fasta_path,
                           output, reassign_exons, canonical_only, output_format)
    return NMDResult.to_dataframe(results)


def annotate_nmd(vcf_path, gtf_path, fasta_path, output=None, reassign_exons=False, canonical_only=False, output_format=None) -> list[NMDResult]:
    """
    Main function for NMD scanner

    Steps:
    1. Read input files (VCF, GTF, FASTA)
    2. Assign exon numbers (optional, recommended for hg19)
    3. Parse and preprocess gene annotations (CDS, exons)
    4. Extract premature termination codons (PTCs) & Evaluate NMD escape rules
    5. Return minimal results (variant_id, transcript_id, gene_id, nmd_escape) or detailed if specified
    6. Save output results

    :param vcf_path: path to the input VCF file
    :param gtf_path: path to the input GTF annotation file
    :param fasta_path: path to the reference FASTA file
    :param output: directory to save the results in
    :param reassign_exons: whether to reassign exon numbers
    :param canonical_only: if True, filter to canonical transcripts only
    :param output_format: output format ('csv' or 'parquet'). If None, inferred from file extension
    :return: DataFrame with NMD results
    """

    # read VCF file (variants)
    logger.info(f"Reading VCF file: {vcf_path}")
    vcf = read_vcf(vcf_path)
    logger.info(f"VCF shape: {vcf.df.shape}")

    # read GTF file (gene annotation)
    logger.info(f"Reading GTF file: {gtf_path}")
    gtf = read_gtf(gtf_path, reassign_exons=reassign_exons,
                   canonical_only=canonical_only)
    logger.info(f"GTF File shape: {gtf.df.shape}")

    # read FASTA file (genome sequence)
    logger.info(f"Reading FASTA file: {fasta_path}")
    fasta = read_fasta(fasta_path)
    logger.info(f"FASTA file loaded with {len(fasta)} sequences")

    # extract CDS regions from the GTF file
    cds = gtf[gtf.Feature == "CDS"]
    cds_df = cds.df

    # extract exon regions from the GTF file and compute exon related metrics:
    # exon length & number of exons contained in each transcript
    exons = gtf[gtf.Feature == "exon"]
    exons_df = exons.df
    exons_df["exon_length"] = exons_df["End"] - exons_df["Start"]

    # Step 1: Annotate CDS
    logger.info("Annotating CDS")
    cds_annotations = annotate_cds(cds_df, vcf, fasta, exons_df)

    # Step 2: Calculate transcript features and evaluate NMD for each
    logger.info(
        "Calculating transcript features and evaluating NMD escape rules")
    results = []
    for cds in tqdm(cds_annotations, desc="Processing variants"):
        # Calculate transcript features
        features = calculate_transcript_features(cds)

        # Evaluate NMD escape
        prediction = evaluate_nmd_escape(cds, features)
        results.append(NMDResult.from_annotation(cds, features, prediction))

    # Write output if output directory is specified
    if output is not None:
        results_df = NMDResult.to_dataframe(results)

        # only keep cases with a premature termination codon (PTC) for output
        results_df = results_df[results_df['alt_is_premature']]

        # in case of missanotaed transcripts (presence of reference PTC), we want to keep only those with an alternative PTC that comes before the "reference PTC"
        results_df = results_df[~((results_df['likely_misannotated']) & (
            results_df['alt_first_stop_pos'] >= results_df['ref_first_stop_pos']))]

        # only write a subset of columns for simplified output
        columns_to_keep = ['variant_id', 'transcript_id', 'gene_id', 'chromosome',
                           'start_variant', 'end_variant', 'ref', 'alt', 'alt_is_premature',
                           *[f.name for f in fields(TranscriptFeatures)],
                           *[f.name for f in fields(NMDPrediction)]]
        results_df = results_df[columns_to_keep]

        # Determine output file path and format
        vcf_base = os.path.splitext(os.path.basename(vcf_path))[0]

        # Infer format from extension if not specified
        if output_format is None:
            if output.endswith('.parquet'):
                output_format = 'parquet'
            elif output.endswith('.csv'):
                output_format = 'csv'
            else:
                output_format = 'csv'  # Default

        if output_format == 'parquet':
            output_file = output if output.endswith('.parquet') else os.path.join(
                output, f"{vcf_base}_nmdscanner.parquet")
            results_df.to_parquet(output_file, index=False)
        else:  # CSV
            output_file = output if output.endswith('.csv') else os.path.join(
                output, f"{vcf_base}_nmdscanner.csv")
            results_df.to_csv(output_file, index=False)

        logger.info(f"Results written to {output_file}")

    return results
