# Import dependencies
import logging
import os
import pandas as pd
from dataclasses import asdict
from tqdm import tqdm
from nmd_scanner.io import read_vcf, read_gtf, read_fasta
from .annotation import calculate_transcript_features, evaluate_nmd_escape, annotate_cds

logger = logging.getLogger(__name__)

def annotate_nmd(vcf_path, gtf_path, fasta_path, output=None, reassign_exons=False, detailed=False):
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
    :param detailed: if True, save all CDS, transcript, and NMD prediction fields; if False, only save IDs and NMD result
    :param only_canonical: if True, filter to canonical transcripts only
    :return: DataFrame with NMD results
    """

    # read VCF file (variants)
    logger.info(f"Reading VCF file: {vcf_path}")
    vcf = read_vcf(vcf_path)
    logger.info(f"VCF shape: {vcf.df.shape}")

    # read GTF file (gene annotation)
    logger.info(f"Reading GTF file: {gtf_path}")
    gtf = read_gtf(gtf_path, reassign_exons=reassign_exons)
    logger.info(f"GTF File shape: {gtf.df.shape}")

    # read FASTA file (genome sequence)
    logger.info(f"Reading FASTA file: {fasta_path}")
    fasta = read_fasta(fasta_path)
    logger.info(f"FASTA file loaded with {len(fasta)} sequences.")

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
    logger.info("Calculating transcript features and evaluating NMD escape rules")
    results = []
    for cds in tqdm(cds_annotations, desc="Processing variants"):
        # Calculate transcript features
        features = calculate_transcript_features(cds)
        
        # Evaluate NMD escape
        prediction = evaluate_nmd_escape(cds, features)
        
        if detailed:
            # Combine all information
            result = {
                **asdict(cds),
                **asdict(features),
                **asdict(prediction)
            }
        else:
            # Minimal output: just IDs and NMD result
            result = {
                'variant_id': cds.variant_id,
                'transcript_id': cds.transcript_id,
                'gene_id': cds.gene_id,
                'chromosome': cds.chromosome,
                'start_variant': cds.start_variant,
                'end_variant': cds.end_variant,
                'ref': cds.ref,
                'alt': cds.alt,
                'alt_is_premature': cds.alt_is_premature,
                **asdict(features),
                **asdict(prediction),
            }
        results.append(result)
    
    nmd_results = pd.DataFrame(results)

    # Write output if output directory is specified
    if output is not None:
        vcf_base = os.path.splitext(os.path.basename(vcf_path))[0]
        output_file = os.path.join(output, f"{vcf_base}_nmdscanner.csv")
        logger.info(f"Writing results to {output_file}")
        nmd_results.to_csv(output_file, index=False)
    
    return nmd_results