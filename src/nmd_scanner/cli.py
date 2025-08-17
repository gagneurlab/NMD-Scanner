# Import dependencies
import argparse
import os
import pandas as pd
from nmd_scanner.scan import read_vcf, read_gtf
from nmd_scanner.rules import extract_ptc, evaluate_nmd_escape_rules
from nmd_scanner.extra_features import add_nmd_features
from pyfaidx import Fasta

def main(vcf_path, gtf_path, fasta_path, output):

    """
    Main function for NMD scanner

    Steps:
    1. Read input files (VCF, GTF, FASTA)
    2. Parse and preprocess gene annotations (CDS, exons)
    3. Get number of exons per transcript
    4. Extract premature termination codons (PTCs) & Evaluate NMD escape rules
    5. Add extra features to output (e.g. 3' & 5'UTR length, downstream & upstream exon counts, etc.)
    6. Return and Save output results

    :param vcf_path: path to the input VCF file
    :param gtf_path: path to the input GTF annotation file
    :param fasta_path: path to the reference FASTA file
    :param output: directory to save the results in
    :return: CSV file summarizing all annotated variants
    """

    # read VCF file (variants)
    print(f"Reading VCF file: {vcf_path}")
    vcf = read_vcf(vcf_path)
    print(f"VCF shape: {vcf.df.shape}")

    # read GTF file (gene annotation)
    print(f"Reading GTF file: {gtf_path}")
    gtf = read_gtf(gtf_path)
    print(f"GTF File shape: {gtf.df.shape}")

    # read FASTA file (genome sequence)
    print(f"Reading FASTA file: {fasta_path}")
    fasta = Fasta(fasta_path)

    # extract CDS regions from the GTF file
    cds = gtf[gtf.Feature == "CDS"]
    cds_df = cds.df

    #transcripts = gtf[gtf.Feature == "transcript"]
    #transcripts_df = transcripts.df

    # extract exon regions from the GTF file and compute exon related metrics:
    # exon length & number of exons contained in each transcript
    exons = gtf[gtf.Feature == "exon"]
    exons_df = exons.df
    exons_df["exon_length"] = exons_df["End"] - exons_df["Start"]
    # I don't think I am using "num_exon_per_transcript" --> leave out or build in --> results in same df
    #exon_counts = exons_df.groupby("transcript_id").size().reset_index(name="num_exons_per_transcript")
    #filtered_df_with_counts = exons_df.merge(exon_counts, on="transcript_id", how="left")
    #filtered_df_with_counts["num_exons_per_transcript"] = filtered_df_with_counts["num_exons_per_transcript"].fillna(0).astype(int)

    # Apply the NMD rules
    print("Extracting PTCs and evaluating NMD escape rules...")
    results = extract_ptc(cds_df, vcf, fasta, exons_df, output) # output of this (+ Zwischenschritte) saved in _final_ptc_analysis.tsv
    nmd_results = results.apply(evaluate_nmd_escape_rules, axis=1, result_type='expand')
    results = pd.concat([results, nmd_results], axis=1)

    # Save intermediate NMD rule output --> can be removed later
    #output_path = os.path.join(output, "6_nmd_rules.tsv")
    #results.to_csv(output_path, sep="\t", index=False)
    #print(f"Save nmd rules results in: {output_path}.")

    # Add additional features (inspired by benchmark dataset)
    extra_features = results.apply(add_nmd_features, axis=1, result_type='expand')
    results = pd.concat([results, extra_features], axis=1)


    # # adjust datatypes --> TODO: adjust in code, not afterwards
    # dtype_mapping = {
    #     "transcript_id": "object",
    #     "variant_id": "object",
    #     "ref_cds_start": "Int64",
    #     "ref_cds_stop": "Int64",
    #     "ref_cds_seq": "string",
    #     "ref_cds_len": "Int64",
    #     "alt_cds_start": "Int64",
    #     "alt_cds_stop": "Int64",
    #     "alt_cds_seq": "string",
    #     "alt_cds_len": "Int64",
    #     "chromosome": "object",
    #     "gene_id": "object",
    #     "strand": "object",
    #     "ref": "object",
    #     "alt": "object",
    #     "start_variant": "Int64",
    #     "end_variant": "Int64",
    #     "ref_cds_info": "object",
    #     "alt_cds_info": "object",
    #     "cds_in_transcript": "boolean",
    #     "ref_start_codon_pos": "Int64",
    #     "ref_start_codon_exon": "Int64",
    #     "ref_last_codon": "object",
    #     "ref_valid_stop": "boolean",
    #     "ref_first_stop_codon": "object",
    #     "ref_first_stop_pos": "Int64",
    #     "ref_num_stop_codons": "Int64",
    #     "ref_all_stop_codons": "object",
    #     "ref_stop_codon_exons": "object",
    #     "ref_is_premature": "boolean",
    #     "alt_start_codon_pos": "Int64",
    #     "alt_start_codon_exon": "Int64",
    #     "alt_last_codon": "object",
    #     "alt_valid_stop": "boolean",
    #     "alt_first_stop_codon": "object",
    #     "alt_first_stop_pos": "Int64",
    #     "alt_num_stop_codons": "Int64",
    #     "alt_all_stop_codons": "object",
    #     "alt_stop_codon_exons": "object",
    #     "alt_is_premature": "boolean",
    #     "start_loss": "boolean",
    #     "stop_loss": "boolean",
    #     "transcript_start": "Int64",
    #     "transcript_end": "Int64",
    #     "transcript_seq": "string",
    #     "transcript_length": "Int64",
    #     "alt_transcript_seq": "string",
    #     "alt_transcript_length": "Int64",
    #     "transcript_exon_info": "object",
    #     "transcript_start_codon_pos": "Int64",
    #     "transcript_start_codon_exon": "Int64",
    #     "transcript_last_codon": "object",
    #     "transcript_valid_stop": "boolean",
    #     "transcript_first_stop_codon": "object",
    #     "transcript_first_stop_pos": "Int64",
    #     "transcript_num_stop_codons": "Int64",
    #     "transcript_all_stop_codons": "object",
    #     "transcript_stop_codon_exons": "object",
    #     "nmd_last_exon_rule": "boolean",
    #     "nmd_50nt_penultimate_rule": "boolean",
    #     "nmd_long_exon_rule": "boolean",
    #     "nmd_start_proximal_rule": "boolean",
    #     "nmd_single_exon_rule": "boolean",
    #     "nmd_escape": "boolean",
    #     "utr3_length": "Int64",
    #     "utr5_length": "Int64",
    #     "total_exon_count": "Int64",
    #     "upstream_exon_count": "Int64",
    #     "downstream_exon_count": "Int64",
    #     "ptc_to_start_codon": "Int64",
    #     "ptc_less_than_150nt_to_start": "boolean",
    #     "ptc_exon_length": "Int64",
    #     "stop_codon_distance_nmd": "Int64"
    # }
    # for col, dtype in dtype_mapping.items():
    #     if col in results.columns:
    #         results[col] = results[col].astype(dtype)


    # Write output
    vcf_base = os.path.splitext(os.path.basename(vcf_path))[0]
    output_file = os.path.join(output, f"{vcf_base}_final_nmd_results.csv")
    print(f"Writing results to {output_file}")
    results.to_csv(output_file, index=False)

    return results

def is_valid_output_path(path):

    """
    Validate if the output path exists or is creatable.
    Allows a file path (where the parent directory must exist) or a directory.
    """

    if os.path.exists(path):
        return True
    parent = os.path.dirname(path)
    return os.path.isdir(parent) if parent else False

if __name__ == '__main__':

    # CLI argument parser
    parser = argparse.ArgumentParser(description="Run NMD pipeline")
    parser.add_argument('--vcf', required=True, help='Path to VCF file')
    parser.add_argument('--gtf', required=True, help='Path to GTF file')
    parser.add_argument('--fasta', required=True, help='Path to FASTA file')
    parser.add_argument('--output', required=True, help='Path to output file or output directory')

    args = parser.parse_args()

    # Check that the output path is valid
    if not is_valid_output_path(args.output):
        raise ValueError(f"Invalid output path: {args.output}")

    # Run the main pipeline
    main(args.vcf, args.gtf, args.fasta, args.output)

