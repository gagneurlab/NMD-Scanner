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

    # read GTF fil (gene annotation)
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
    output_path = os.path.join(output, "6_nmd_rules.tsv")
    results.to_csv(output_path, sep="\t", index=False)
    print(f"Save nmd rules results in: {output_path}.")

    # Add additional features (inspired by benchmark dataset)
    extra_features = results.apply(add_nmd_features, axis=1, result_type='expand')
    results = pd.concat([results, extra_features], axis=1)


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

