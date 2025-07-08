import argparse
from nmd_scanner.scan import *
from nmd_scanner.rules import *
from nmd_scanner.analyze_gtf import *

def main(vcf_path, gtf_path, fasta_path, output):

    # read VCF file
    print(f"Reading VCF file: {vcf_path}")
    vcf = read_vcf(vcf_path)
    print(f"VCF shape: {vcf.df.shape}")

    # read GTF
    print(f"Reading GTF file: {gtf_path}")
    gtf = read_gtf(gtf_path)
    print(f"GTF File shape: {gtf.df.shape}")

    # read FASTA
    print(f"Reading FASTA file: {fasta_path}")
    fasta = read_fasta(fasta_path)
    print(f"FASTA File shape: {fasta.shape}")

    # analyze GTF and extract features
    cds = gtf[gtf.Feature == "CDS"]
    cds_df = cds.df

    #transcripts = gtf[gtf.Feature == "transcript"]
    #transcripts_df = transcripts.df

    exons = gtf[gtf.Feature == "exon"]
    exons_df = exons.df
    exons_df["exon_length"] = exons_df["End"] - exons_df["Start"]
    exon_counts = exons_df.groubpy("transcript_id").size().reset_index(name="num_exons_per_transcript")

    filtered_df_with_counts = exons_df.merge(exon_counts, on="transcript_id", how="left")
    filtered_df_with_counts["num_exons_per_transcript"] = filtered_df_with_counts["num_exons_per_transcript"].fillna(0).astype(int)


    # Apply the NMD rules
    print("Extracting PTCs and evaluating NMD escape rules...")
    results = extract_PTC(cds_df, vcf, fasta, exons_df)
    nmd_results = results.apply(evaluate_nmd_escape_rules, axis=1, result_type='expand')
    results = pd.concat([results, nmd_results], axis=1)

    # Determine output file path
    if os.path.isdir(output):
        os.makedirs(output, exist_ok=True)
        vcf_base = os.path.splitext(os.path.basename(vcf_path))[0]
        output_file = os.path.join(output, f"{vcf_base}_nmd_results.csv")
    else:
        # Ensure directory exists
        os.makedirs(os.path.dirname(output), exist_ok=True)
        output_file = output

    # Write output
    print(f"Writing results to {output_file}")
    results.to_csv(output_file, index=False)

    return results

    #filtered_df_with_counts = single_exon_rule(filtered_df_with_counts)
    #filtered_df_with_counts = long_exon_rule(filtered_df_with_counts, exon_length=400)
    #filtered_df_with_counts = last_exon_rule(filtered_df_with_counts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run NMD pipeline")
    parser.add_argument('--vcf', required=True, help='Path to VCF file')
    parser.add_argument('--gtf', required=True, help='Path to GTF file')
    parser.add_argument('--fasta', required=True, help='Path to FASTA file')
    parser.add_argument('--output', required=False, help='Path to output file or output directory')

    args = parser.parse_args()
    main(args.vcf, args.gtf, args.fasta, args.output)

