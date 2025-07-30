import argparse
from nmd_scanner.scan import *
from nmd_scanner.rules import *
from nmd_scanner.analyze_gtf import *
from nmd_scanner.extra_features import add_nmd_features
from pyfaidx import Fasta

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
    # fasta = read_fasta(fasta_path)
    fasta = Fasta(fasta_path)
    # print(f"FASTA File shape: {fasta.df.shape}")

    # analyze GTF and extract features
    cds = gtf[gtf.Feature == "CDS"]
    cds_df = cds.df

    #transcripts = gtf[gtf.Feature == "transcript"]
    #transcripts_df = transcripts.df

    exons = gtf[gtf.Feature == "exon"]
    exons_df = exons.df
    exons_df["exon_length"] = exons_df["End"] - exons_df["Start"]
    exon_counts = exons_df.groupby("transcript_id").size().reset_index(name="num_exons_per_transcript")

    filtered_df_with_counts = exons_df.merge(exon_counts, on="transcript_id", how="left")
    filtered_df_with_counts["num_exons_per_transcript"] = filtered_df_with_counts["num_exons_per_transcript"].fillna(0).astype(int)


    # Apply the NMD rules
    print("Extracting PTCs and evaluating NMD escape rules...")
    results = extract_ptc(cds_df, vcf, fasta, exons_df, output) # output of this (+ Zwischenschritte) saved in _final_ptc_analysis.tsv
    nmd_results = results.apply(evaluate_nmd_escape_rules, axis=1, result_type='expand')
    results = pd.concat([results, nmd_results], axis=1)

    # save as file --> can be removed later
    output_path = os.path.join(output, "nmd_rules.tsv")
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
    if os.path.exists(path):
        return True
    parent = os.path.dirname(path)
    return os.path.isdir(parent) if parent else False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run NMD pipeline")
    parser.add_argument('--vcf', required=True, help='Path to VCF file')
    parser.add_argument('--gtf', required=True, help='Path to GTF file')
    parser.add_argument('--fasta', required=True, help='Path to FASTA file')
    parser.add_argument('--output', required=True, help='Path to output file or output directory')

    args = parser.parse_args()

    # Validate output path before running main
    if not is_valid_output_path(args.output):
        raise ValueError(f"Invalid output path: {args.output}")

    main(args.vcf, args.gtf, args.fasta, args.output)

