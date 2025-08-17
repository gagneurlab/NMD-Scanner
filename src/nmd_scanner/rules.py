# Import dependencies
import pandas as pd
import pyranges as pr
import os
from Bio.Seq import Seq
from nmd_scanner import catch_sequence

# Main extract PTC script:
def extract_ptc(cds_df, vcf, fasta, exons_df, output):

    """
    Main function for extracting reference coding sequence, alternative coding sequence by incorporating the variant, analyzing for premature termination codons (PTCs),
    start and stop loss, and getting the transcript information.

    :param cds_df: CDS entries from the GTF file (Dataframe)
    :param vcf: Parsed VCF variant entries (PyRanges object)
    :param fasta: Reference genome sequence (pyfaidx.Fasta object)
    :param exons_df: All exonic entries from the GTF file (DataFrame)
    :param output: Output directory path for intermediate results (str)
    :return: analyze_transcript_df: Annotated dataframe with ref and alt CDS information, PTC analysis, start & stop loss analysis and transcript information
    """

    # Adjust the last 3 CDS positions to include stop codons
    cds_df_adj = adjust_last_cds_for_stop_codon(cds_df)
    print("Adjusting last CDS for stop codon: done.")

    # Intersect variants with CDS regions
    intersection_cds_vcf = pr.PyRanges(cds_df_adj).join(vcf, how=None, suffix="_variant").df
    print("Joining variants with cds entries: done.")

    ##########################################################################################
    # TODO: fix minus strand variants (only for TCGA and MMRF VCF!)
    # Fix REF and ALT for minus-strand CDSs
    # mask_minus_strand = intersection_cds_vcf["Strand"] == "-"
    # intersection_cds_vcf.loc[mask_minus_strand, "Ref"] = intersection_cds_vcf.loc[mask_minus_strand, "Ref"].apply(
    #    lambda seq: str(Seq(seq).reverse_complement()))
    # intersection_cds_vcf.loc[mask_minus_strand, "Alt"] = intersection_cds_vcf.loc[mask_minus_strand, "Alt"].apply(
    #    lambda seq: str(Seq(seq).reverse_complement()))
    ##########################################################################################

    # intersection = intersection_test.copy()
    # df3["Exon_CDS_seq"] = df3.apply(lambda row: fasta[row["Chromosome"]][row["Start"]:row["End"]].seq.upper(), axis=1)
    #intersection_cds_vcf["Exon_CDS_seq"] = [
    #    fasta[chrom][start:end].seq.upper()
    #    for chrom, start, end in zip(intersection_cds_vcf["Chromosome"], intersection_cds_vcf["Start"], intersection_cds_vcf["End"])
    #]

    # Fetch reference CDS sequence for each variant region
    print("Begin creating exon CDS sequence.")
    intersection_cds_vcf = catch_sequence.add_exon_cds_sequence(intersection_cds_vcf, fasta) # for faster access
    print("Creating exon CDS sequence: done.")

    # Apply variant to CDS and compute alternative CDS sequence and lengths
    intersection_cds_vcf[["Exon_CDS_length", "Exon_Alt_CDS_seq", "Exon_Alt_CDS_length"]] = intersection_cds_vcf.apply(
        apply_variant_edge_aware_with_lengths,
        axis=1
    )
    print("Creating exon CDS and alt CDS sequence: done.")

    ##### New ######
    # Filter out Variants with a reference mismatch
    mismatched_rows = intersection_cds_vcf[intersection_cds_vcf["Exon_Alt_CDS_seq"].isna()]
    print(f"\n[Warning] Skipping {len(mismatched_rows)} variants due to reference mismatches:")
    if not mismatched_rows.empty:
        print(mismatched_rows[["transcript_id", "Chromosome", "Start_variant", "End_variant", "Ref", "Alt"]].to_string(
            index=False))
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
    cds_df_adj = catch_sequence.add_exon_cds_sequence(cds_df_adj, fasta) # for faster access

    # save preliminary result
    #output_path = os.path.join(output, "2_cds_df_adj.tsv")
    #cds_df_adj.to_csv(output_path, sep="\t", index=False)
    #print(f"Creating exon CDS sequence for all exons for transcripts in df3: done. Saved in: {output_path}")


    # get full reference CDS per transcript by stiching exon CDS regions, plus alternative CDS with CDS exon information for both ref and alt
    results_df = create_reference_cds(intersection_cds_vcf, cds_df_adj)
    print("Create reference CDS: done.")

    # make intermediate output file of results_df
    #output_path = os.path.join(output, "3_create_reference_CDS.tsv")
    #results_df.to_csv(output_path, sep="\t", index=False)
    #print(f"Creating {output_path}: done.")

    # Get transcript sequence for relevant transcripts (speed up process) + length and transcript exon information (Tuple: exon number & exon length)
    exons_df = exons_df[exons_df["transcript_id"].isin(relevant_transcripts)].copy()
    exon_seqs = get_transcript_sequence(exons_df, fasta)
    print("Get transcript sequence: done.")

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

    results_df["cds_in_transcript"] = results_df.apply(check_cds_in_transcript, axis=1)

    # TODO: Analyze reference and alternative CDS for start / stop codons
    analysis_df = analyze_sequence(results_df)
    loss_df = start_stop_loss(analysis_df) # instead of loss_analysis_df (test for start stop loss)
    print("Analyzing sequence: done.")


    # Annotate transcript information (transcript start, end, sequence, length, exon info) in case of start or stop loss
    # transcript sequences are in: exon_seqs_subset
    print("Annotating transcript information in case of start/stop loss.")
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
    loss_df["alt_transcript_seq"] = loss_df.apply(
        lambda row: splice_alt_cds_into_transcript(row, row["transcript_seq"])
        if pd.notnull(row["transcript_seq"]) else None,
        axis=1
    )
    loss_df["alt_transcript_length"] = loss_df["alt_transcript_seq"].apply(
        lambda x: len(x) if pd.notnull(x) else None
    )


    # Add exon information to dataframe
    transcript_exon_info = exon_seqs.set_index("transcript_id")["transcript_exon_info"].to_dict()
    loss_df["transcript_exon_info"] = loss_df["transcript_id"].map(transcript_exon_info)


    # Analyze transcript sequence (e.g., frame, length, stop codon position, etc.) in case of start or stop loss
    analyze_transcript_df = analyze_transcript(loss_df)


    # save result
    #output_path = os.path.join(output, "5_final_ptc_analysis.tsv")
    #analyze_transcript_df.to_csv(output_path, sep="\t", index=False)
    #print(f"Save results in: {output_path}.")

    return analyze_transcript_df


# Functions used for extracting PTC:

def adjust_last_cds_for_stop_codon(df, exon_col="exon_number", transcript_col="transcript_id"):

    """
    Adjusts the genomic coordinates of the last CDS exon in each transcript by adding 3 positions, thus to include the stop codon.
    :param df: Dataframe containing CDS annotation
    :param exon_col: The name of the column that indicated the exon number, so we can find out which is the last CDS snippet
    :param transcript_col: The name of the column that indicates the transcript ID
    :return: Modified pandas DataFrame where the last exon of each transcript is extended by 3 bases to include the stop codon.
    """

    df = df.copy()

    # Search for the last exon in a transcript
    df[exon_col] = df[exon_col].astype(int)  # because otherwise e.g. 10 smaller than 9
    df.sort_values(by=[transcript_col, exon_col], inplace=True)
    last_exon_idx = df.groupby(transcript_col)[exon_col].idxmax()

    # Adjust start or end based on strand
    for idx in last_exon_idx:
        strand = df.at[idx, "Strand"]
        if strand == "+":
            df.at[idx, "End"] += 3
        elif strand == "-":
            df.at[idx, "Start"] -= 3

    return df

def apply_variant_edge_aware_with_lengths(row):

    """
    Applies a variant to a CDS exon sequence, taking into account not only SNVs but also partial overlaps at exon
    boundaries and computing the alternative sequence.

    :param row: A single row from the DataFrame (A pandas.Series) containing among others CDS Start and End, Variant
                Start and End, Ref, Alt, Exon_CDS_seq (Original CDS sequence as string)
    :return: The input pandas.Series with additional information:
             Exon_CDS_length (length of the original CDS),
             Exon_Alt_CDS_seq (alternative CDS after applying the variant / None if invalid),
             Exon_Alt_CDS_length (length of the alternative CDS / None if invalid).
    """

    cds_seq = list(row["Exon_CDS_seq"])
    strand = row["Strand"]
    ref = row["Ref"]
    alt = row["Alt"]
    cds_start = int(row["Start"])
    cds_end = int(row["End"])
    var_start = int(row["Start_variant"])
    var_end = int(row["End_variant"])

    # Determine the overlap between variant and this CDS region
    overlap_start = max(var_start, cds_start)
    overlap_end = min(var_end, cds_end)

    # If there is no overlap between the variant and this CDS region
    if overlap_start >= overlap_end:
        return pd.Series({
            "Exon_CDS_length": len(cds_seq),
            "Exon_Alt_CDS_seq": None,
            "Exon_Alt_CDS_length": None
        })

    # Position of overlap within the CDS
    cds_index = overlap_start - cds_start
    overlap_len = overlap_end - overlap_start

    # Offset of the overlapping region within the variant
    ref_offset = overlap_start - var_start
    ref_in_cds = ref[ref_offset:ref_offset + overlap_len]
    alt_in_cds = alt[ref_offset:ref_offset + overlap_len]

    # Determine if there’s leftover alt outside CDS (insertions at end)
    extra_alt = ""
    if len(alt) > len(ref):
        # Limit extra_alt to what corresponds to CDS overlap
        extra_start = ref_offset + overlap_len
        if var_end > cds_end:
            # Only include alt bases that map to CDS
            remaining_cds_len = cds_end - overlap_end
            extra_alt = alt[extra_start:extra_start + remaining_cds_len]
        else:
            extra_alt = alt[extra_start:]

    # Confirm that the reference matches
    cds_ref_part = "".join(cds_seq[cds_index:cds_index + overlap_len])
    if cds_ref_part != ref_in_cds.upper(): # reference mismatch
        return pd.Series({
            "Exon_CDS_length": len(cds_seq),
            "Exon_Alt_CDS_seq": None,
            "Exon_Alt_CDS_length": None
        })

    # Build alternative sequence
    alt_seq = []
    alt_seq.extend(cds_seq[:cds_index]) # Copy the CDS up to the variant position

    if len(ref) == len(alt): # Substitution
        alt_seq.extend(list(alt_in_cds))
    elif len(alt) > len(ref): # Insertion
        alt_seq.extend(list(alt_in_cds))
        alt_seq.extend(list(extra_alt))
    elif len(ref) > len(alt): # Deletion
        alt_seq.extend(list(alt_in_cds))

    # Add remaining CDS sequence after variant
    alt_seq.extend(cds_seq[cds_index + overlap_len:])

    return pd.Series({
        "Exon_CDS_length": len(cds_seq),
        "Exon_Alt_CDS_seq": "".join(alt_seq),
        "Exon_Alt_CDS_length": len(alt_seq)
    })

def create_reference_cds(intersection_cds_vcf, cds_df_test):

    """
    Constructs the whole CDS sequence (multiple exons) for transcripts affected by a variant, both in their reference
    and alternative form.
    :param intersection_cds_vcf: DataFrame containing variant-CDS intersection and corresponding alternative CDS sequences
                                 includes: transcript_id, Exon_CDS_seq + length, Exon_ALT_CDS_seq + length
    :param cds_df_test: Reference exon-level CDS data for all transcripts with exon_number
                        includes: transcript_id, exon_number, Start, End, Strand, Exon_CDS_seq
    :return: DataFrame with one row per variant-transcript pair, containing full reference and alternative CDS + lengths,
             exon-wise CDS information as tuple (exon number, exon-wise CDS length)
    """

    results = []

    for transcript_id, var_df in intersection_cds_vcf.groupby("transcript_id"):  # Only transcripts with a variant

        # 1. Get reference exons
        ref_exons = cds_df_test[cds_df_test["transcript_id"] == transcript_id].copy()
        ref_exons = ref_exons.sort_values("Start")

        # Get reference CDS sequence start und stop position for finding position in transcript sequence
        ref_cds_start = ref_exons["Start"].min()
        ref_cds_stop = ref_exons["End"].max()

        # Join reference exon sequences to form full CDS sequence
        ref_seq = "".join(ref_exons["Exon_CDS_seq"].tolist())

        ######
        # Since I sometimes get errors in the following code snippet because of NaN values,
        # we print them but leave them in our dataframe for now
        #nan_rows = ref_exons[ref_exons["Exon_CDS_seq"].isna()]
        #if not nan_rows.empty:
        #    print(f"\n[Warning] Found {len(nan_rows)} NaN Exon_CDS_seq entries in transcript: {transcript_id}")
        #    print(nan_rows.to_string(index=False))

        # Collect exon numbers and lengths (for tracking exon contribution later on)
        # ref_cds_lengths = [len(seq) for seq in ref_exons["Exon_CDS_seq"].tolist()]
        ref_cds_info = [
            (row["exon_number"], len(row["Exon_CDS_seq"]))
            for _, row in ref_exons.iterrows()
        #    if pd.notna(row["Exon_CDS_seq"])  # sometimes we get NaN errors
        ]
        ######

        # Get strand info (all should be the same within transcript)
        strand = ref_exons["Strand"].iloc[0]

        for variant, cds_df in var_df.groupby(["Chromosome", "Start_variant", "End_variant", "Ref", "Alt"],
                                              observed=True):

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

            ######
            # Since I sometimes get errors in the following code snippet because of NaN values,
            # we print them but leave them in our dataframe for now (same as before)
            #nan_alt_rows = alt_exons[alt_exons["Exon_CDS_seq"].isna()]
            #if not nan_alt_rows.empty:
            #    print(f"\n[Warning] NaN Exon_CDS_seq values found in alt_exons for variant in transcript: {transcript_id}")
            #    print(nan_alt_rows.to_string(index=False))

            # alt_cds_lengths = [len(seq) for seq in alt_exons["Exon_CDS_seq"].tolist()]
            alt_cds_info = [
                (row["exon_number"], len(row["Exon_CDS_seq"]))
                for _, row in alt_exons.iterrows()
            #    if pd.notna(row["Exon_CDS_seq"])  # sometimes we get NaN errors
            ]
            ######

            # Get alternative CDS sequence start and stop position for finding position in transcript sequence
            alt_cds_start = alt_exons["Start"].min()
            alt_cds_stop = alt_exons["End"].max()

            alt_seq = "".join(alt_exons["Exon_CDS_seq"].tolist())

            # Apply reverse complement if on minus strand
            if strand == "-":
                ref_seq_final = str(Seq(ref_seq).reverse_complement())
                alt_seq_final = str(Seq(alt_seq).reverse_complement())
            else:
                ref_seq_final = ref_seq
                alt_seq_final = alt_seq

            # Append to results
            results.append({
                "transcript_id": transcript_id,
                "variant_id": var_row["ID"],

                "ref_cds_start": ref_cds_start,
                "ref_cds_stop": ref_cds_stop,
                "ref_cds_seq": ref_seq_final,
                "ref_cds_len": len(ref_seq_final),

                "alt_cds_start": alt_cds_start,
                "alt_cds_stop": alt_cds_stop,
                "alt_cds_seq": alt_seq_final,
                "alt_cds_len": len(alt_seq_final),

                "chromosome": var_row["Chromosome"],
                "gene_id": var_row["gene_id"],
                "strand": strand,
                "ref": var_row["Ref"],
                "alt": var_row["Alt"],
                "start_variant": var_row["Start_variant"],
                "end_variant": var_row["End_variant"],

                "ref_cds_info": ref_cds_info,
                "alt_cds_info": alt_cds_info,

                # "ref_cds_lengths": [length for exon_num, length in ref_cds_info],
                # "alt_cds_lengths": [length for exon_num, length in alt_cds_info]
            })

    results_df = pd.DataFrame(results)
    return(results_df)

def get_transcript_sequence(exons_df, fasta):

    """
    Construct full transcript sequences by concatenating the exon sequences from the FASTA genome reference, grouped by transcript.
    Get transcript length and transcript information as well.
    :param exons_df: DataFrame containing exon-level annotations from the GTF file.
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
            print(f"Unknown strand for {transcript_id}")
            continue

        # Sort by exon start coordinate (strand not considered here yet)
        group_sorted = group.sort_values(by="Start").copy()

        seq_parts = [] # to accumulate exon sequences
        starts = [] # for overall transcript start
        ends = [] # for overall transcript end
        exon_info = [] # for tracking exon_number and length

        # fetch exon sequence and metadata
        for _, row in group_sorted.iterrows():
            chrom = row["Chromosome"]
            start = int(row["Start"])
            end = int(row["End"])
            exon_number = row["exon_number"]

            starts.append(start)
            ends.append(end)

            # Fetch exon sequence from fasta reference genome
            exon_seq = fasta[chrom][start:end] #.seq
            exon_seq_str = str(exon_seq).upper()
            seq_parts.append(exon_seq_str)

            exon_info.append((exon_number, len(exon_seq_str)))

        # join exon sequences into a full transcript sequence
        joined_seq = "".join(seq_parts)

        # Apply reverse complement for minus strand transcripts
        if strand == "-":
            joined_seq = str(Seq(joined_seq).reverse_complement())
            exon_info = exon_info[::-1]

        exon_data.append({
            "Chromosome": chrom,
            "transcript_id": transcript_id,
            "start": min(starts),
            "end": max(ends),
            "strand": strand,
            "transcript_sequence": joined_seq,
            "transcript_length": len(joined_seq),
            "transcript_exon_info": exon_info
        })

    exon_seqs = pd.DataFrame(exon_data)
    return(exon_seqs)

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

    :param results_df: DataFrame containing CDS sequences and exon information for both reference and alternative sequences, per variant
    :return: DataFrame with added annotation columns for reference and alternative sequence separately:
             such as start codon position / exon, last codon and its validity as stop codon, first in-frame stop codon + position,
             number and information of all available stop codons, premature stop codon flag
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
                codon = seq[i:i + 3]
                if codon == start_codon and start_pos is None: # first start codon position
                    start_pos = i
                if codon in valid_stop_codons: # record all stop codons with their positions and exons
                    stop_codons.append((i, codon))
                    stop_exons.append(get_exon(i, exon_info))  # for exon number

            last_codon = seq[-3:]
            is_valid_stop = last_codon in valid_stop_codons
            first_stop_pos = stop_codons[0][0] if stop_codons else None
            first_stop = stop_codons[0][1] if stop_codons else None
            is_premature = first_stop_pos is not None and first_stop_pos < len(seq) - 3
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

    # Start codon loss: reference sequence has a start codon, alternative sequence does not or the position is changed
    df["start_loss"] = (
                               (df["ref_start_codon_pos"].notna()) & df["alt_start_codon_pos"].isna()
                       ) | (
                               df["ref_start_codon_pos"] != df["alt_start_codon_pos"]
                       )

    # Stop codon loss: reference sequence had a valid stop codon, the alternative sequence does not or the position is changed
    df["stop_loss"] = (
                              (df["ref_valid_stop"] == True) & (df["alt_valid_stop"] != True)
                      ) | (
                              df["ref_last_codon"] != df["alt_last_codon"]  # Or take this out?
                      )

    return df

def splice_alt_cds_into_transcript(row, transcript_seq):

    """
    Splice the alternative CDS sequence into the full transcript sequence to create the alternative transcript
    :param row: A pd.Series row containing "ref_cds_seq" (Reference CDS) and "alt_cds_seq" (Alternative / Variant-modified CDS)
    :param transcript_seq: Full transcript sequence
    :return: Modified (alternative) transcript sequence with the alternative CDS spliced in the correct position,
             or None if no match is found
    """

    # Step 1: search for ref_cds_seq match in the transcript, and replace that with the alt_cds_seq
    ref_cds_seq = row["ref_cds_seq"].upper()
    alt_cds_seq = row["alt_cds_seq"].upper()

    # Find the ref CDS in the transcript sequence
    ref_start_idx = transcript_seq.find(ref_cds_seq)

    if ref_start_idx == -1:
        return None  # Cannot find ref CDS, alignment problem

    ref_end_idx = ref_start_idx + len(ref_cds_seq)

    # Replace the reference CDS with the variant-modified / alternative one
    new_transcript_seq = (
        transcript_seq[:ref_start_idx] +
        alt_cds_seq +
        transcript_seq[ref_end_idx:]
    )

    return new_transcript_seq

def analyze_transcript(results_df):

    """
    Analyze the alternative transcript sequence in cases of start or stop codons loss due to mutations.
    Scan for new in-frame start or stop codons in the alternative transcript sequence.
    :param results_df: DataFrame containing transcript sequence data and annotations, including start_loss and stop_loss flags
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
        cds_start = row["alt_cds_start"] - row["transcript_start"]  # start analysis at cds start position
        # cds_stop = cds_start + len(seq)

        exon_info = row["transcript_exon_info"]  # for exon number

        # Skip rows with invalid or too-short sequences
        if not isinstance(seq, str) or len(seq) < 3:
            continue

        start_pos = None
        start_exon = None  # for exon number
        stop_codons = []
        stop_exons = []  # for exon number

        # only analyze rows flagged with start or stop codon loss: skip the others and fill with None values
        # Skips only if we have both start_loss = FALSE and stop_loss = FALSE. If one is true, then don't skip.
        if not (row["start_loss"] or row["stop_loss"]):
            continue

        # START LOSS rescue search
        if row["start_loss"]:
            # Walk through sequence starting at CDS start with +1 positions until start codon is found
            for i in range(cds_start, len(seq) - 2):
                codon = seq[i:i + 3]
                if codon == start_codon:
                    start_pos = i

                    start_exon = get_exon(start_pos, exon_info)  # for exon number

                    # From new start codon, scan codons in frame
                    for j in range(i, len(seq) - 2, 3):
                        codon2 = seq[j:j + 3]
                        if codon2 in valid_stop_codons:
                            stop_codons.append((j, codon2))
                            stop_exons.append(get_exon(j, exon_info))  # for exon number

                    break

        # STOP LOSS rescue search
        elif row["stop_loss"]:
            # Start at cds_stop, scan codons in frame to end --> lets start at cds start so we are in frame
            for i in range(cds_start, len(seq) - 2, 3):
                codon = seq[i:i + 3]
                if codon == start_codon and start_pos is None:
                    start_pos = i
                    start_exon = get_exon(start_pos, exon_info)  # for exon number
                if codon in valid_stop_codons:
                    stop_codons.append((i, codon))
                    stop_exons.append(get_exon(i, exon_info))  # for exon number

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

def evaluate_nmd_escape_rules(row):

    """
    Evaluate whether a premature stop codon in a transcript is likely to escape nonsense-mediated decay (NMD) based on
    established biological rules. This function applies five NMD escape rules to determine if a premature termination
    codon (PTC) is likely to escape degradation:
    1. Last exon rule: The PTC is in the last exon
    2. 50nt penultimate rule: The PTC is within 50 nucleotides upstream of the last exon junction
    3. Long exon rule: The PTC is in an exon with >407 nucleotides
    4. Start proximal rule: The PTC is within 150 nucleotides of the start codon
    5. Single exon rule: The transcript where the PTC lays consists only of a single exon
    A PTC is considered to escape NMD if it satisfies any of the above rules.

    :param row: A row of the DataFrame including alt_is_premature (bool), alt_first_stop_pos (int),
                alt_stop_codon_exons (list[int]), transcript_exon_info (list[tuple[exon_number (int), exon_length (int)]]),
                alt_start_codon_pos (int)
    :return: A dictionary with boolean flags for each rule and overall NMD escape
    """

    # Only relevant for premature stop codons
    if not row.get("alt_is_premature"):
        return {
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": False
        }

    # Extract relevant data
    stop_pos = row.get("alt_first_stop_pos")
    stop_exons = row.get("alt_stop_codon_exons") or []
    exon_info = row.get("transcript_exon_info") or []
    start_pos = row.get("alt_start_codon_pos")

    # Preprocess exon info
    sorted_exons = sorted(exon_info, key=lambda x: x[0])
    exon_length_map = {exon_num: length for exon_num, length in sorted_exons}
    exon_offsets = {}
    offset = 0
    for exon_num, length in sorted_exons:
        exon_offsets[exon_num] = (offset, offset + length)
        offset += length

    # Single exon rule
    rule_single_exon = len(sorted_exons) == 1

    # Last exon rule
    if not rule_single_exon:
        last_exon = sorted_exons[-1][0] if sorted_exons else None
        rule_last_exon = bool(stop_exons and max(stop_exons) == last_exon)
    else:
        rule_last_exon = False

    # 50nt from penultimate exon end
    if len(sorted_exons) >= 2:
        penultimate_exon_num, penultimate_len = sorted_exons[-2]
        pen_start, pen_end = exon_offsets.get(penultimate_exon_num, (None, None))
        rule_50nt_penultimate = pen_end is not None and (stop_pos >= pen_end - 50) and (stop_pos < pen_end)
    else:
        rule_50nt_penultimate = False

    # Long exon rule (with exon longer than >407nt)
    rule_long_exon = any(exon_length_map.get(exon, 0) > 407 for exon in stop_exons)

    # Start-proximal rule (closer than 150nt from the start codon)
    rule_start_proximal = start_pos is not None and stop_pos is not None and (stop_pos - start_pos) < 150 and (stop_pos - start_pos) >= 0

    # NMD escape if any rule is true
    escape = rule_last_exon or rule_50nt_penultimate or rule_long_exon or rule_start_proximal or rule_single_exon

    return {
        "nmd_last_exon_rule": rule_last_exon,
        "nmd_50nt_penultimate_rule": rule_50nt_penultimate,
        "nmd_long_exon_rule": rule_long_exon,
        "nmd_start_proximal_rule": rule_start_proximal,
        "nmd_single_exon_rule": rule_single_exon,
        "nmd_escape": escape
    }
