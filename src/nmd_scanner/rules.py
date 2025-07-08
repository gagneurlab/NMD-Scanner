import pandas as pd
import pyranges as pr

# def single_exon_rule(filtered_df_with_counts):
#     filtered_df_with_counts["NMD_is_single_exon"] = filtered_df_with_counts["num_exons_per_transcript"] == 1
#     return filtered_df_with_counts
#
# def long_exon_rule(filtered_df_with_counts, exon_length = 400):
#     filtered_df_with_counts["NMD_is_long_exon"] = filtered_df_with_counts["exon_length"] >= exon_length
#     return filtered_df_with_counts
#
# def last_exon_rule(filtered_df_with_counts):
#     # Ensure exon_number is numeric
#     filtered_df_with_counts["exon_number"] = pd.to_numeric(filtered_df_with_counts["exon_number"], errors="coerce")
#
#     # Compute the max exon number per transcript
#     max_exon = filtered_df_with_counts.groupby("transcript_id")["exon_number"].transform("max")
#
#     # Mark whether the exon is the last (=max) one
#     filtered_df_with_counts["NMD_is_last_exon"] = filtered_df_with_counts["exon_number"] == max_exon
#
#     return filtered_df_with_counts

# Main extract PTC script:
def extract_PTC(cds_df, hg38_example, fasta, exons_df):
    #intersection = transcripts.join(hg38_example, how=None, suffix="_variant")
    #df = intersection.df

    # Drop duplicates to get one strand per transcript_id
    #transcript_strands = df[["transcript_id", "Strand"]].drop_duplicates()

    # Optionally, convert to a dictionary: {transcript_id: strand}
    #transcript_strand_dict = dict(zip(transcript_strands["transcript_id"], transcript_strands["Strand"]))

    cds_df_test = adjust_last_cds_for_stop_codon(cds_df)

    intersection_test = pr.PyRanges(cds_df_test).join(hg38_example, how=None,
                                                      suffix="_variant").df  # test with two transcript example

    df3 = intersection_test.copy()
    df3["Exon_CDS_seq"] = df3.apply(lambda row: fasta[row["Chromosome"]][row["Start"]:row["End"]].seq.upper(), axis=1)

    df3[["Exon_CDS_length", "Exon_Alt_CDS_seq", "Exon_Alt_CDS_length"]] = df3.apply(
        apply_variant_edge_aware_with_lengths,
        axis=1
    )

    # get reference sequence
    results_df = create_reference_CDS(df3, cds_df_test)

    # get transcript sequence
    exon_seqs = get_transcript_sequence(exons_df, fasta)

    # check if the CDS sequences we generated before are contained in the transcript sequence,
    # so I know whether the transcript sequence was computed correctly
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

    analysis_df = analyze_sequence(results_df)
    loss_df = start_stop_loss(analysis_df) # instead of loss_analysis_df (test for start stop loss)


    # 2. Step: loop through the dataframe and check for start_loss or stop_loss and handle it
    # transcript sequences are in: exon_seqs_subset

    # add information such as transcript start and end
    transcript_starts = exon_seqs.set_index("transcript_id")["start"].to_dict()
    loss_df["transcript_start"] = loss_df["transcript_id"].map(transcript_starts)

    transcript_ends = exon_seqs.set_index("transcript_id")["end"].to_dict()
    loss_df["transcript_end"] = loss_df["transcript_id"].map(transcript_ends)

    # add trancript sequence to loss_df
    transcript_sequences = exon_seqs.set_index("transcript_id")[
        "transcript_sequence"].to_dict()  # create map of transcript-id to transcript sequence
    loss_df["transcript_seq"] = loss_df["transcript_id"].map(transcript_sequences)

    # Add transcript length to dataframe
    transcript_lengths = exon_seqs.set_index("transcript_id")["transcript_length"].to_dict()
    loss_df["transcript_length"] = loss_df["transcript_id"].map(transcript_lengths)

    # apply splice_alt_cds_into_transcript
    loss_df["alt_transcript_seq"] = loss_df.apply(
        lambda row: splice_alt_cds_into_transcript(row, row["transcript_seq"])
        if pd.notnull(row["transcript_seq"]) else None,
        axis=1
    )

    loss_df["alt_transcript_length"] = loss_df["alt_transcript_seq"].apply(
        lambda x: len(x) if pd.notnull(x) else None
    )

    # Add exon info to dataframe
    transcript_exon_info = exon_seqs.set_index("transcript_id")["transcript_exon_info"].to_dict()
    loss_df["transcript_exon_info"] = loss_df["transcript_id"].map(transcript_exon_info)

    # Analyze transcript sequence in case of start or stop loss
    analyze_transcript_df = analyze_transcript(loss_df)


    return analyze_transcript_df # (?)


# Functions used for extracting PTC:
def adjust_last_cds_for_stop_codon(df, exon_col="exon_number", transcript_col="transcript_id"):
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

    # To create the CDS sequence

    cds_seq = list(row["Exon_CDS_seq"])
    strand = row["Strand"]
    ref = row["Ref"]
    alt = row["Alt"]
    cds_start = int(row["Start"])
    cds_end = int(row["End"])
    var_start = int(row["Start_variant"])
    var_end = int(row["End_variant"])

    # Overlap of the variant with the CDS
    overlap_start = max(var_start, cds_start)
    overlap_end = min(var_end, cds_end)

    if overlap_start >= overlap_end:
        return pd.Series({
            "Exon_CDS_length": len(cds_seq),
            "Exon_Alt_CDS_seq": None,
            "Exon_Alt_CDS_length": None
        })

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
    if cds_ref_part != ref_in_cds.upper():
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

    # Add the rest of the CDS
    alt_seq.extend(cds_seq[cds_index + overlap_len:])

    return pd.Series({
        "Exon_CDS_length": len(cds_seq),
        "Exon_Alt_CDS_seq": "".join(alt_seq),
        "Exon_Alt_CDS_length": len(alt_seq)
    })


def create_reference_CDS(df3, cds_df_test):
    results = []

    for transcript_id, var_df in df3.groupby("transcript_id"):  # Only transcripts with a variant

        # 1. Get reference exons
        ref_exons = cds_df_test[cds_df_test["transcript_id"] == transcript_id].copy()
        ref_exons = ref_exons.sort_values("Start")

        # Get reference CDS sequence start und stop position for finding position in transcript sequence
        ref_cds_start = ref_exons["Start"].min()
        ref_cds_stop = ref_exons["End"].max()

        # Join reference sequence
        ref_seq = "".join(ref_exons["Exon_CDS_seq"].tolist())

        ######
        # TODO: save exon lengths for ref seq --> what we do here is to save the cds chunk lengths I think
        # ref_cds_lengths = [len(seq) for seq in ref_exons["Exon_CDS_seq"].tolist()]
        ref_cds_info = [
            (row["exon_number"], len(row["Exon_CDS_seq"]))
            for _, row in ref_exons.iterrows()
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
            # TODO: save exon lengths for alt seq --> what we do here is to save the cds chunk lengths I think
            # alt_cds_lengths = [len(seq) for seq in alt_exons["Exon_CDS_seq"].tolist()]
            alt_cds_info = [
                (row["exon_number"], len(row["Exon_CDS_seq"]))
                for _, row in alt_exons.iterrows()
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
    exon_data = []

    for transcript_id, group in exons_df.groupby("transcript_id"):
        strand = group.iloc[0]["Strand"]

        if strand not in ["+", "-"]:
            print(f"Unknown strand for {transcript_id}")
            continue

        # Sort by exon start coordinate — not affected by strand yet
        group_sorted = group.sort_values(by="Start").copy()

        seq_parts = []
        starts = []
        ends = []
        exon_info = []

        for _, row in group_sorted.iterrows():
            chrom = row["Chromosome"]
            start = int(row["Start"])
            end = int(row["End"])
            exon_number = row["exon_number"]

            starts.append(start)
            ends.append(end)

            # Fetch exon sequence from fasta reference genome
            exon_seq = fasta[chrom][start:end].seq
            exon_seq_str = str(exon_seq).upper()
            seq_parts.append(exon_seq_str)

            exon_info.append((exon_number, len(exon_seq_str)))

        joined_seq = "".join(seq_parts)

        # reverse complement for minus strand
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

    # Create DataFrame
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
    valid_stop_codons = {"TAA", "TAG", "TGA"}
    start_codon = "ATG"

    df = results_df.copy()

    # new columns: extend old dataframe by these
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

            if not isinstance(seq, str) or len(seq) < 3:
                continue

            start_pos = None
            stop_codons = []
            stop_exons = []  # for exon number

            # Scan in codons (step=3)
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i:i + 3]
                if codon == start_codon and start_pos is None:
                    start_pos = i
                if codon in valid_stop_codons:
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
    df = df.copy()

    df["start_loss"] = (
                               (df["ref_start_codon_pos"].notna()) & df["alt_start_codon_pos"].isna()
                       ) | (
                               df["ref_start_codon_pos"] != df["alt_start_codon_pos"]
                       )

    df["stop_loss"] = (
                              (df["ref_valid_stop"] == True) & (df["alt_valid_stop"] != True)
                      ) | (
                              df["ref_last_codon"] != df["alt_last_codon"]  # Or take this out?
                      )

    return df

def splice_alt_cds_into_transcript(row, transcript_seq):

    # Step 1: search for ref_cds_seq match in the transcript, and replace that with the alt_cds_seq

    ref_cds_seq = row["ref_cds_seq"].upper()
    alt_cds_seq = row["alt_cds_seq"].upper()

    # Find the ref CDS in the transcript sequence
    ref_start_idx = transcript_seq.find(ref_cds_seq)

    if ref_start_idx == -1:
        return None  # Cannot find ref CDS, alignment problem

    ref_end_idx = ref_start_idx + len(ref_cds_seq)

    # Replace the CDS with the variant-modified one
    new_transcript_seq = (
        transcript_seq[:ref_start_idx] +
        alt_cds_seq +
        transcript_seq[ref_end_idx:]
    )

    return new_transcript_seq


def analyze_transcript(results_df):
    valid_stop_codons = {"TAA", "TAG", "TGA"}
    start_codon = "ATG"

    df = results_df.copy()

    # Add new columns
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

        if not isinstance(seq, str) or len(seq) < 3:
            continue

        start_pos = None
        start_exon = None  # for exon number
        stop_codons = []
        stop_exons = []  # for exon number

        if not (row["start_loss"] or row["stop_loss"]):
            # If no start loss nor stop loss: skip the row and fill with None values
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

# TODO: test analysis of rules
def evaluate_nmd_escape_rules(row):

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

    # Long exon rule
    rule_long_exon = any(exon_length_map.get(exon, 0) > 400 for exon in stop_exons)

    # Start-proximal rule
    rule_start_proximal = start_pos is not None and stop_pos is not None and (stop_pos - start_pos) < 150

    # Escape if any rule is true
    escape = rule_last_exon or rule_50nt_penultimate or rule_long_exon or rule_start_proximal or rule_single_exon

    return {
        "nmd_last_exon_rule": rule_last_exon,
        "nmd_50nt_penultimate_rule": rule_50nt_penultimate,
        "nmd_long_exon_rule": rule_long_exon,
        "nmd_start_proximal_rule": rule_start_proximal,
        "nmd_single_exon_rule": rule_single_exon,
        "nmd_escape": escape
    }