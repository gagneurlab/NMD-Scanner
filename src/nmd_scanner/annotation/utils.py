# Import dependencies
import pandas as pd
from Bio.Seq import Seq
from tqdm import tqdm


# Functions used for extracting PTC:

def adjust_last_cds_for_stop_codon(df, transcript_col="transcript_id"):

    """
    Adjusts the genomic coordinates of the last CDS exon in each transcript by adding 3 positions,
    thus to include the stop codon.
    Plus strand: extend last exon (largest start position) at the END (+3 to End)
    Minus strand: extend last exon (smallest start position) at the START (-3 from Start)
    :param df: Dataframe containing CDS annotation
    :param exon_col: The name of the column that indicated the exon number, so we can find out which is the last CDS snippet
    :param transcript_col: The name of the column that indicates the transcript ID
    :return: Modified pandas DataFrame where the last exon of each transcript is extended by 3 bases to include the stop codon.
    """

    df = df.copy()

    adjusted_idx = []
    for tx, group in df.groupby(transcript_col):
        strand = group["Strand"].iloc[0]

        if strand == "+":
            # last exon has max Start position
            idx = group["Start"].idxmax()
            df.at[idx, "End"] += 3
        elif strand == "-":
            # last exon has min Start position
            idx = group["Start"].idxmin()
            df.at[idx, "Start"] -= 3

        adjusted_idx.append(idx)

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

    # Special handling for deletions (Ref = N and Alt = <DEL>)
    if ref == "N" and alt == "<DEL>":
        # Clip deletion to the CDS region (only remove overlap part)

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

        cds_index_start = overlap_start - cds_start
        cds_index_end = overlap_end - cds_start

        alt_seq = []
        alt_seq.extend(cds_seq[:cds_index_start])  # keep sequence before deletion
        alt_seq.extend(cds_seq[cds_index_end:])  # keep sequence after deletion

        return pd.Series({
            "Exon_CDS_length": len(cds_seq),
            "Exon_Alt_CDS_seq": "".join(alt_seq),
            "Exon_Alt_CDS_length": len(alt_seq)
        })

    # Special handling for duplications (Ref = N and Alt = <DUP>)
    if ref == "N" and alt == "<DUP>":

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

        cds_index_start = overlap_start - cds_start
        cds_index_end = overlap_end - cds_start

        alt_seq = []
        alt_seq.extend(cds_seq[:cds_index_start])  # sequence up to the end of the overlap
        alt_seq.extend(cds_seq[cds_index_start:cds_index_end]) # overlap-region (original, in CDS)
        alt_seq.extend(cds_seq[cds_index_start:cds_index_end])  # duplicate the overlapped region
        alt_seq.extend(cds_seq[cds_index_end:])  # keep sequence after duplication

        return pd.Series({
            "Exon_CDS_length": len(cds_seq),
            "Exon_Alt_CDS_seq": "".join(alt_seq),
            "Exon_Alt_CDS_length": len(alt_seq)
        })

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

    transcript_groups = list(intersection_cds_vcf.groupby("transcript_id"))
    for transcript_id, var_df in tqdm(transcript_groups, desc="Processing transcripts"):  # Only transcripts with a variant

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

        ## old option:
        # ref_cds_info = [
        #    (row["exon_number"], len(row["Exon_CDS_seq"]))
        #    for _, row in ref_exons.iterrows()
        # ]

        ## new option:
        ref_cds_info = sorted([
            (row["exon_number"], len(row["Exon_CDS_seq"]))
            for _, row in ref_exons.iterrows()
        ], key=lambda x: x[0])

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

            ## old option:
            # alt_cds_info = [
            #    (row["exon_number"], len(row["Exon_CDS_seq"]))
            #    for _, row in alt_exons.iterrows()
            # ]

            ## new option:
            alt_cds_info = sorted([
                (row["exon_number"], len(row["Exon_CDS_seq"]))
                for _, row in alt_exons.iterrows()
            ], key=lambda x: x[0])
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
            # Start at cds_start, scan codons in frame to end --> lets start at cds start so we are in frame
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

def calculate_utr_lengths(strand, ref_cds_info, transcript_exon_info, ref_cds_start, transcript_start, ref_cds_stop, transcript_end):

    ref_cds_info = ref_cds_info or []
    transcript_exon_info = transcript_exon_info or []

    if not ref_cds_info or not transcript_exon_info:
        return {"utr5_length": None, "utr3_length": None}

    # Convert to dicts for easier lookup
    transcript_exon_dict = {int(k): int(v) for k, v in transcript_exon_info}
    cds_exons_dict = {int(exon): int(length) for exon, length in ref_cds_info}

    # Handle single exon
    if len(transcript_exon_dict) == 1:
        if strand == "+":
            utr5 = ref_cds_start - transcript_start
            utr3 = transcript_end - ref_cds_stop
        else:
            utr5 = transcript_end - ref_cds_stop
            utr3 = ref_cds_start - transcript_start

        utr5 = utr5 if utr5 >= 0 else None
        utr3 = utr3 if utr3 >= 0 else None

        return {"utr5_length": utr5, "utr3_length": utr3}

    cds_exon_nums = sorted(cds_exons_dict.keys())
    tx_exon_nums = sorted(transcript_exon_dict.keys())

    utr5 = 0
    utr3 = 0

    for exon in tx_exon_nums:
        exon_len = transcript_exon_dict[exon]
        cds_len = cds_exons_dict.get(exon, 0)
        utr_len = exon_len - cds_len

        if utr_len <= 0:
            continue

        if exon in cds_exons_dict:
            # Exon overlaps CDS, partial UTR
            if strand == "+":
                if exon == cds_exon_nums[0]:
                    utr5 += utr_len
                elif exon == cds_exon_nums[-1]:
                    utr3 += utr_len
            else:
                if exon == cds_exon_nums[0]:
                    utr3 += utr_len
                elif exon == cds_exon_nums[-1]:
                    utr5 += utr_len
        else:
            # Exon is outside CDS
            if strand == "+":
                if exon < cds_exon_nums[0]:
                    utr5 += exon_len
                elif exon > cds_exon_nums[-1]:
                    utr3 += exon_len
            else:
                if exon > cds_exon_nums[-1]:
                    utr5 += exon_len
                elif exon < cds_exon_nums[0]:
                    utr3 += exon_len

    utr5 = utr5 if utr5 >= 0 else None
    utr3 = utr3 if utr3 >= 0 else None
    return {"utr5_length": utr5, "utr3_length": utr3}

def calculate_exon_features(transcript_exon_info, alt_stop_codon_exons, alt_is_premature):

    """
    Calculate exon-related features:
    - total_exon_count: always computed if transcript_exon_info is available
    - upstream_exon_count / downstream_exon_count: only computed if a PTC exists
    """

    exon_info = transcript_exon_info or []
    stop_exons = alt_stop_codon_exons or []

    total_exons = len(exon_info)

    if not alt_is_premature or not stop_exons or not exon_info:
        return {
            "total_exon_count": total_exons if total_exons > 0 else None,
            "upstream_exon_count": None,
            "downstream_exon_count": None
        }

    # Take the PTC exon closest to CDS start
    ptc_exon = min(int(e) for e in stop_exons)

    # get exon numbers from transcript_exon_info
    exon_numbers = [int(e[0]) for e in exon_info]

    # If the PTC exon is not in transcript → cannot compute
    if ptc_exon not in exon_numbers:
        return {
            "total_exon_count": int(total_exons),
            "upstream_exon_count": None,
            "downstream_exon_count": None
        }

    upstream = sum(1 for e in exon_numbers if e < ptc_exon)
    downstream = sum(1 for e in exon_numbers if e > ptc_exon)

    return {
        "total_exon_count": int(total_exons),
        "upstream_exon_count": int(upstream),
        "downstream_exon_count": int(downstream)
    }

def calculate_ptc_to_start_distance(alt_is_premature, alt_start_codon_pos, alt_first_stop_pos):

    if not alt_is_premature:
        return None

    start = alt_start_codon_pos
    stop = alt_first_stop_pos

    if start is None or stop is None:
        return None

    # PTC codon position: is PTC_to_start_codon / 3 --> leave it out
    #offset = stop - start
    #return offset // 3 if offset >= 0 else None

    if stop <= start:
        return None

    return stop-start # distance between the PTC to start codon in nt

def calculate_ptc_exon_length(alt_is_premature, alt_stop_codon_exons, transcript_exon_info):
    """
    Return the length of the exon containing the first premature stop codon (PTC).
    """

    if not alt_is_premature:
        return None

    stop_exons = alt_stop_codon_exons or []
    exon_info = transcript_exon_info or []

    if not stop_exons or not exon_info:
        return None

    # First PTC exon = smallest exon number (transcript-order, strand-corrected)
    ptc_exon = min(int(e) for e in stop_exons)

    exon_dict = {int(e): int(length) for e, length in exon_info}
    return exon_dict.get(ptc_exon)

def calculate_stop_codon_dist(ref_first_stop_pos, alt_first_stop_pos):

    """
    Calculate the distance between the reference stop codon and the alternative stop codon.
    Positive means the PTC is upstream of the reference stop codon.
    """

    ref_stop = ref_first_stop_pos
    alt_stop = alt_first_stop_pos

    if ref_stop is None or alt_stop is None:
        return None

    return ref_stop - alt_stop

def evaluate_nmd_escape_rules(alt_is_premature, alt_first_stop_pos, alt_stop_codon_exons, 
                              transcript_exon_info, alt_start_codon_pos, total_exon_count, 
                              downstream_exon_count, ptc_exon_length):

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

    :return: A dictionary with boolean flags for each rule and overall NMD escape
    """

    # Only relevant for premature stop codons
    if not alt_is_premature:
        return {
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": False
        }

    # Extract relevant data
    stop_pos = alt_first_stop_pos
    exon_info = transcript_exon_info or []
    start_pos = alt_start_codon_pos

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

    return {
        "nmd_last_exon_rule": rule_last_exon,
        "nmd_50nt_penultimate_rule": rule_50nt_penultimate,
        "nmd_long_exon_rule": rule_long_exon,
        "nmd_start_proximal_rule": rule_start_proximal,
        "nmd_single_exon_rule": rule_single_exon,
        "nmd_escape": escape
    }

def calculate_ptc_to_downstream_ej(alt_is_premature, alt_stop_codon_exons, alt_cds_info, alt_first_stop_pos):
    """
    Calculate distance from PTC to the downstream exon junction (next exon start/end depending on strand).
    Returns None if not applicable.
    """

    # only calculate if we have PTC
    if not alt_is_premature:
        return None

    stop_exons = alt_stop_codon_exons or []

    # Assuming that the PTC cannot be outside the CDS, since it needs to come before the original stop codon
    exon_info = alt_cds_info or []

    ptc_pos = alt_first_stop_pos

    if not stop_exons or not exon_info or ptc_pos is None:
        return None

    # Choose the PTC exon (smallest number, closer to start)
    ptc_exon = min(stop_exons)

    # Sum lengths of exons up to and including ptc_exon
    cumulative_length = 0
    for exon_num, length in exon_info:
        cumulative_length += length
        if exon_num == ptc_exon:
            break

    # Distance from PTC to downstream exon junction
    distance = cumulative_length - ptc_pos
    return distance

def add_likely_misannotated_flag(cds_in_transcript, ref_start_codon_pos, ref_valid_stop):

    """
    Flag rows that look inconsistent between CDS and transcript annotations and might be likely misannotated.
    A row is flagged as likely misannotated if any of these conditions apply:
        cds_in_transcript = False (the assembled CDS is not found in the transcript sequence)
        ref_start_codon_pos is defined and not 0 (reference CDS has a start codon not at the very start)
        ref_valid_stop is False (the last reference codon is not a valid stop codon)

    :return: A boolean flag. True if any condition above is met and thus the row is likely misannotated, False otherwise.
    """

    # Add likely_misannotated flag: when
    # "cds_in_transcript" is FALSE
    # "ref_start_codon_pos" is not 0
    # "ref_valid_stop" is FALSE

    # if any of these are missing entirely, flag as likely misannotated
    if cds_in_transcript is None or ref_start_codon_pos is None or ref_valid_stop is None:
        return True

    flag = (
        (cds_in_transcript is False) or
        ((ref_start_codon_pos is not None) and (ref_start_codon_pos != 0)) or
        (ref_valid_stop is False)
    )

    return(flag)

def evaluate_nmd_escape_rules(alt_is_premature, alt_first_stop_pos, alt_stop_codon_exons, 
                              transcript_exon_info, alt_start_codon_pos, total_exon_count, 
                              downstream_exon_count, ptc_exon_length):

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

    :return: A dictionary with boolean flags for each rule and overall NMD escape
    """

    # Only relevant for premature stop codons
    if not alt_is_premature:
        return {
            "nmd_last_exon_rule": False,
            "nmd_50nt_penultimate_rule": False,
            "nmd_long_exon_rule": False,
            "nmd_start_proximal_rule": False,
            "nmd_single_exon_rule": False,
            "nmd_escape": False
        }

    # Extract relevant data
    stop_pos = alt_first_stop_pos
    exon_info = transcript_exon_info or []
    start_pos = alt_start_codon_pos

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

    return {
        "nmd_last_exon_rule": rule_last_exon,
        "nmd_50nt_penultimate_rule": rule_50nt_penultimate,
        "nmd_long_exon_rule": rule_long_exon,
        "nmd_start_proximal_rule": rule_start_proximal,
        "nmd_single_exon_rule": rule_single_exon,
        "nmd_escape": escape
    }