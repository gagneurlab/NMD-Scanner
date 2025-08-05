Summary of the NMD-Scanner Script:

1. Command Line argument Parser
2. Check that Output Path is valid
3. Reads genomic data files (VCF for variants, GTF for annotations, FASTA for sequences)
4. Extracts coding sequences and exons + exon length
5. Identifies premature termination codons (PTCs) —> extract_ptc()
    1. Adjust the last 3 CDS positions to include stop codons in the sequence —> adjust_last_cds_for_stop_codon()
    2. Intersect variants with CDS regions
    3. (in TCGA & MMRF only: adjust minus strand variants)
    4. Fetch reference CDS sequence for each variant region (on variant level: only CDS where a variant is located on) —> catch_sequence.add_exon_cds_sequence()
    5. Apply variant to CDS and compute alternative CDS sequence and get lengths of Ref-CDS and Alt-CDS (on variant level: only CDS where a variant is located on) —> apply_variant_edge_aware_with_lengths()
    6. Filter out variants with a reference mismatch and print those
    7. Limit to relevant transcripts (which are in the variant-CDS-intersection-DataFrame) for faster processing
    8. Fetch reference sequence for all CDS entries per (relevant) transcripts —> we get cds_df (GTF filtered for CDS) with Exon_CDS_seq in the end
    9. create_reference_cds()
        1. Get full reference CDS per transcript by stitching exon CDS regions together plus length
        2. Get full alternative CDS with length
        3. get CDS exon information (exon number, CDS exon length) for ref and alt
    10. Get whole transcript sequences plus information of relevant transcripts —> get_transcript_sequence()
    11. Validate that the CDS sequence we generated before is present inside the transcript sequence generated in step 10, to make sure the transcript sequence was computed correctly
    12. Analyze reference and alternative CDS for start and stop codons, their positions, and potential premature termination codons (PTCs) —> analyze_sequence()
    13. Compare reference to alternative sequence start and stop codons and check if the codons were lost in the alternative sequence —> start_stop_loss()
    14. In case of start or stop loss:
        1. Annotate transcript information (transcript start, end, sequence, length, exon information)
        2. create alternative transcript sequence + length by exchanging the reference CDS by the alternative CDS—> splice_alt_cds_into_transcript()
        3. Add transcript exon information to the dataframe
        4. Analyze the transcript sequence for length, start / stop codon positions, etc., basically same analysis as we did for reference and alternative CDS sequence in step 12) —> analyze_transcript()
    15. Return to cli.py and call evaluate_nmd_escape_rules(): Evaluates whether a premature stop codon in a transcript is likely to escape nonsense-mediated decay (NMD) based on established biological rules. This function applies five NMD escape rules to determine if a premature termination codon (PTC) is likely to escape degradation. Returns dictionary
    16. Join Dictionary containing NMD rules with results from extract_ptc() and save output
    17. Compute some additional features such as UTR lengths, total / downstream / upstream exon count, and other positional information of the PTC —> extra_features.py : add_nmd_features()
    18. Join Output with additional features with our Original Dataframe (summarizing all annotated variants) and save output


Output files:

1_variant_exon_output.tsv: exon variant merge result, saved after step 5.6
2_cds_df_adj.tsv: reference sequence for entire CDS per (relevant) transcripts, saved after step 5.8
3_create_reference_CDS.tsv: full ref and alt CDS sequence plus length and CDS exon information, saved after step 5.9
4_transcript_sequences.tsv: full transcript sequences for relevant transcript plus start, end, strand, transcript length, transcript exon information, saved after step 10
5_final_ptc_analysis.tsv: Dataframe with ref & alt & transcript sequence with length / start / end / exon information / start & stop codon information etc., saved after step 14.4
6_nmd_rules.tsv: File with all information for ref and alt CDS sequences and transcript sequences per variant + NMD rules in case of PTC, saved in step 16
final_nmd_results.csv: File with all features, saved in step 18

- [x] cli.py
- [x] scan.py
- [x] rules.py
- [x] catch_sequence.py
- [ ] analyze_gtf.py (where called?) —> not used
- [x] extra_features.py


Output_features (TODO: need to revise this):

transcript_id, variant_id, chromosome, strand, ref, alt, start_variant, end_variant
————————————————————————————————————————————————————————————
for ref_cds and alt_cds:
- start, stop, seq, len
- info: List of tuples (exon_number, exon_lengths)
————————————————————————————————————————————————————————————
- cds_in_transcript: computed to check if the CDS sequence is in the transcript sequence
————————————————————————————————————————————————————————————
analyzing ref_ and alt_  (CDS)
- start_codon_pos
- start_codon_exon
- last_codon
- valid_stop
- first_stop_codon
- first_stop_pos
- num_stop_codons
- all_stop_codons
- stop_codon_exons
- is_premature
————————————————————————————————————————————————————————————
- start loss
- stop loss
————————————————————————————————————————————————————————————
transcript:
- start, end, seq, len
- alt_transcript_seq
- alt_transcript_length
- transcript_:
    - exon_info
    - start_codon_pos
    - start_codon_exon
    - last_codon
    - valid_stop
    - first_stop_codon
    - first_stop_pos
    - num_stop_codons
    - all_stop_codons
    - stop_codon_exons
————————————————————————————————————————————————————————————
NMD rules:
- Last exon rule: The PTC is in the last exon
- 50nt penultimate rule: The PTC is within 50 nucleotides upstream of the last exon junction
- Long exon rule: The PTC is in an exon with >407 nucleotides
- Start proximal rule: The PTC is within 150 nucleotides of the start codon
- Single exon rule: The transcript where the PTC lays consists only of a single exon
- NMD escape: A PTC is considered to escape NMD if it satisfies any of the above rules.
————————————————————————————————————————————————————————————
extra features:
- utr3_length
- utr5_length
- total_exon_count
- upstream_exon_count
- downstream_exon_count
- ptc_to_start_codon
- ptc_less_than_150nt_to_start
- ptc_exon_length
