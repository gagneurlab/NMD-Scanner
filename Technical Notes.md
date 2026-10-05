Summary of the NMD-Scanner Script:

1. Command Line argument Parser
2. Check that Output Path is valid
3. Reads genomic data files (VCF for variants, GFF3 for annotations, FASTA for sequences)
4. Extracts coding regions and exons + exon length -> scan.read_annotation(). The coding regions are CDS rows that include the stop codon, as a GFF3 CDS does -> scan.read_gff3()
5. Identifies premature termination codons (PTCs) —> extract_ptc()
    1. Check that the CDS rows are coding regions (column has_stop_codon)
    2. Intersect variants with CDS regions -> join_variants_to_cds()
    3. (in TCGA & MMRF only: adjust minus strand variants)
    4. Fetch reference CDS sequence for each variant region (on variant level: only CDS where a variant is located on) —> catch_sequence.add_exon_cds_sequence()
    5. Apply variant to CDS and compute alternative CDS sequence and get lengths of Ref-CDS and Alt-CDS (on variant level: only CDS where a variant is located on) —> apply_variant_edge_aware_with_lengths()
    6. Filter out variants with a reference mismatch and print those
    7. Limit to relevant transcripts (which are in the variant-CDS-intersection-DataFrame) for faster processing
    8. Fetch reference sequence for all CDS entries per (relevant) transcripts —> we get cds_df (annotation filtered for CDS) with Exon_CDS_seq in the end
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
- [x] extra_features.py


Output_features (TODO: need to revise this):

transcript_id, variant_id, chromosome, strand, ref, alt, start_variant, end_variant
————————————————————————————————————————————————————————————
for ref_cds and alt_cds:
- start, stop, seq, len
- info: List of tuples (exon_number, exon_lengths)
————————————————————————————————————————————————————————————
- cds_in_transcript: computed to check if the CDS sequence is in the transcript sequence
- has_stop_codon: whether the coding region ends in an annotated stop codon, i.e. whether the transcript has stop_codon rows. An Ensembl GFF3 has no stop_codon rows; there, the last 3 CDS bases in the FASTA decide -> scan.read_gff3(). Without an annotated stop codon (e.g. cds_end_NF), valid_stop and stop loss are False, every in-frame stop codon is premature, and utr3_length and stop_codon_distance are empty
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
- cds_start_in_transcript, cds_end_in_transcript: position of the ref CDS (with stop codon) in the transcript sequence, 0-based half-open
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
NMD rules (figures at the end of this file):
- Last exon rule: The PTC is in the last exon
- 50nt penultimate rule: The PTC is within 50 nucleotides upstream of the last exon junction
- Long exon rule: The PTC is in an exon with >407 nucleotides
- Start proximal rule: The PTC is within 150 nucleotides of the start codon
- Single exon rule: The transcript where the PTC lays consists only of a single exon
- NMD escape: A PTC is considered to escape NMD if it satisfies any of the above rules.
————————————————————————————————————————————————————————————
extra features (figures at the end of this file):
- utr3_length
- utr5_length
- total_exon_count
- upstream_exon_count
- downstream_exon_count
- ptc_to_start_codon
- ptc_less_than_150nt_to_start
- ptc_exon_length
- ptc_to_intron: distance in nt from the PTC to the 3' end of the PTC exon. For an internal exon, that end is the downstream exon junction. For the last exon, it is the transcript end, so the distance is the length of the 3' UTR that the PTC creates.

Figures of the features and the NMD rules:

The figures show a transcript 5' to 3' and are not to scale. `u` is UTR, `=` is CDS, `[...]` is an exon, `|` between two exons is an exon junction, and `*` is the PTC. `*--->|` is the distance from the PTC to an exon junction or to the transcript end, and `<--->` is a length. The numbers are positions in CDS coordinates, as alt_first_stop_pos: 0 is the first base of the start codon, and the 5' UTR has negative positions. After an indel, the positions are in alt CDS coordinates, so the exon junctions move with the PTC (see exon_end_in_alt_cds() in extra_features.py). The figures and the exon numbers follow the transcript. On the minus strand, the genomic coordinates run the other way.

Features of a PTC in exon 2 of 3. The CDS starts at transcript position 50 and has 360 nt with its stop codon, so the stop codon of the reference starts at 357. The figure gives each feature its value:

```
        exon 1: 150 nt       exon 2: 200 nt              exon 3: 250 nt
    5' [uuuuu==========]|[========*===========]|[============uuuuuuuuuuuuuuuuuuu] 3'
       -50   0          100       180          300          357                 550
       <----->  utr5_length = 50
             <-------------------->  ptc_to_start_codon = 180
                          <------------------>  ptc_exon_length = 200
                                  *----------->|  ptc_to_intron = 120
                                  <------------------------->  stop_codon_distance = 177
                                                             <------------------>  utr3_length = 190
```

The exon counts are total_exon_count = 3, upstream_exon_count = 1 and downstream_exon_count = 1. ptc_less_than_150nt_to_start is False, because ptc_to_start_codon is 180.

ptc_to_intron runs to the 3' end of the PTC exon. Where that end lies depends on the PTC exon:

PTC in an internal exon: to the downstream exon junction.

```
5' [uuu=====]|[=====*======]|[=========uuuuuuu] 3'
                    *------>|
```

PTC in the last CDS exon, followed by an exon with only 3' UTR: to the exon junction in the 3' UTR, not to the CDS end.

```
5' [uuu=====]|[============]|[====*===uuu]|[uuuuuuuuu] 3'
                                  *------>|
                                  *-->|  not to the CDS end
```

PTC in the last exon: to the transcript end. The value is the length of the 3' UTR that the PTC creates.

```
5' [uuu=====]|[============]|[====*======uuuuuuuuu] 3'
                                  *-------------->|
```

PTC in a single exon transcript: to the transcript end, as in the last exon.

```
5' [uuuu=========*=======uuuuuu] 3'
                 *------------>|
```

NMD rules. nmd_escape is True if one of the rules is True.

nmd_last_exon_rule: the PTC lies in the last exon, i.e. downstream_exon_count is 0. See the figure "PTC in the last exon" above.

nmd_single_exon_rule: the transcript has one exon, i.e. total_exon_count is 1. See the figure "PTC in a single exon transcript" above.

nmd_50nt_penultimate_rule: the PTC lies 1 to 50 nt upstream of the last exon junction, the 3' end of the penultimate exon.

```
5' [uuu=====]|[=======*====]|[=========uuuuuuu] 3'
                      *---->|  1 to 50 nt
```

If the last exon holds only 3' UTR, the last exon junction lies in the 3' UTR. The rule then measures to that junction, as in the figure "PTC in the last CDS exon" above.

nmd_long_exon_rule: the PTC exon has more than 407 nt, UTR included.

```
5' [uuu=====]|[=========*==============================]|[=====uuuuuu] 3'
               <-------------------------------------->  ptc_exon_length > 407
```

nmd_start_proximal_rule: the PTC lies less than 150 nt downstream of the start codon, i.e. alt_first_stop_pos - alt_start_codon_pos < 150.

```
5' [uuu=====*===]|[============]|[=========uuuuuuu] 3'
       <---->  < 150 nt
```
