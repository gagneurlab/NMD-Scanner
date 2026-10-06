Summary of the NMD-Scanner Script:

1. Command Line argument Parser
2. Check that Output Path is valid
3. Reads genomic data files (VCF for variants, GFF3 for annotations, FASTA for sequences)
4. Extracts coding regions and exons + exon length -> scan.read_annotation(). The coding regions are CDS rows that include the stop codon, as a GFF3 CDS does. Their columns has_start_codon and has_stop_codon say whether the transcript has an annotated start and stop codon -> scan.read_gff3()
5. Identifies premature termination codons (PTCs) —> extract_ptc()
    1. Check that the CDS rows are coding regions (columns has_start_codon and has_stop_codon)
    2. Find the equivalent placements of each variant and check its whole REF against the reference genome -> place_variants(). Join the placements with the CDS regions and the splice dinucleotides at their exon edges -> join_variant_windows(). The section "Variants at exon boundaries" below explains the placements.
    3. (in TCGA & MMRF only: adjust minus strand variants)
    4. Filter out variants with a reference mismatch and print those
    5. Apply each variant to the coding rows of each transcript, and compute the alternative sequence of each coding row and the change of the UTRs next to the coding region -> apply_variants(), variant_placement.place_in_transcript(). If the variant destroys a splice site or leaves an exon boundary ambiguous, the alternative transcript is unknown. The row stays, with unknown_reason set and without an alternative sequence.
    6. Drop the variant-transcript pairs that touch neither a coding row nor the splice dinucleotides at its exon edges
    7. Limit to relevant transcripts (which are in the variant-CDS-intersection-DataFrame) for faster processing
    8. Fetch reference sequence for all CDS entries per (relevant) transcripts —> we get cds_df (annotation filtered for CDS) with Exon_CDS_seq in the end
    9. create_reference_cds()
        1. Get full reference CDS per transcript by stitching exon CDS regions together plus length
        2. Get full alternative CDS with length
        3. get CDS exon information (exon number, CDS exon length) for ref and alt
        4. take the GFF3 phase (column Frame) of the 5'-most CDS row as cds_frame, the start of the first complete codon
    10. Get whole transcript sequences plus information of relevant transcripts —> get_transcript_sequence()
    11. Validate that the CDS sequence we generated before is present inside the transcript sequence generated in step 10, to make sure the transcript sequence was computed correctly
    12. Analyze reference and alternative CDS for start and stop codons, their positions, and potential premature termination codons (PTCs) —> analyze_sequence()
    13. Compare reference to alternative sequence start and stop codons and check if the codons were lost in the alternative sequence —> start_stop_loss()
    14. In case of start or stop loss:
        1. Annotate transcript information (transcript start, end, sequence, length, exon information)
        2. create alternative transcript sequence + length by exchanging the reference CDS by the alternative CDS, and the UTR bases next to it that the variant changes -> splice_alt_cds_into_transcript()
        3. Add transcript exon information to the dataframe
        4. Analyze the transcript sequence for length, start / stop codon positions, etc., basically same analysis as we did for reference and alternative CDS sequence in step 12) —> analyze_transcript()
        5. For every row, classify the first in-frame stop codon of the alternative transcript, read from the first complete codon of the CDS into the 3' UTR -> analyze_transcript(). The result replaces alt_is_premature and stop_loss of steps 12 and 13, before the scan of step 4 reads stop_loss. The section "Stop codon classification" below gives the cases.
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


## Output columns

Each row of the result is one variant in one transcript whose coding region the variant touches. The coding region is the CDS plus the stop codon. A variant touches it if one of its placements changes a base of the coding region, or of a splice dinucleotide at one of its exon edges (see [Variants at exon boundaries](#variants-at-exon-boundaries)). A variant gives no row if it touches no coding region, or if its ALT is a symbolic allele or a breakend (`rules.drop_symbolic_alleles`). It gives no row either if its REF does not match the reference genome.

The tables below list the 82 columns in output order, as `OUTPUT_COLUMN_KINDS` in `nmd_scanner.schema` does. There is one table for each of its three parts: the PTC columns, the NMD features and the NMD rules.

### Kinds and dtypes

The kind of a column sets its pandas dtype (`schema.KIND_DTYPES`) and its Parquet type (`cli.parquet_schema`):

| Kind | pandas dtype | Parquet type |
|---|---|---|
| int | `Int64` | `int64` |
| bool | `boolean` | `bool` |
| string | `string`, with python storage | `string` |
| pair_list | `object`: a list of (exon_number, length) tuples | `list<list<int64>>`: each inner list is [exon_number, length] |
| int_list | `object`: a list of exon numbers | `list<int64>` |
| stop_codon_list | `object`: a list of (position, codon) tuples, e.g. (5442, "TGA") | `list<struct<position: int64, codon: string>>` |

A null is pd.NA in an int, bool or string column. In a list column, it is None or NaN, so test for it with `pd.isna`. Parquet stores a null as null, and CSV writes it as an empty field.

### Positions and terms

- Genomic positions are 0-based half-open: a start is the first base, and an end is the base after the last one. They do not depend on the strand, so on the minus strand a start is the 3' end.
- CDS positions count from the 5' base of the ref or alt CDS, which is position 0. Transcript positions count the same way from the 5' base of the transcript sequence. Both run 5' to 3', also on the minus strand.
- In-frame means in the frame of the first complete codon of the CDS, at CDS position `cds_frame`. `cds_frame` is the GFF3 phase of the 5'-most CDS row. It is 1 or 2 only for a CDS that lacks its 5' end, e.g. of a `cds_start_NF` transcript, whose first 1 or 2 bases belong to no complete codon. Every codon scan starts at the first complete codon: the scan of the ref and the alt CDS, the stop codon classification, and the scan of the alt transcript.
- The start codon (`ref_start_codon_pos`, `alt_start_codon_pos`) is the annotated start codon at CDS position 0 (`has_start_codon`). It can be a non-ATG codon such as CTG, and an in-frame ATG downstream of it is an internal Met. A `cds_start_NF` transcript has no annotated start codon: its true start lies upstream of the CDS, at an unknown distance. For the rows of such a transcript, both positions and `ptc_to_start_codon` are null, `start_loss` and `nmd_start_proximal_rule` are False, and `likely_misannotated` is True.
- The codon scans know only the stop codons TAA, TAG and TGA, also on the mitochondrial chromosome.
- A PTC row is a row with `alt_is_premature` True. Its PTC is the first in-frame stop codon of the alt CDS, at `alt_first_stop_pos`. The PTC exon is the exon that holds the PTC.

### Variants at exon boundaries

A VCF record names one placement of a variant. An indel in a repeat has more placements that give the same alt sequence, because it shifts through the repeat. A delins whose REF and ALT differ in length has two placements: its bases matched from the left, and from the right. An SNV or MNV has one placement. `rules.place_variants` and `variant_placement.place_in_transcript` use all placements of a variant.

A splice dinucleotide is the pair of intron bases next to an exon edge: the donor after an exon, and the acceptor before it. A transcript end has none. Each placement maps each exon edge base for base to a position in the alt sequence. So a splice dinucleotide that the ALT bases form elsewhere does not count. A placement is valid if the reference splice dinucleotide stays next to each of its positions. If all valid placements agree, the exon edges lie there. Otherwise the alt transcript is unknown, and `unknown_reason` says why:

- `splice_site_destroyed`: no placement is valid. E.g. an SNV in a splice dinucleotide, or a deletion over an exon edge that leaves no GT after the exon. It is also the reason if the edge rules of the coding region below take the two edges of a coding row from different placements, so that the row would end before it starts.
- `exon_boundary_ambiguous`: the valid placements put an exon edge at different positions. E.g. AG inserted at CAG|AG, where the exon can start at either AG. A transcript end has no splice dinucleotide, so placements that put it at different positions make it ambiguous too. An SNV or MNV has one placement, so it is never ambiguous.

An edge of the coding region inside an exon has no splice dinucleotide. There, the valid placements decide:

- At the stop codon, the edge takes the placement shifted farthest into the 3' UTR. So an indel that fits into the 3' UTR leaves the CDS unchanged, e.g. an insertion right after the stop codon. Such a variant gives no row, unless another placement reaches into the coding region. Then its row has an unchanged CDS.
- At the start codon, the edge takes the 5'-most position, over all placements, where an ATG starts in the alt sequence, because a scanning ribosome starts at the first ATG. Without such a position, it takes the placement shifted farthest into the 5' UTR. So an insertion right before the start codon changes the CDS only if its bases start with ATG.
- Of the two placements of a delins, the one with its length change at the UTR end counts as shifted into the UTR.

Where an edge of the coding region is also an exon edge, the splice site decides which bases enter the mRNA, and these rules decide which of them are coding. The alt bases that these rules leave outside the coding region go into the UTRs of `alt_transcript_seq`, and the deleted UTR bases leave it. So a length change of the 5' UTR moves the alt CDS in the alt transcript. `alt_cds_start_in_transcript` gives its position.

An unknown alt transcript has no alt CDS. Its row is null in the alt CDS columns (`alt_cds_start` to `alt_cds_len`, `alt_cds_info`), in the alt codon columns (`alt_start_codon_pos` to `alt_is_premature`), in `start_loss`, `stop_loss`, `alt_transcript_seq`, `alt_transcript_length`, `alt_cds_start_in_transcript` and the 9 scan columns of the alt transcript. The NMD features are null too, except `utr3_length`, `utr5_length`, `total_exon_count` and `likely_misannotated`, which describe the ref transcript. Every NMD rule is null.

### Stop codon classification

`alt_is_premature` and `stop_loss` come from the first in-frame stop codon of `alt_transcript_seq`, read from the first complete codon of the alt CDS, at `alt_cds_start_in_transcript` + `cds_frame`, through the CDS into the 3' UTR (`rules.analyze_transcript`). The codon scan of the alt CDS misses this stop codon after a frameshift: the last 3 nt of the alt CDS are then out of frame, and the stop codon can lie in the 3' UTR. The first stop codon is compared with the annotated stop codon, both at their positions in `alt_transcript_seq`:

- Upstream of the annotated stop codon: a PTC. `alt_is_premature` is True, and `stop_loss` is False.
- At the annotated stop codon: neither. Both are False. An insertion inside the stop codon, e.g. TAA>TGAA, is such a case.
- Downstream of the annotated stop codon, or no stop codon up to the transcript end (nonstop): a stop loss. `alt_is_premature` is False, and `stop_loss` is True.

`rules.annotated_stop_in_alt` locates the annotated stop codon in `alt_transcript_seq` from the longest common prefix and suffix of `transcript_seq` and `alt_transcript_seq`. So the result does not depend on how the VCF describes the variant. An indel in a repeat that reaches into the stop codon has two placements that matter. It is placed 3'-most, as HGVS does, and also 5'-most if that placement ends at or before the stop codon. A first stop codon at either position is the annotated stop codon. A variant that removes the start of the stop codon maps it to the position anchored on the unchanged sequence downstream. So the annotated stop codon counts as lost only if every placement of the variant removes it. A delins is not split into an SNV plus an indel: ATAG>T at the last sense codon and the stop codon of GTA TAG TAG CAT gives GTT TAG CAT. That is a stop loss, although the protein is unchanged.

Without an annotated stop codon (`has_stop_codon` False), a first stop codon inside the alt CDS is a PTC, and one past its end is neither. `stop_loss` is then False.

Two kinds of rows keep the flags from the CDS: `alt_is_premature` from the codon scan of the alt CDS, and `stop_loss` from `ref_valid_stop` and `alt_valid_stop`.

- Rows where `alt_transcript_seq` or `cds_start_in_transcript` is null, except rows with `unknown_reason`. On those, both flags are null.
- Rows whose `transcript_seq`, read in frame from `cds_start_in_transcript` + `cds_frame`, does not read its first stop codon at the annotated stop codon. The annotated stop codon is then out of frame, e.g. because the CDS length does not fit `cds_frame`, or it is no stop codon, or an in-frame stop codon lies upstream of it. That upstream stop codon is a selenocysteine TGA or a misannotation. The pipeline does not read selenocysteine annotations, so it cannot tell the two apart.

### PTC columns

`extract_ptc` returns these 65 columns. A codon scan of the ref and alt CDS gives the `ref_*` and `alt_*` codon columns. `alt_is_premature` and `stop_loss` come from the alt transcript instead (see [Stop codon classification](#stop-codon-classification)).

A second scan, of `alt_transcript_seq`, gives the last 9 columns, from `transcript_start_codon_pos` on. It runs only if `start_loss` or `stop_loss` is True and `alt_transcript_seq` is not null. Otherwise the 9 columns are null, and the table says "not scanned". After a start loss, the scan takes the first ATG at or after the scan start, in any frame, and reads the stop codons in the frame of that ATG. After a stop loss without a start loss, it reads the codons in the frame of the scan start. The scan start is the first complete codon of the alt CDS in the alt transcript, at `alt_cds_start_in_transcript` + `cds_frame`. Since `stop_loss` comes from the alt transcript, the scan also runs after a frameshift whose first stop codon lies in the 3' UTR. After a stop loss without a start loss, `transcript_first_stop_pos` is the first in-frame stop codon from the scan start. Unless the row keeps the flags from the CDS, that is the new stop codon downstream of the annotated one, or null for a nonstop.

| Column | Kind | Meaning | Null when |
|---|---|---|---|
| `transcript_id` | string | Transcript ID from the annotation | never |
| `variant_id` | string | ID of the VCF record, `.` if it has none | never |
| `ref_cds_start` | int | Genomic start of the coding region: the smallest Start of its CDS rows | never |
| `ref_cds_stop` | int | Genomic end of the coding region: the largest End of its CDS rows | never |
| `ref_cds_seq` | string | Sequence of the ref CDS, 5' to 3', stop codon included | never |
| `ref_cds_len` | int | Length of `ref_cds_seq` | never |
| `alt_cds_start` | int | Same as `ref_cds_start`: the variant does not move the bounds of the coding region | `unknown_reason` is set |
| `alt_cds_stop` | int | Same as `ref_cds_stop` | `unknown_reason` is set |
| `alt_cds_seq` | string | Sequence of the alt CDS: `ref_cds_seq` with the variant applied. At an edge of the coding region, the placements of the variant decide which alt bases are coding (see [Variants at exon boundaries](#variants-at-exon-boundaries)) | `unknown_reason` is set |
| `alt_cds_len` | int | Length of `alt_cds_seq` | `unknown_reason` is set |
| `chromosome` | string | Chromosome name, as in the input files | never |
| `gene_id` | string | Gene ID from the annotation | never |
| `strand` | string | Strand of the transcript, `+` or `-` | never |
| `has_start_codon` | bool | Whether the CDS starts with an annotated start codon, at CDS position 0. It can be a non-ATG codon such as CTG. A GENCODE GFF3 marks it with `start_codon` rows, and a `cds_start_NF` transcript has none. An Ensembl GFF3 has no `start_codon` rows: a CDS has a start codon if it starts with ATG and its 5'-most CDS row has phase 0 (`scan.read_gff3`) | never |
| `has_stop_codon` | bool | Whether the coding region ends in an annotated stop codon. A GENCODE GFF3 marks it with `stop_codon` rows. In an Ensembl GFF3, the last 3 CDS bases in the FASTA decide (`scan.read_gff3`) | never |
| `cds_frame` | int | GFF3 phase (column 8) of the 5'-most CDS row: the number of bases before the first complete codon, 0 to 2. It is 1 or 2 only for a CDS that lacks its 5' end, e.g. of a `cds_start_NF` transcript. Every codon scan starts at CDS position `cds_frame` (see [Positions and terms](#positions-and-terms)) | never |
| `ref` | string | REF allele of the VCF record | never |
| `alt` | string | ALT allele of the VCF record | never |
| `start_variant` | int | Genomic start of the variant: VCF POS minus 1 | never |
| `end_variant` | int | Genomic end of the variant: `start_variant` plus the length of REF | never |
| `ref_cds_info` | pair_list | (exon_number, length) of the ref CDS part in each exon, in exon number order | never |
| `alt_cds_info` | pair_list | (exon_number, length) of the alt CDS part in each exon, in exon number order | `unknown_reason` is set |
| `cds_in_transcript` | bool | Whether `ref_cds_seq` occurs in `transcript_seq`. False if the transcript has no exon rows | never |
| `ref_start_codon_pos` | int | CDS position of the annotated start codon of the ref CDS: 0 | `has_start_codon` is False; the CDS has fewer than 3 nt |
| `ref_start_codon_exon` | int | Exon number of `ref_start_codon_pos` | as `ref_start_codon_pos` |
| `ref_last_codon` | string | Last 3 nt of the ref CDS | the CDS has fewer than 3 nt |
| `ref_valid_stop` | bool | Whether `has_stop_codon` is True and `ref_last_codon` is a stop codon | the CDS has fewer than 3 nt |
| `ref_first_stop_codon` | string | First in-frame stop codon of the ref CDS | no in-frame stop codon; the CDS has fewer than 3 nt |
| `ref_first_stop_pos` | int | CDS position of `ref_first_stop_codon` | as `ref_first_stop_codon` |
| `ref_num_stop_codons` | int | Number of in-frame stop codons in the ref CDS | the CDS has fewer than 3 nt |
| `ref_all_stop_codons` | stop_codon_list | (CDS position, codon) of each in-frame stop codon | the CDS has fewer than 3 nt |
| `ref_stop_codon_exons` | int_list | Exon number of each in-frame stop codon, in the order of `ref_all_stop_codons` | the CDS has fewer than 3 nt |
| `ref_is_premature` | bool | Whether the first in-frame stop codon starts before the last 3 nt of the CDS, which hold the annotated stop codon. If `has_stop_codon` is False, any in-frame stop codon is premature. False without an in-frame stop codon | the CDS has fewer than 3 nt |
| `alt_start_codon_pos` | int | CDS position of the annotated start codon in the alt CDS: 0 if the alt CDS starts with the same 3 nt as the ref CDS | `unknown_reason` is set; as `ref_start_codon_pos`; the variant changes the start codon |
| `alt_start_codon_exon` | int | As `ref_start_codon_exon`, for the alt CDS | as `alt_start_codon_pos` |
| `alt_last_codon` | string | As `ref_last_codon`, for the alt CDS | `unknown_reason` is set; as `ref_last_codon` |
| `alt_valid_stop` | bool | As `ref_valid_stop`, for the alt CDS | `unknown_reason` is set; as `ref_valid_stop` |
| `alt_first_stop_codon` | string | As `ref_first_stop_codon`, for the alt CDS. On a PTC row, this is the PTC | `unknown_reason` is set; as `ref_first_stop_codon` |
| `alt_first_stop_pos` | int | As `ref_first_stop_pos`, for the alt CDS | `unknown_reason` is set; as `ref_first_stop_codon` |
| `alt_num_stop_codons` | int | As `ref_num_stop_codons`, for the alt CDS | `unknown_reason` is set; as `ref_num_stop_codons` |
| `alt_all_stop_codons` | stop_codon_list | As `ref_all_stop_codons`, for the alt CDS | `unknown_reason` is set; as `ref_all_stop_codons` |
| `alt_stop_codon_exons` | int_list | As `ref_stop_codon_exons`, for the alt CDS | `unknown_reason` is set; as `ref_stop_codon_exons` |
| `alt_is_premature` | bool | Whether the first in-frame stop codon of the alt transcript lies upstream of the annotated stop codon (see [Stop codon classification](#stop-codon-classification)). Without an annotated stop codon: whether it lies inside the alt CDS. True marks a PTC row. On a row that keeps the flags from the CDS: as `ref_is_premature`, for the alt CDS | `unknown_reason` is set; the row keeps the flags from the CDS, and the alt CDS has fewer than 3 nt |
| `start_loss` | bool | Whether the variant changes the annotated start codon: `ref_start_codon_pos` is 0, and `alt_start_codon_pos` is null. The start codon can be a non-ATG codon such as CTG. A transcript without an annotated start codon, e.g. a `cds_start_NF` one, has no start codon to lose, so `start_loss` is False there | `unknown_reason` is set |
| `stop_loss` | bool | Whether the first in-frame stop codon of the alt transcript lies downstream of the annotated stop codon, or the alt transcript has none (see [Stop codon classification](#stop-codon-classification)). A change to another stop codon, e.g. TAA>TAG, is no loss. False without an annotated stop codon. On a row that keeps the flags from the CDS: whether `ref_valid_stop` is True and `alt_valid_stop` is not | `unknown_reason` is set |
| `transcript_start` | int | Genomic start of the transcript: the smallest Start of its exon rows | the transcript has no exon rows |
| `transcript_end` | int | Genomic end of the transcript: the largest End of its exon rows | as `transcript_start` |
| `transcript_seq` | string | Sequence of the ref transcript: its exons, spliced, 5' to 3' | as `transcript_start` |
| `transcript_length` | int | Length of `transcript_seq` | as `transcript_start` |
| `cds_start_in_transcript` | int | Transcript position of the 5' CDS base | the transcript has no exon rows; the 5' CDS base lies outside the exons |
| `cds_end_in_transcript` | int | Transcript position after the 3' CDS base. The coding region ends there, after the stop codon | as `cds_start_in_transcript` |
| `alt_transcript_seq` | string | `transcript_seq` with `ref_cds_seq` replaced by `alt_cds_seq`, and with the change of the UTRs next to it: the alt bases that the edge rules leave outside the coding region, and the deleted UTR bases (see [Variants at exon boundaries](#variants-at-exon-boundaries)) | `unknown_reason` is set; `cds_start_in_transcript` is null; `transcript_seq` does not hold `ref_cds_seq` at `cds_start_in_transcript`, or the ref bases of a UTR change next to it |
| `alt_transcript_length` | int | Length of `alt_transcript_seq` | as `alt_transcript_seq` |
| `alt_cds_start_in_transcript` | int | Position of the 5' base of the alt CDS in `alt_transcript_seq`. It differs from `cds_start_in_transcript` if the variant changes the length of the 5' UTR. The stop codon classification and the scan of the alt transcript read from here | as `alt_transcript_seq` |
| `transcript_exon_info` | pair_list | (exon_number, length) of each exon of the transcript, 5' to 3' | the transcript has no exon rows |
| `transcript_start_codon_pos` | int | Transcript position of the start codon of the scan. After a start loss: the ATG that the scan found. After a stop loss without a start loss: the annotated start codon, at `alt_cds_start_in_transcript` | not scanned; after a start loss: no ATG found; after a stop loss: `has_start_codon` is False |
| `transcript_start_codon_exon` | int | Exon number of `transcript_start_codon_pos` | as `transcript_start_codon_pos` |
| `transcript_last_codon` | string | Last 3 nt of `alt_transcript_seq` | not scanned |
| `transcript_valid_stop` | bool | Whether `transcript_last_codon` is a stop codon | not scanned |
| `transcript_first_stop_codon` | string | First stop codon that the scan found | not scanned; no stop codon found |
| `transcript_first_stop_pos` | int | Transcript position of `transcript_first_stop_codon` | as `transcript_first_stop_codon` |
| `transcript_num_stop_codons` | int | Number of stop codons that the scan found | not scanned |
| `transcript_all_stop_codons` | stop_codon_list | (transcript position, codon) of each stop codon that the scan found | not scanned |
| `transcript_stop_codon_exons` | int_list | Exon number of each stop codon that the scan found | not scanned |
| `unknown_reason` | string | Why the alt transcript is unknown: `splice_site_destroyed` or `exon_boundary_ambiguous` (see [Variants at exon boundaries](#variants-at-exon-boundaries)). The section also lists the columns that are null on such a row | the alt transcript is known |

### NMD features

`add_nmd_features` adds these 11 columns. The [figures](#figures-of-the-features-and-the-nmd-rules) below draw them on an example transcript, and each case of `ptc_to_intron`.

| Column | Kind | Meaning | Null when |
|---|---|---|---|
| `utr3_length` | int | Length of the 3' UTR of the ref transcript: `transcript_length` minus `cds_end_in_transcript` | `has_stop_codon` is False; `cds_end_in_transcript` is null |
| `utr5_length` | int | Length of the 5' UTR of the ref transcript: `cds_start_in_transcript`. It ends at CDS position 0, not at the start codon | `cds_start_in_transcript` is null |
| `total_exon_count` | int | Number of exons in `transcript_exon_info` | the transcript has no exon rows |
| `upstream_exon_count` | int | Number of exons upstream of the PTC exon | not a PTC row; the PTC exon is not in `transcript_exon_info` |
| `downstream_exon_count` | int | Number of exons downstream of the PTC exon | as `upstream_exon_count` |
| `ptc_to_start_codon` | int | Distance in nt from the start codon to the PTC. Without a start loss: from the annotated start codon, `alt_first_stop_pos - alt_start_codon_pos`. After a start loss, translation starts at the ATG of the scan: `transcript_first_stop_pos - transcript_start_codon_pos`, both positions in `alt_transcript_seq` | not a PTC row; without a start loss: the transcript has no annotated start codon; after a start loss: the scan found no ATG or no stop codon |
| `ptc_less_than_150nt_to_start` | bool | Whether `ptc_to_start_codon` is less than 150. False if `ptc_to_start_codon` is null | `unknown_reason` is set |
| `ptc_exon_length` | int | Length of the PTC exon, UTR included | not a PTC row; the PTC exon is not in `transcript_exon_info` |
| `stop_codon_distance` | int | Distance in nt from the first in-frame stop codon of the alt transcript to the annotated stop codon, both at their positions in `alt_transcript_seq`. Unless the row keeps the flags from the CDS, it compares the same two stop codons as the [stop codon classification](#stop-codon-classification), so its sign gives the class. Positive: a PTC upstream of the annotated stop codon. On a PTC row, the first stop codon is the PTC. 0: the annotated stop codon, also after an insertion inside it (TAA>TGAA) or an in-frame indel right before it. Negative: a stop loss, with the new stop codon downstream, e.g. -3 for ATAG>T in GTA TAG TAG CAT. On a row that keeps the flags from the CDS, the first stop codon is the one at `alt_first_stop_pos`. Without `alt_transcript_seq`, both positions are in alt CDS coordinates, and the annotated stop codon starts at `alt_cds_len - 3` | `unknown_reason` is set; `has_stop_codon` is False; a nonstop: the alt transcript has no in-frame stop codon; on a row that keeps the flags from the CDS, the alt CDS has no in-frame stop codon |
| `ptc_to_intron` | int | Distance in nt from the PTC to the 3' end of the PTC exon. That end is the downstream exon junction, or the transcript end for the last exon | not a PTC row; the PTC exon is not in `transcript_exon_info`; `cds_start_in_transcript` is null |
| `likely_misannotated` | bool | True if `cds_in_transcript` is False, `ref_start_codon_pos` is not 0, or `ref_valid_stop` is False. A null in one of these 3 columns gives True too, e.g. `ref_start_codon_pos` of a transcript without an annotated start codon | never |

### NMD rules

`evaluate_nmd_escape_rules` adds these 6 columns. On a row with `unknown_reason`, every rule is null. On any other row that is not a PTC row, every rule is False. A rule whose input is null is False too. The [figures](#figures-of-the-features-and-the-nmd-rules) below define the rules in detail.

| Column | Kind | Meaning | Null when |
|---|---|---|---|
| `nmd_last_exon_rule` | bool | The PTC lies in the last exon: `downstream_exon_count` is 0 | `unknown_reason` is set |
| `nmd_50nt_penultimate_rule` | bool | The PTC lies 1 to 50 nt upstream of the last exon junction | `unknown_reason` is set |
| `nmd_long_exon_rule` | bool | The PTC exon has more than 407 nt: `ptc_exon_length` > 407 | `unknown_reason` is set |
| `nmd_start_proximal_rule` | bool | The PTC lies less than 150 nt downstream of the start codon, which after a start loss is the ATG of the scan. It equals `ptc_less_than_150nt_to_start` | `unknown_reason` is set |
| `nmd_single_exon_rule` | bool | The transcript has one exon: `total_exon_count` is 1 | `unknown_reason` is set |
| `nmd_escape` | bool | One of the 5 rules above is True | `unknown_reason` is set |

## Figures of the features and the NMD rules

The figures show a transcript 5' to 3' and are not to scale. `u` is UTR, `=` is CDS, `[...]` is an exon, `|` between two exons is an exon junction, and `*` is the PTC. `*--->|` is the distance from the PTC to an exon junction or to the transcript end, and `<--->` is a length. The numbers are positions in CDS coordinates, as alt_first_stop_pos: 0 is the 5' base of the CDS, here also of the start codon, and the 5' UTR has negative positions. After an indel, the positions are in alt CDS coordinates, so the exon junctions move with the PTC (see exon_end_in_alt_cds() in extra_features.py). The figures and the exon numbers follow the transcript. On the minus strand, the genomic coordinates run the other way.

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

nmd_start_proximal_rule: the PTC lies less than 150 nt downstream of the start codon, i.e. ptc_to_start_codon < 150. The start codon is the annotated one, and without one the rule is False. After a start loss, the start codon is the ATG that the scan found, and the distance runs to the first in-frame stop codon after it: transcript_first_stop_pos - transcript_start_codon_pos < 150.

```
5' [uuu=====*===]|[============]|[=========uuuuuuu] 3'
       <---->  < 150 nt
```
