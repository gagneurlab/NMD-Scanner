# Input Defects

tl;dr: "Technical Notes.md" assumes a well-formed GFF3 and matching input files. This file lists the input that
NMD-Scanner rejects, the misannotations that it accepts, and the null cases that they add to the output columns.

## Rejected input

- `scan.read_gff3` raises a ValueError for an exon or CDS record whose strand is not + or -, for a CDS record whose
  phase is not 0, 1 or 2, and for a CDS record that does not lie inside an exon record of its transcript. The error
  names the transcript of the first such record. A transcript without exon records is not checked.
- `scan.read_gff3` raises a ValueError for an Ensembl transcript without an mRNA record: "Could not resolve
  transcript_id via the ID/Parent hierarchy for N record(s)." Its exon and CDS records get their `transcript_id` from
  that record only.
- `rules.extract_ptc` raises a ValueError if the CDS records lack the column `has_start_codon`, `has_stop_codon`,
  `transcript_type` or `transcript_tags`, i.e. are not the coding regions that `scan.read_annotation` returns. It also
  raises one if the FASTA has no sequence for a chromosome with variants and CDS records, and the error names these
  chromosomes.

## Misannotations

NMD-Scanner reads a transcript with one of these misannotations. The null cases below say which columns turn null.

- An annotated start codon that is a stop codon, such as TAG. Translation cannot start on a stop codon.
  `likely_misannotated` does not flag it, because its start codon check only asks for an annotated start codon at CDS
  position 0. A variant that leaves it unchanged gives a PTC variant whose PTC is this start codon.
- An annotated stop codon that is out of frame, because the CDS length does not fit `cds_frame`, because it is no stop
  codon, or because an in-frame stop codon lies upstream of it. The variant-transcript pairs of the transcript keep the
  flags from the CDS ("Technical Notes.md", section "Stop codon classification").
- A transcript without exon records. Its variant-transcript pairs have no transcript sequence, and `cds_in_transcript`
  is False.
- A GENCODE transcript without a transcript record, or a transcript record without the type of the transcript:
  `transcript_type` in GENCODE, `biotype` in Ensembl. `transcript_type` and `transcript_tags` come from the transcript
  record only, and the exon and CDS records do not fill them in.
- An exon record that the GFF3 holds twice. `transcript_seq` then holds the exon twice. After an indel in that exon, the
  exon lengths do not add up to `alt_transcript_length`.

## Null cases

These clauses add to the "Null when" column of the tables in "Technical Notes.md". "as `x`" stands for the clauses
of column `x` in both files.

| Column                    | Null when                                                                                                                        |
| ------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `transcript_type`         | the GENCODE transcript has no transcript record; the transcript record has no `transcript_type` (GENCODE) or `biotype` (Ensembl) |
| `transcript_tags`         | the GENCODE transcript has no transcript record                                                                                  |
| `transcript_start`        | the transcript has no exon records                                                                                               |
| `transcript_end`          | as `transcript_start`                                                                                                            |
| `transcript_seq`          | as `transcript_start`                                                                                                            |
| `transcript_length`       | as `transcript_start`                                                                                                            |
| `cds_start_in_transcript` | the transcript has no exon records; the 5' CDS base lies outside the exons                                                       |
| `cds_end_in_transcript`   | as `cds_start_in_transcript`                                                                                                     |
| `alt_transcript_seq`      | `cds_start_in_transcript` is null                                                                                                |
| `transcript_exons`        | the transcript has no exon records                                                                                               |
| `alt_transcript_exons`    | the exon lengths do not add up to `alt_transcript_length`                                                                        |
| `total_exon_count`        | the transcript has no exon records                                                                                               |
| `ptc_to_start_codon`      | the annotated start codon is a stop codon, such as TAG                                                                           |
