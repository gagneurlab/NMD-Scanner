# Input Defects

tl;dr: "Technical Notes.md" assumes a well-formed GFF3 and matching input files. This file lists the input that
NMD-Scanner rejects, the misannotations that it accepts, and the null cases that they add to the output columns.

## Rejected input

- `scan.read_gff3` raises a ValueError for an exon or CDS row whose strand is not + or -, for a CDS row whose phase
  is not 0, 1 or 2, and for a CDS row that does not lie inside an exon row of its transcript. The error names the
  transcript of the first such row. A transcript without exon rows is not checked.
- `rules.extract_ptc` raises a ValueError if the CDS rows lack the column `has_start_codon` or `has_stop_codon`,
  i.e. are no coding regions. It also raises one if the FASTA has no sequence for a chromosome with variants and CDS
  rows, and the error names these chromosomes.

## Misannotations

NMD-Scanner reads a transcript with one of these misannotations. The null cases below say which columns turn null.

- An annotated start codon that is a stop codon, such as TAG. Translation cannot start on a stop codon.
  `likely_misannotated` does not flag it, because its start codon check only asks for an annotated start codon at
  CDS position 0. A variant that leaves it unchanged gives a PTC row whose PTC is this start codon.
- An annotated stop codon that is out of frame, because the CDS length does not fit `cds_frame`, because it is no
  stop codon, or because an in-frame stop codon lies upstream of it. The rows of the transcript keep the flags from
  the CDS ("Technical Notes.md", section "Stop codon classification").
- A transcript without exon rows. Its rows have no transcript sequence, and `cds_in_transcript` is False.
- An exon row that the GFF3 holds twice. `transcript_seq` then holds the exon twice. After an indel in that exon,
  the exon lengths do not add up to `alt_transcript_length`.

## Null cases

These clauses add to the "Null when" column of the tables in "Technical Notes.md". "as `x`" stands for the clauses
of column `x` in both files.

| Column | Null when |
|---|---|
| `transcript_start` | the transcript has no exon rows |
| `transcript_end` | as `transcript_start` |
| `transcript_seq` | as `transcript_start` |
| `transcript_length` | as `transcript_start` |
| `cds_start_in_transcript` | the transcript has no exon rows; the 5' CDS base lies outside the exons |
| `cds_end_in_transcript` | as `cds_start_in_transcript` |
| `alt_transcript_seq` | `cds_start_in_transcript` is null |
| `transcript_exons` | the transcript has no exon rows |
| `alt_transcript_exons` | the exon lengths do not add up to `alt_transcript_length` |
| `total_exon_count` | the transcript has no exon rows |
| `ptc_to_start_codon` | the annotated start codon is a stop codon, such as TAG |
