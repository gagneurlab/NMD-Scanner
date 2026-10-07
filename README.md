# NMD variant effect prediction

The NMD-Scanner is a Python-based variant effect annotation tool that predicts the likelihood of transcript degradation through nonsense-mediated decay (NMD).
It reconstructs reference and alternative coding sequences as well as transcript sequences in some cases, identifies premature termination codons (PTCs), and evaluates canonical and non-canonical NMD escape rules.
It can handle single-nucleotide variants, multiple base substitutions, long and short deletions and duplications as well as frameshift variants.

## Features
- Reads variants from VCF and gene annotations from GFF3 (GENCODE or Ensembl flavor), gzip-compressed or not
- Reconstructs reference and alternative CDS, reference transcript sequence and (in some cases) the alternative transcript sequences with metadata
- Detects start / stop-loss and premature termination codons (PTCs) with the exact position in the CDS and in which exon it lies
- Computes different NMD-related features:
  - Total, upstream and downstream exon count
  - Distance of PTC to original stop codon (0 for the original stop codon, negative for the new stop codon after a stop loss, empty for a nonstop)
  - Distance of PTC to start codon
  - Transcript length
  - 3' and 5' UTR lengths
- Evaluates five canonical NMD escape rules:
  - Last exon rule
  - 50nt penultimate rule
  - Long exon rule
  - Start-proximal rule
  - Single-exon rule
- Writes the results as CSV or Parquet, or returns them as a pandas DataFrame, with the same 77 columns and dtypes for every input, or 73 without the sequence columns

[Technical Notes](Technical%20Notes.md) defines the features and the NMD escape rules, with figures.

## Installation
Requires Python >= 3.12.

From [PyPI](https://pypi.org/project/nmd-scanner/):
```bash
pip install nmd-scanner
```

`polars-bio` reads the VCF and GFF3 files and `pyarrow` writes Parquet output. Both come with the package, so Parquet output needs no extra install. `polars-bio` and its dependencies add about 700 MB to the install.

## Usage

### Option 1: Annotating a VCF on the command line

After the install, the `nmd-scanner` command is available:
```bash
nmd-scanner --vcf input.vcf --annotation annotation.gff3.gz --fasta reference.fa --output results/input.csv

# write Parquet instead of CSV
nmd-scanner --vcf input.vcf --annotation annotation.gff3.gz --fasta reference.fa --output results/input.parquet

# option: fix exon numbering (recommended for hg19)
nmd-scanner --vcf input.vcf --annotation annotation.gff3.gz --fasta reference.fa --output results/input.csv --reassign_exons

# option: leave out the 4 sequence columns
nmd-scanner --vcf input.vcf --annotation annotation.gff3.gz --fasta reference.fa --output results/input.parquet --no-sequences
```

The equivalent `python -m nmd_scanner.cli ...` invocation also works without installing the console script.

Arguments:
- `--vcf`: Path to input VCF, plain or gzip-compressed, with one ALT allele per record. Split multi-allelic records first, e.g. with `bcftools norm -m-`. Records with a symbolic, breakend, `.` or `*` ALT are skipped with a warning.
- `--annotation`: Path to the gene annotation in GFF3 (GENCODE or Ensembl flavor), optionally gzip-compressed, with the suffix `.gff3` or `.gff`.
- `--fasta`: Path to reference genome FASTA.
- `--output`: Path to the output file. Extension selects the format: `.csv` for CSV, `.parquet` or `.pq` for Parquet. The parent directory must already exist; the file is overwritten if present.
- `--reassign_exons`: (flag) Recompute exon numbers (recommended for hg19)
- `--no-sequences`: (flag) Leave out the 4 sequence columns `ref_cds_seq`, `alt_cds_seq`, `transcript_seq` and `alt_transcript_seq`. The output then has 73 columns. On the pathogenic ClinVar variants of chr22, these 4 columns make up 92% of the compressed Parquet bytes and about 216 of 272 MiB in pandas.

The chromosome names must match in the VCF, the GFF3 and the FASTA, e.g. all `chr1` or all `1`.

Output:
- The file specified by `--output`, containing:
  - reconstructed reference / alternative CDS and transcript sequences (+ metadata); `--no-sequences` leaves out the sequences
  - PTC detection and start / stop-loss flags
  - NMD escape rules
  - extra features such as UTR lengths, exon counts, distances, etc.
  - `unknown_reason`: empty if the alt transcript is known. Otherwise it says why the alt transcript is unknown: `splice_site_destroyed` or `exon_boundary_ambiguous`. The alt columns, the start / stop-loss flags, the PTC features and the NMD escape rules are empty then.

`nmd_scanner.schema` lists the 77 output columns and their dtypes. `nmd_scanner.schema.output_column_kinds(sequences=False)` gives the 73 columns without the sequences. The columns and dtypes are the same for every input, also for a result without rows. [Output columns](Technical%20Notes.md#output-columns) in the Technical Notes gives the meaning of each column and says when it is null.

### Option 2: Import as a python module
Instead of running the entire pipeline, you can import NMD-Scanner in Python and call only specific components.
This is useful if you want to 
- only reconstruct transcript / CDS sequences
- only compute NMD escape rules
- integrate NMD-Scanner into a larger workflow
- build custom features

To get the result table of the CLI as a `pandas.DataFrame` without writing it, call `annotate`. It takes the inputs and options of the CLI, except the output path; `sequences=False` is the flag `--no-sequences`. It does not configure logging or write files:
```python
import nmd_scanner

results = nmd_scanner.annotate("input.vcf", "annotation.gff3.gz", "reference.fa", reassign_exons=False)
results["my_key"] = "sample_1"  # add your own columns
results.to_csv("results.csv", index=False)
```
`DataFrame.to_csv` writes a list column, e.g. `transcript_exons`, as its Python repr. The CLI writes it as JSON instead (see [Kinds and dtypes](Technical%20Notes.md#kinds-and-dtypes)).

To convert the table to a `pyarrow.Table`, call `to_arrow`. Each column gets the Arrow type of its kind, so the types are the same for every input, also for a result without rows. The CLI writes this table for Parquet output. `to_arrow` takes the output columns, with or without the sequences; a column of your own, such as `my_key` above, raises a KeyError:
```python
import pyarrow.parquet as pq

results = nmd_scanner.annotate("input.vcf", "annotation.gff3.gz", "reference.fa")
pq.write_table(nmd_scanner.to_arrow(results), "results.parquet")
```

`polars-bio` shows tqdm progress bars on stderr, e.g. one for every file it reads. To turn them off, set `TQDM_DISABLE=1` in the environment before Python starts. Importing `nmd_scanner` imports `polars-bio`, which sets `POLARS_FORCE_NEW_STREAMING` in `os.environ` if it is not set, and adds 4 filters to the `warnings` module.

For reconstructing reference and alternative coding and transcript sequences, PTC detection and start / stop-loss information:
```python
import pandas as pd
from pyfaidx import Fasta

import nmd_scanner

vcf = nmd_scanner.read_vcf("input.vcf")
fasta = Fasta("reference.fa")
# exon rows and coding regions: CDS rows that include the stop codon, with the columns has_start_codon and
# has_stop_codon.
# Optional: reassign_exons=True recomputes the exon numbers (recommended for hg19).
annotation = nmd_scanner.read_annotation("annotation.gff3.gz", fasta, reassign_exons=False)

cds_df = annotation[annotation["Feature"] == "CDS"]
exons_df = annotation[annotation["Feature"] == "exon"].copy()
exons_df["exon_length"] = exons_df["End"] - exons_df["Start"]

results = nmd_scanner.extract_ptc(cds_df, vcf, fasta, exons_df)
```

`results` has one row per VCF record and transcript where the variant of the record touches the coding region (CDS plus stop codon) or the splice dinucleotide next to one of its exon edges. An indel in a repeat counts with every equivalent placement, not only the one in the VCF. If the variant destroys a splice site or leaves an exon boundary ambiguous, the alt transcript is unknown. Then `unknown_reason` says why, and the alt columns and the predictions are empty. [Variants at exon boundaries](Technical%20Notes.md#variants-at-exon-boundaries) in the Technical Notes gives the rules.

Add the extra NMD-related features (utr lengths, exon counts, ptc-related features) and the NMD escape rules (last exon rule, 50 nt penultimate rule, long exon rule, start proximal rule, single exon rule, nmd escape) to the above computed results. The columns, column order and dtypes are the same for every input, also for a result without rows (see `nmd_scanner.schema`):
```python
results = nmd_scanner.add_features_and_rules(results)
```

The last column, `nmd_model_status`, says whether the NMD efficiency model can score a row: `ok`, or the first reason why not, e.g. `no_annotated_stop`. `nmd_scanner.schema.MODEL_INPUTS` lists the 19 model inputs in model order. [Model status](Technical%20Notes.md#model-status) in the Technical Notes lists all values:
```python
from nmd_scanner.schema import MODEL_INPUTS

scorable = results[results["nmd_model_status"] == "ok"]
X = scorable[MODEL_INPUTS]  # no null value
```

To work on single rows, `nmd_scanner.add_nmd_features` and `nmd_scanner.evaluate_nmd_escape_rules` are public too. Run the features **before** the escape rules: `evaluate_nmd_escape_rules` reads the exon-count and ptc-exon-length columns produced by the features. The row functions do not add `nmd_model_status`; `add_features_and_rules` adds it.

## License
All source code in this repository is licensed under the [MIT License](https://github.com/gagneurlab/NMD-Scanner/blob/main/LICENSE).

## Citation 
Schröder, C.H. (2025). *Enhanced Aberrant Gene Expression Prediction across Human Tissues*.
Master's Thesis, Technical University of Munich / Ludwig-Maximilians-Universität München.
