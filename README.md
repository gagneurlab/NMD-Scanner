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
  - Distance of PTC to original stop codon
  - Distance of PTC to start codon
  - Transcript length
  - 3' and 5' UTR lengths
- Evaluates five canonical NMD escape rules:
  - Last exon rule
  - 50nt penultimate rule
  - Long exon rule
  - Start-proximal rule
  - Single-exon rule
- Writes the results as CSV or Parquet, or returns them as a pandas DataFrame, with the same 78 columns and dtypes for every input

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
```

The equivalent `python -m nmd_scanner.cli ...` invocation also works without installing the console script.

Arguments:
- `--vcf`: Path to input VCF, plain or gzip-compressed (SNVs / Indels supported; frameshifts handled). It needs its header, at least the `##fileformat` and `#CHROM` lines, and one ALT allele per record: split multi-allelic records first, e.g. with `bcftools norm -m-`. QUAL, FILTER and INFO are not read. A record whose ALT is a symbolic allele (e.g. `<DEL>`, `<DUP>`, `<INS>`, `<INV>`, `<CNV>`) or a breakend is skipped, and a warning gives their count: NMD-Scanner cannot apply structural variants yet.
- `--annotation`: Path to gene annotation file in GFF3, optionally gzip-compressed, with the suffix `.gff3` or `.gff`. Both GENCODE and Ensembl GFF3 flavors are supported.
- `--fasta`: Path to reference genome FASTA. It also shows whether a CDS ends in a stop codon.
- `--output`: Path to the output file. Extension selects the format: `.csv` for CSV, `.parquet` or `.pq` for Parquet. The parent directory must already exist; the file is overwritten if present.
- `--reassign_exons`: (flag) Recompute exon numbers (recommended for hg19)

The chromosome names must match in the VCF, the GFF3 and the FASTA, e.g. all `chr1` or all `1`.

The coding region of a transcript is its CDS plus the stop codon. A GFF3 CDS includes the stop codon. Ensembl
GFF3 has no `stop_codon` rows; whether a transcript ends in a stop codon comes from the last 3 CDS bases in the
FASTA. Ensembl GFF3 has no `cds_end_NF` tag either. So a `cds_end_NF` transcript gets a stop codon if its CDS
ends in stop codon bases (13 transcripts in Ensembl 108, none on chr22).

Output:
- The file specified by `--output`, containing:
  - reconstructed reference / alternative CDS and transcript sequences (+ metadata)
  - PTC detection and start / stop-loss flags
  - NMD escape rules
  - extra features such as UTR lengths, exon counts, distances, etc.

`nmd_scanner.schema` lists the 78 output columns and their dtypes. The columns and dtypes are the same for every input, also for a result without rows.

### Option 2: Import as a python module
Instead of running the entire pipeline, you can import NMD-Scanner in Python and call only specific components.
This is useful if you want to 
- only reconstruct transcript / CDS sequences
- only compute NMD escape rules
- integrate NMD-Scanner into a larger workflow
- build custom features

To get the result table of the CLI as a `pandas.DataFrame` without writing it, call `annotate`. It takes the inputs and options of the CLI, except the output path. It does not configure logging or write files:
```python
import nmd_scanner

results = nmd_scanner.annotate("input.vcf", "annotation.gff3.gz", "reference.fa", reassign_exons=False)
results["my_key"] = "sample_1"  # add your own columns
results.to_csv("results.csv", index=False)
```

`polars-bio` shows tqdm progress bars on stderr, e.g. one for every file it reads. To turn them off, set `TQDM_DISABLE=1` in the environment before Python starts. Importing `nmd_scanner` imports `polars-bio`, which sets `POLARS_FORCE_NEW_STREAMING` in `os.environ` if it is not set, and adds 4 filters to the `warnings` module.

For reconstructing reference and alternative coding and transcript sequences, PTC detection and start / stop-loss information:
```python
import pandas as pd
from pyfaidx import Fasta

import nmd_scanner

vcf = nmd_scanner.read_vcf("input.vcf")
fasta = Fasta("reference.fa")
# exon rows and coding regions: CDS rows that include the stop codon, with the column has_stop_codon.
# The FASTA shows whether a CDS ends in a stop codon.
# Optional: reassign_exons=True recomputes the exon numbers (recommended for hg19).
annotation = nmd_scanner.read_annotation("annotation.gff3.gz", fasta, reassign_exons=False)

cds_df = annotation[annotation["Feature"] == "CDS"]
exons_df = annotation[annotation["Feature"] == "exon"].copy()
exons_df["exon_length"] = exons_df["End"] - exons_df["Start"]

results = nmd_scanner.extract_ptc(cds_df, vcf, fasta, exons_df)
```

Add the extra NMD-related features (utr lengths, exon counts, ptc-related features) and the NMD escape rules (last exon rule, 50 nt penultimate rule, long exon rule, start proximal rule, single exon rule, nmd escape) to the above computed results. The columns, column order and dtypes are the same for every input, also for a result without rows (see `nmd_scanner.schema`):
```python
results = nmd_scanner.add_features_and_rules(results)
```

To work on single rows, `nmd_scanner.add_nmd_features` and `nmd_scanner.evaluate_nmd_escape_rules` are public too. Run the features **before** the escape rules: `evaluate_nmd_escape_rules` reads the exon-count and ptc-exon-length columns produced by the features.

## License
All source code in this repository is licensed under the [MIT License](https://github.com/gagneurlab/NMD-Scanner/blob/main/LICENSE).

## Citation 
Schröder, C.H. (2025). *Enhanced Aberrant Gene Expression Prediction across Human Tissues*.
Master's Thesis, Technical University of Munich / Ludwig-Maximilians-Universität München.
