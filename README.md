# NMD variant effect prediction

The NMD-Scanner is a Python-based variant effect annotation tool that predicts the likelihood of transcript degradation through nonsense-mediated decay (NMD).
It reconstructs reference and alternative coding sequences as well as transcript sequences in some cases, identifies premature termination codons (PTCs), and evaluates canonical and non-canonical NMD escape rules.
It can handle single-nucleotide variants, multiple base substitutions, long and short deletions and duplications as well as frameshift variants.

## Features
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
- Outputs all annotations as a structured DataFrame (CSV)

## Installation
```bash
git clone https://github.com/gagneurlab/NMD-Scanner.git
cd NMD-Scanner
pip install -e .
```

## Usage

### Option 1: Using the command-line tool
After installation, you can use the `nmd_scanner` command directly:

```bash
# Basic usage (outputs CSV by default)
nmd_scanner --vcf input.vcf --gtf annotation.gtf --fasta reference.fa --output results/

# Output as Parquet
nmd_scanner --vcf input.vcf --gtf annotation.gtf --fasta reference.fa --output results/ --output-format parquet

# With canonical transcripts only
nmd_scanner --vcf input.vcf --gtf annotation.gtf --fasta reference.fa --output results/ --canonical-only

# With exon numbering fix (recommended for hg19)
nmd_scanner --vcf input.vcf --gtf annotation.gtf --fasta reference.fa --output results/ --reassign-exons
```

Alternatively, you can run it as a Python module:
```bash
python -m nmd_scanner.cli --vcf input.vcf --gtf annotation.gtf --fasta reference.fa --output results/
```

Arguments:
- `--vcf`: Path to input VCF (SNVs / Indels supported; frameshifts handled)
- `--gtf`: Path to gene annotation (GTF)
- `--fasta`: Path to reference genome FASTA
- `--output`: Path to an existing directory or file (supports .csv and .parquet extensions)
- `--output-format`: Output format: 'csv' or 'parquet' (default: inferred from file extension or 'csv')
- `--canonical-only`: (flag) Only process canonical transcripts
- `--reassign-exons`: (flag) Recompute exon numbers (useful for hg19)
- `-v, --verbose`: Increase verbosity (use -vv for debug output)

Output:
- A CSV or Parquet file named `<vcf_basename>_nmdscanner.csv` (or `.parquet`) saved to `--output`, containing:
  - Variant information (ID, chromosome, position, ref, alt)
  - Transcript and gene IDs
  - PTC detection (`alt_is_premature`)
  - NMD escape rules (last exon, 50nt penultimate, long exon, start proximal, single exon)
  - NMD efficiency prediction
  - Transcript features (UTR lengths, exon counts, distances, etc.)

Note: Only variants with premature stop codons (`alt_is_premature=True`) are included in the output file.

### Option 2: Import as a Python module

#### Get DataFrame Results

The simplest way to use NMD-Scanner as a library:

```python
import nmd_scanner

# Returns a pandas DataFrame with all annotations
df = nmd_scanner.annotate_nmd_pandas(
    vcf_path="input.vcf",
    gtf_path="annotation.gtf",
    fasta_path="reference.fa",
    output="results/",  # Optional: save to file
    canonical_only=True,  # Optional: only canonical transcripts
    reassign_exons=False  # Optional: fix exon numbering (for hg19)
)

# Access results
print(f"Total variants analyzed: {len(df)}")
ptc_variants = df[df['alt_is_premature']]
print(f"Variants with PTCs: {len(ptc_variants)}")
print(f"NMD escape rate: {ptc_variants['nmd_escape'].mean():.2%}")
```

#### Get Structured Objects

For programmatic access with type safety:

```python
import nmd_scanner

# Returns a list of NMDResult objects
results = nmd_scanner.annotate_nmd(
    vcf_path="input.vcf",
    gtf_path="annotation.gtf",
    fasta_path="reference.fa"
)

# Access individual results
for result in results:
    if result.cds_annotation.alt_is_premature:
        print(f"Variant: {result.cds_annotation.variant_id}")
        print(f"Transcript: {result.cds_annotation.transcript_id}")
        print(f"NMD escape: {result.nmd_prediction.nmd_escape}")
        print(f"NMD efficiency: {result.nmd_prediction.nmd_efficiency}")
        print(f"Downstream exons: {result.transcript_features.downstream_exon_count}")

# Convert to DataFrame if needed
df = nmd_scanner.NMDResult.to_dataframe(results)
```

## License
All source code in this repository is licensed under the [MIT License](./License).

## Citation 
Schröder, C.H. (2025). *Enhanced Aberrant Gene Expression Prediction across Human Tissues*.
Master's Thesis, Technical University of Munich / Ludwig-Maximilians-Universität München.