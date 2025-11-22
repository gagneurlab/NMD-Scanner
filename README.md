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
pip install .
```

## Usage

### Option 1: Annotating a VCF on the command line
```bash
# if running the script directly
python -m nmd_scanner.cli --vcf input.vcf --gtf annotation.gtf --fasta reference.fa --output results/

# option: fix exon numbering (recommended for hg19)
python -m nmd_scanner.cli --vcf input.vcf --gtf annotation.gtf --fasta reference.fa --output results/ --reassign_exons
```

Arguments:
- `--vcf`: Path to input VCF (SNVs / Indels supported; frameshifts handled)
- `--gtf`: Path to gene annotation (GTF)
- `--fasta`: Path to reference genome FASTA
- `--output`: Path to an existing directory (or a file path whose parent exists)
- `--reassign_exons`: (flag) Recompute exon numbers (useful for hg19)

Output:
- A CSV named `<vcf_basename>_final_nmd_results.csv` saved to `--output`, containing:
  - reconstructed reference / alternative CDS and transcript sequences(+ metadata)
  - PTC detection and start / stop-loss flags
  - NMD escape rules
  - extra features such as UTR lengths, exon counts, distances, etc.)

### Option 2: Import as a python moduele
Instead of running the entire pipeline, you can import NMD-Scanner in Python and call only specific components.
This is useful if you want to 
- only reconstruct transcript / CDS sequences
- only compute NMD escape rules
- integrate NMD-Scanner into a larger workflow
- build custom features

For reconstructing reference and alternative coding and transcript sequences, PTC detection and start / stop-loss information:
```python
import pandas as pd
import pyranges as pr
from pyfaidx import Fasta

import nmd_scanner

vcf = nmd_scanner.read_vcf("input.vcf")
gtf_pr = nmd_scanner.read_gtf("annotation.gtf")
fasta = Fasta("reference.fa")

# Optional: fix exon numbering (recommended for hg19)
gtf_pr = nmd_scanner.compute_exon_numbers(gtf_pr)

gtf_df = gtf_pr.df
cds_df = gtf_df[gtf_df["Feature"] == "CDS"]
exons_df = gtf_df[gtf_df["Feature"] == "exon"].copy()
exons_df["exon_length"] = exons_df["End"] - exons_df["Start"]

results = extract_ptc(cds_df, vcf, fasta, exons_df, output="tmp/")
```

Add NMD escape rules (last exon rule, 50 nt penultimate rule, long exon rule, start proximal rule, single exon rule, nmd escape) to the above computed results:
```python
nmd_results = results.apply(nmd_scanner.evaluate_nmd_escape_rules, axis=1, result_type='expand')
results = pd.concat([results, nmd_results], axis=1)
```

Add extra NMD-related features (utr lengths, exon counts, ptc-related features) to above computed results:
```python
extra_features = results.apply(nmd_scanner.add_nmd_features, axis=1, result_type='expand')
results = pd.concat([results, extra_features], axis=1)
```

## License
All source code in this repository is licensed under the [MIT License](./License).

## Citation 
Schröder, C.H. (2025). *Enhanced Aberrant Gene Expression Prediction across Human Tissues*.
Master's Thesis, Technical University of Munich / Ludwig-Maximilians-Universität München.





[//]: # (## Description)

[//]: # (Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.)

[//]: # ()
[//]: # (## Installation)

[//]: # (Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.)

[//]: # ()
[//]: # (## Usage)

[//]: # (Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.)

[//]: # ()
[//]: # (## Support)

[//]: # (Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.)

[//]: # ()
[//]: # (## Contributing)

[//]: # (State if you are open to contributions and what your requirements are for accepting them.)

[//]: # ()
[//]: # (For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.)

[//]: # ()
[//]: # (You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.)

[//]: # ()
[//]: # (## Authors and acknowledgment)

[//]: # (Show your appreciation to those who have contributed to the project.)

[//]: # ()
[//]: # (## License)

[//]: # (For open source projects, say how it is licensed.)

[//]: # ()
[//]: # (## Project status)

[//]: # (If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.)
