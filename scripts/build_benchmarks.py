"""
Build the NMD efficiency benchmarks of nmd_scanner.benchmark from their pinned public sources, and write each one to
<out-dir>/<name>.parquet.

Run it with:

    uv run --extra scripts python scripts/build_benchmarks.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa \
        --out-dir benchmarks/

--dataset builds only the named benchmarks, e.g. --dataset gtex --dataset geuvadis. By default, it builds all.
"""

import argparse
import logging
from pathlib import Path

from pyfaidx import Fasta

from nmd_scanner.benchmark import DATASETS
from nmd_scanner.benchmark.sources import read_transcripts

logger = logging.getLogger("build_benchmarks")

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--gff3", required=True, type=Path, help="GENCODE GFF3 of GRCh38, optionally gzip-compressed")
parser.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA with chr names and a .fai index")
parser.add_argument("--out-dir", required=True, type=Path, help="directory for the benchmark tables")
parser.add_argument(
    "--downloads",
    type=Path,
    help="directory for the downloaded study tables, by default <out-dir>/downloads. A study table that is there"
    " already is read from there, with the same sha256 check.",
)
parser.add_argument(
    "--dataset", action="append", choices=list(DATASETS), help="a benchmark to build; repeat it for more"
)
args = parser.parse_args()
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
args.out_dir.mkdir(parents=True, exist_ok=True)
downloads = args.downloads or args.out_dir / "downloads"
transcripts = read_transcripts(args.gff3)
fasta = Fasta(str(args.fasta))
for name in args.dataset or DATASETS:
    dataset = DATASETS[name]
    table = dataset.build(dataset.download(downloads), transcripts, fasta)
    path = args.out_dir / f"{name}.parquet"
    table.write_parquet(path)
    logger.info("Wrote %d rows to %s", len(table), path)
