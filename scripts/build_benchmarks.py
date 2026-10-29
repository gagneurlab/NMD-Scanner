"""
Build the NMD efficiency benchmarks of nmd_scanner.benchmark from their pinned public sources, and write each one to
<out-dir>/<name>.parquet.

Run it with:

\b
    uv run --extra scripts python scripts/build_benchmarks.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa \
        --out-dir benchmarks/

--dataset builds only the named benchmarks, e.g. --dataset gtex --dataset geuvadis. By default, it builds all.
"""

import logging
import sys
from pathlib import Path
from types import SimpleNamespace

import click
from pyfaidx import Fasta

from nmd_scanner.benchmark import DATASETS
from nmd_scanner.benchmark.sources import read_transcripts

logger = logging.getLogger("build_benchmarks")


@click.command(help=__doc__)
@click.option(
    "--gff3",
    required=True,
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    help="GENCODE GFF3 of GRCh38, optionally gzip-compressed",
)
@click.option(
    "--fasta",
    required=True,
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    help="GRCh38 FASTA with chr names and a .fai index",
)
@click.option(
    "--out-dir",
    required=True,
    type=click.Path(path_type=Path, file_okay=False),
    help="directory for the benchmark tables",
)
@click.option(
    "--downloads",
    type=click.Path(path_type=Path, file_okay=False),
    help="directory for the downloaded study tables, by default <out-dir>/downloads. A study table that is there"
    " already is read from there, with the same sha256 check.",
)
@click.option(
    "--dataset", multiple=True, type=click.Choice(list(DATASETS)), help="a benchmark to build; repeat it for more"
)
def command(**params: object) -> None:
    pass


try:
    args = SimpleNamespace(**command.make_context("build_benchmarks.py", sys.argv[1:]).params)
except click.exceptions.Exit as exit_:
    sys.exit(exit_.exit_code)
except click.ClickException as error:
    error.show()
    sys.exit(error.exit_code)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
args.out_dir.mkdir(parents=True, exist_ok=True)
downloads = args.downloads or args.out_dir / "downloads"
transcripts_df = read_transcripts(args.gff3)
fasta = Fasta(str(args.fasta))
for name in args.dataset or DATASETS:
    dataset = DATASETS[name]
    table_df = dataset.build(dataset.download(downloads), transcripts_df, fasta)
    path = args.out_dir / f"{name}.parquet"
    table_df.write_parquet(path)
    logger.info("Wrote %d rows to %s", len(table_df), path)
