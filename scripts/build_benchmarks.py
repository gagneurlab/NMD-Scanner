"""
Build the NMD efficiency benchmarks of nmd_scanner.benchmark from their pinned public sources, and write each one to
<out-dir>/<name>.parquet.

Run it from the repo root with:

\b
    uv run --extra scripts python scripts/build_benchmarks.py

--dataset builds only the named benchmarks, e.g. --dataset gtex --dataset geuvadis. By default, it builds all.
"""

import logging
from pathlib import Path

import click
import polars as pl
from pyfaidx import Fasta

from nmd_scanner import read_annotation
from nmd_scanner.benchmark import DATASETS

logger = logging.getLogger("build_benchmarks")


@click.command(help=__doc__)
@click.option(
    "--gff3",
    default="data/gencode.v42.annotation.gff3.gz",
    show_default=True,
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    help="GENCODE GFF3 of GRCh38, optionally gzip-compressed",
)
@click.option(
    "--fasta",
    default="data/GRCh38.fa",
    show_default=True,
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    help="GRCh38 FASTA with chr names and a .fai index",
)
@click.option(
    "--out-dir",
    default="data/benchmarks",
    show_default=True,
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
    "--dataset",
    "datasets",
    multiple=True,
    type=click.Choice(list(DATASETS)),
    help="a benchmark to build; repeat it for more",
)
def main(gff3: Path, fasta: Path, out_dir: Path, downloads: Path | None, datasets: tuple[str, ...]) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    out_dir.mkdir(parents=True, exist_ok=True)
    downloads = downloads or out_dir / "downloads"
    genome = Fasta(str(fasta))
    transcripts_df = (
        pl.from_pandas(read_annotation(gff3, genome)[["transcript_id", "Strand", "gene_id"]])
        .select(
            transcript=pl.col("transcript_id").str.split(".").list.first(),
            strand="Strand",
            gene=pl.col("gene_id").str.split(".").list.first(),
        )
        .unique()
    )
    if transcripts_df.is_empty():
        raise ValueError(f"{gff3} has no exon or CDS row")
    for name in datasets or DATASETS:
        dataset = DATASETS[name]
        table_df = dataset.build(dataset.download(downloads), transcripts_df, genome)
        path = out_dir / f"{name}.parquet"
        table_df.write_parquet(path)
        logger.info("Wrote %d rows to %s", len(table_df), path)


if __name__ == "__main__":
    main()
