# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
"""
Compute the numbers of the section "Model validation" of the main README: how well nmd_efficiency_rf.onnx and
nmd_escape alone rank the measured NMD efficiency in 4 cohorts.

Run it as a script with:

    uv run --group scripts python scripts/validate_model.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa \
        --train-dir train/ --out-dir validate/

--train-dir is the --out-dir of the run of scripts/train_model.py that trained nmd_efficiency_rf.onnx.

It is also a notebook (see scripts/README.md). In a notebook, the arguments come from the environment variable
VALIDATE_MODEL_ARGS.

Cohorts:
- TCGA: the out-of-fold predictions in oof_predictions.parquet of --train-dir.
- Geuvadis: germline nonsense variants in lymphoblastoid cell lines (Iha et al. 2025), Zenodo record 16666299. Their
  NMD efficiency is -log2(VAF_mean / 0.5). Only the variants of the study's accurate annotation count: no rescued MNV,
  no nontranslating gene, no mixed NMD isoforms.
- GTEx: germline variants in the tissues of GTEx v8, Supplementary Data 2 of Kim et al. 2024. The table has one row
  per variant and tissue (variant-tissue rows). The validation uses the mean of its column NMD_efficiency for each
  variant over its tissues. In 4,102 of the 17,816 variant-tissue rows with ASE_altered above 0,
  NMD_efficiency differs from -log2(ASE_altered / 0.5). The difference moves the Spearman correlations by only
  about 0.001.
- MMRF-TARGET: somatic variants in tumors, MMRF_TARGET_dataset.csv of NMDEff (https://github.com/hjkng/nmdeff).

Steps:
1. Download the tables of Geuvadis, GTEx and MMRF-TARGET and check their sha256. Zenodo has the Geuvadis tables only
   inside one zip of 3.2 GB, so fetch just the 2 members by HTTP range requests. With --geuvadis-dir, read the 2
   members from a local directory instead, with the same check.
2. GTEx and MMRF-TARGET give the alleles in HGVSc on the transcript strand. Convert them to the forward strand with
   add_alleles of make_benchmark_vcfs.py. Geuvadis gives genomic alleles. Fail if a REF does not match the FASTA.
3. Run nmd_scanner.annotate on the unique variants of the 3 cohorts. It adds the prediction of the model as
   nmd_efficiency_prediction. Join each measured row to the NMD-Scanner row of its transcript and variant.
4. Keep the rows with nmd_efficiency_prediction_status "ok" and a variant allele fraction above 0 in the RNA. Log the
   rows that this drops, and why.
5. Per cohort, compute the Spearman correlation of the predictions and of (not nmd_escape) with the measured NMD
   efficiency. The 95% CIs are percentiles of N_RESAMPLES bootstrap resamples of the genes. A resample draws genes
   with replacement and takes all rows of each drawn gene. Rows of one gene share a transcript structure, so they are
   not independent. The gene of a row is the gene of its transcript in the GFF3. The gain, model minus rules, has a
   paired CI from the same resamples. Each cohort draws from its own generator, seeded with SEED.
6. Per cohort, compute the mean measured NMD efficiency of the rows with nmd_escape False. For GTEx, also compute it
   over the variant-tissue rows, in total and per tissue with at least MIN_TISSUE_ROWS of them.
7. For MMRF-TARGET, compute the numbers of step 5 also for the rows on transcripts of the TCGA training rows, and for
   the other rows.

Outputs in --out-dir: variants.vcf (the input of nmd_scanner.annotate), validation.json, and figures/ with each
figure as PNG and PDF.
"""

import argparse
import io
import json
import logging
import os
import shlex
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import polars as pl
from figures import (
    AXIS,
    BLUE,
    ORANGE,
    SECONDARY_INK,
    calibration,
    calibration_plot,
    in_notebook,
    save_figure,
)
from make_benchmark_vcfs import (
    KEY,
    add_alleles,
    add_key,
    check_sha256,
    download,
    join_features,
    read_nmdeff,
    read_transcripts,
    write_vcf,
)
from plotnine import (
    aes,
    annotate,
    coord_flip,
    element_blank,
    facet_wrap,
    geom_abline,
    geom_boxplot,
    geom_hline,
    geom_point,
    geom_segment,
    geom_text,
    geom_vline,
    ggplot,
    labs,
    scale_color_identity,
    scale_y_continuous,
    stat_summary,
    theme,
)
from pyfaidx import Fasta
from scipy.stats import rankdata

import nmd_scanner
from nmd_scanner.model import MODEL_PATH

logger = logging.getLogger("validate_model")

# %%
# Supplementary Data 1 to 12 of Kim et al. 2024 (Commun Biol 7:1461), one sheet each
GTEX_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-024-07136-y/MediaObjects/"
    "42003_2024_7136_MOESM3_ESM.xlsx"
)
GTEX_SHA256 = "637c92b17095ef6147b74ae0ab6549f674adc6a67ee04c06d483fa7331ff34c6"
GTEX_SHEET = "Supplementary Data 2"
# Zenodo record 16666299 of Iha et al. 2025 has one file, NMD-rules.zip, of ZENODO_SIZE bytes
ZENODO_URL = "https://zenodo.org/api/records/16666299/files/NMD-rules.zip/content"
ZENODO_SIZE = 3_189_701_188
NONSENSE_LIST = "NMD-rules/Datasets/Nonsense_list.txt"
ISOFORMS = "NMD-rules/Datasets/IsoformLevel.NMDannotation.txt"
# zipfile checks the CRC32 of each member, and read_member checks the sha256
GEUVADIS_SHA256 = {
    NONSENSE_LIST: "5a248be88a6b2a287c2baf12e542d638257a06daadc21a4c440af812518bce26",
    ISOFORMS: "bc16f2b869149f07e0a5d7ce4173df4fbca9ea9c4fab2dd9f3c07a82c4011a7d",
}
# A Geuvadis variant with one of these flags in IsoformLevel.NMDannotation.txt is not in the accurate annotation
ACCURATE_ANNOTATION_FLAGS = ["resMNV", "nonTranslating", "Multiisoforms"]
# VAF_DNA of a heterozygous germline variant
GERMLINE_VAF_DNA = 0.5

SEED = 0
N_RESAMPLES = 2000
# The markdown cell before gtex_tissues and scripts/README.md state this number
MIN_TISSUE_ROWS = 100
# The names of the cohorts in the figures
COHORT_NAMES = {
    "tcga": "TCGA, out of fold",
    "geuvadis": "Geuvadis",
    "gtex": "GTEx, mean per variant",
    "mmrf_target": "MMRF-TARGET",
}
# The scores of interval_plot: the key in a result of correlations, the direct label, the color, and the offset from the
# position of the result
MODEL_AND_RULES = [("model", "Model", BLUE, 0.15), ("rules", "Canonical rules", ORANGE, -0.15)]
GAIN = [("gain", "", BLUE, 0.0)]
# The keys of the MMRF-TARGET rows by transcript in validation.json, and their names in the log and the figure
TRANSCRIPT_GROUPS = {
    "in_tcga_training": "Transcript in the TCGA training rows",
    "not_in_tcga_training": "Transcript not in the TCGA training rows",
}


class RemoteFile:
    """
    A read-only file at a URL that answers HTTP range requests. Each read is one request, so zipfile can read single
    members of a remote zip.
    """

    def __init__(self, url: str, size: int) -> None:
        self.url = url
        self.size = size
        self.position = 0

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.position

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        self.position = offset + {io.SEEK_SET: 0, io.SEEK_CUR: self.position, io.SEEK_END: self.size}[whence]
        return self.position

    def read(self, n: int = -1) -> bytes:
        end = self.size if n < 0 else min(self.position + n, self.size)
        if end <= self.position:
            return b""
        request = urllib.request.Request(self.url, headers={"Range": f"bytes={self.position}-{end - 1}"})
        with urllib.request.urlopen(request, timeout=120) as response:
            # A server that ignores the range answers 200 and sends the whole file
            if response.status != 206:
                raise ValueError(f"{self.url} answered a range request with HTTP {response.status}")
            data = response.read()
        self.position += len(data)
        return data


def read_member(member: str, geuvadis_dir: Path | None) -> io.BytesIO:
    """Read a member of the Zenodo zip from geuvadis_dir, or else by HTTP range requests, and check its sha256."""
    if geuvadis_dir is None:
        with zipfile.ZipFile(RemoteFile(ZENODO_URL, ZENODO_SIZE)) as archive:
            data = archive.read(member)
    else:
        data = (geuvadis_dir / Path(member).name).read_bytes()
    return io.BytesIO(check_sha256(data, GEUVADIS_SHA256[member], member))


def read_geuvadis(geuvadis_dir: Path | None, fasta: Fasta) -> pl.DataFrame:
    """
    Read the nonsense SNVs of Nonsense_list.txt, one row per variant and MANE Select transcript. Add the flags of
    IsoformLevel.NMDannotation.txt, which has the evaluated variants. REF and ALT are genomic. Fail if a REF does not
    match the FASTA.
    """
    table = pl.read_csv(read_member(NONSENSE_LIST, geuvadis_dir), separator="\t", infer_schema_length=None).rename(
        {"CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt", "TranscriptID": "transcript"}
    )
    fasta_ref = [
        fasta[chrom][pos - 1 : pos - 1 + len(ref)].seq.upper()
        for chrom, pos, ref in zip(table["chrom"], table["pos"], table["ref"], strict=True)
    ]
    if table["ref"].to_list() != fasta_ref:
        raise ValueError(f"rows of {NONSENSE_LIST} have a REF that does not match the FASTA")
    flags = pl.read_csv(
        read_member(ISOFORMS, geuvadis_dir), separator="\t", columns=["Variant_ID", *ACCURATE_ANNOTATION_FLAGS]
    )
    return table.join(flags, on="Variant_ID", how="left", validate="m:1", maintain_order="left")


def scored_rows(
    cohort: str, measured: pl.DataFrame, scored: pl.DataFrame, transcripts: pl.DataFrame, vaf_rna: str
) -> pl.DataFrame:
    """
    Join each measured row to the NMD-Scanner row of its transcript and variant. Keep the rows with
    nmd_efficiency_prediction_status "ok" and a variant allele fraction above 0 in the RNA, which is the column vaf_rna.
    Log the rows that this drops, and why.
    """
    rows = join_features(measured, scored, transcripts, cohort, logger)
    expressed = rows.filter(pl.col(vaf_rna) > 0)
    if len(expressed) < len(rows):
        logger.info(
            "%s: %d of the %d rows with nmd_efficiency_prediction_status ok have no %s above 0: %d lack it, %d have 0",
            cohort,
            len(rows) - len(expressed),
            len(rows),
            vaf_rna,
            rows[vaf_rna].null_count(),
            (rows[vaf_rna] == 0).sum(),
        )
    return expressed


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def percentiles(values: np.ndarray) -> list[float]:
    return [float(v) for v in np.percentile(values, [2.5, 97.5])]


def correlations(rows: pl.DataFrame) -> dict:
    """
    Spearman correlation of the prediction and of the rules with NMD_efficiency, and their difference, with bootstrap
    CIs. The rules score is (not nmd_escape), so that both scores rank high NMD efficiency first. Each resample draws
    the genes of the column gene with replacement, and it takes all rows of each drawn gene.
    """
    y = rows["NMD_efficiency"].to_numpy()
    model = rows["prediction"].to_numpy()
    rules = (~rows["nmd_escape"]).cast(pl.Float64).to_numpy()
    # The positions of the rows of each gene, in the order of the genes
    genes = rows.with_row_index().group_by("gene").agg("index").sort("gene")
    gene_rows = [positions.to_numpy() for positions in genes["index"]]
    rng = np.random.default_rng(SEED)
    model_resampled = np.empty(N_RESAMPLES)
    rules_resampled = np.empty(N_RESAMPLES)
    for resample in range(N_RESAMPLES):
        i = np.concatenate([gene_rows[gene] for gene in rng.integers(0, len(gene_rows), len(gene_rows))])
        model_resampled[resample] = spearman(model[i], y[i])
        rules_resampled[resample] = spearman(rules[i], y[i])
    model_value = spearman(model, y)
    rules_value = spearman(rules, y)
    return {
        "n": len(y),
        "genes": len(gene_rows),
        "model": model_value,
        "model_ci": percentiles(model_resampled),
        "rules": rules_value,
        "rules_ci": percentiles(rules_resampled),
        "gain": model_value - rules_value,
        "gain_ci": percentiles(model_resampled - rules_resampled),
    }


def mean_nmd_triggering(rows: pl.DataFrame) -> float:
    """The mean NMD_efficiency of the rows with nmd_escape False."""
    # np.mean, because polars mean sums in another order and can differ in the last digit
    return float(np.mean(rows.filter(~pl.col("nmd_escape"))["NMD_efficiency"].to_numpy()))


def estimate(value: float, ci: list[float]) -> str:
    return f"{value:.2f} ({ci[0]:.2f} to {ci[1]:.2f})"


def facet_frame(panels: dict[str, pl.DataFrame]) -> pl.DataFrame:
    """
    Concatenate the frames of panels for facet_wrap, with the columns of all frames. The column panel has their keys as
    categories, in order.
    """
    frames = [frame.with_columns(panel=pl.lit(title)) for title, frame in panels.items()]
    return pl.concat(frames, how="diagonal_relaxed").cast({"panel": pl.Enum(list(panels))})


def intervals(results: list[dict], scores: list[tuple[str, str, str, float]]) -> pl.DataFrame:
    """
    One row per result and score for interval_plot, with the value of the score and its 95% CI as low and high. Result
    i has the position -i plus the offset of the score. The rows of the last result, at the bottom, have the direct
    labels.
    """
    return pl.DataFrame(
        [
            {
                "position": -i + offset,
                "value": result[score],
                "low": result[f"{score}_ci"][0],
                "high": result[f"{score}_ci"][1],
                "color": color,
                "label": label if i == len(results) - 1 else "",
            }
            for i, result in enumerate(results)
            for score, label, color, offset in scores
        ]
    )


def interval_plot(rows: pl.DataFrame, labels: list[str]) -> ggplot:
    """
    A dot at each value of the intervals in rows, with a line from low to high, and the direct label right of high.
    The y axis shows labels from the top, one per result.
    """
    return (
        ggplot(rows, aes(y="position", color="color"))
        + geom_segment(aes(x="low", xend="high", yend="position"), size=1)
        + geom_point(aes(x="value"), size=2.5)
        + geom_text(aes(x="high", label="label"), ha="left", nudge_x=0.01, color=SECONDARY_INK, size=9)
        + scale_color_identity()
        + scale_y_continuous(breaks=-np.arange(len(labels)), labels=labels, expand=(0, 0.4))
        + theme(axis_title_y=element_blank(), panel_grid_major_y=element_blank())
    )


# %%
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--gff3", required=True, type=Path, help="GENCODE GFF3 of GRCh38, optionally gzip-compressed")
parser.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA with chr names and a .fai index")
parser.add_argument(
    "--train-dir",
    required=True,
    type=Path,
    help="--out-dir of the train_model.py run that trained the model; it has the out-of-fold predictions",
)
parser.add_argument("--out-dir", required=True, type=Path, help="directory for the outputs")
parser.add_argument(
    "--geuvadis-dir",
    type=Path,
    help="local directory with Nonsense_list.txt and IsoformLevel.NMDannotation.txt of NMD-rules.zip; skips the"
    " range requests to Zenodo. The sha256 of both files is still checked.",
)
# In a notebook, sys.argv holds the arguments of the kernel
args = parser.parse_args(shlex.split(os.environ.get("VALIDATE_MODEL_ARGS", "")) if in_notebook() else None)
logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s")
logger.setLevel(logging.INFO)
args.out_dir.mkdir(parents=True, exist_ok=True)
trained = args.train_dir / "nmd_efficiency_rf.onnx"
if trained.read_bytes() != MODEL_PATH.read_bytes():
    raise ValueError(f"{trained} differs from {MODEL_PATH}, so its out-of-fold predictions are of another model")

# %%
transcripts = read_transcripts(args.gff3)
fasta = Fasta(str(args.fasta))
geuvadis = read_geuvadis(args.geuvadis_dir, fasta)
# The first line of the sheet is its title. polars reads the sheet with fastexcel.
gtex = pl.read_excel(download(GTEX_URL, GTEX_SHA256), sheet_name=GTEX_SHEET, read_options={"header_row": 1})
gtex = add_alleles(gtex, transcripts, fasta).rename({"Tissue_type": "tissue"})
mmrf = add_alleles(read_nmdeff("MMRF_TARGET_dataset.csv"), transcripts, fasta)
vcf = args.out_dir / "variants.vcf"
variants = pl.concat([table.select("chrom", "pos", "ref", "alt") for table in (geuvadis, gtex, mmrf)])
write_vcf(variants.unique(keep="first", maintain_order=True), vcf)
features = add_key(
    pl.from_arrow(nmd_scanner.to_arrow(nmd_scanner.annotate(vcf, args.gff3, args.fasta, sequences=False)))
)
scored = features.select(*KEY, "nmd_efficiency_prediction_status", "nmd_escape", prediction="nmd_efficiency_prediction")

geuvadis = scored_rows("geuvadis", geuvadis, scored, transcripts, "VAF_mean")
accurate = geuvadis.select(~pl.any_horizontal(pl.col(ACCURATE_ANNOTATION_FLAGS).fill_null(False))).to_series()
logger.info(
    "geuvadis: %d of the %d rows are outside the study's accurate annotation (flags %s)",
    (~accurate).sum(),
    len(geuvadis),
    " ".join(ACCURATE_ANNOTATION_FLAGS),
)
# np.log2, because polars log(2) can differ from it in the last digit
geuvadis = geuvadis.filter(accurate).with_columns(NMD_efficiency=-np.log2(pl.col("VAF_mean") / GERMLINE_VAF_DNA))
# ASE_altered is VAF_RNA. For VAF_RNA 0, the table caps NMD_efficiency at 10.
gtex = scored_rows("gtex", gtex, scored, transcripts, "ASE_altered")
# group_by has no fixed order, and the bootstrap of correlations depends on the order of the rows
gtex_variants = (
    gtex.group_by(KEY).agg(pl.col("NMD_efficiency").mean(), pl.col("prediction", "nmd_escape").first()).sort(KEY)
)
cohorts = {
    # Every row of tcga_dataset.csv has VAF_RNA above 0
    "tcga": pl.read_parquet(args.train_dir / "oof_predictions.parquet"),
    "geuvadis": geuvadis,
    "gtex": gtex_variants,
    "mmrf_target": scored_rows("mmrf_target", mmrf, scored, transcripts, "VAF_RNA"),
}
# correlations resamples the genes
gene = pl.col("transcript").replace_strict(transcripts["transcript"], transcripts["gene"])
cohorts = {cohort: rows.with_columns(gene=gene) for cohort, rows in cohorts.items()}

# %%
validation = {}
for cohort, rows in cohorts.items():
    result = {**correlations(rows), "mean_nmd_triggering": mean_nmd_triggering(rows)}
    logger.info(
        "%s, %d rows of %d genes: model %s, rules %s, gain %s, NMD-triggering mean %.2f",
        cohort,
        result["n"],
        result["genes"],
        estimate(result["model"], result["model_ci"]),
        estimate(result["rules"], result["rules_ci"]),
        estimate(result["gain"], result["gain_ci"]),
        result["mean_nmd_triggering"],
    )
    validation[cohort] = result
tissue_means = {
    tissue: mean_nmd_triggering(rows)
    for (tissue,), rows in gtex.group_by("tissue", maintain_order=True)
    if len(rows) >= MIN_TISSUE_ROWS
}
tissue_means = dict(sorted(tissue_means.items(), key=lambda item: item[1]))
validation["gtex_variant_tissue_rows"] = {
    "n": len(gtex),
    "tissues": gtex["tissue"].n_unique(),
    "mean_nmd_triggering": mean_nmd_triggering(gtex),
    "min_rows_per_tissue": MIN_TISSUE_ROWS,
    "mean_nmd_triggering_per_tissue": tissue_means,
}
mmrf_rows = cohorts["mmrf_target"]
seen = mmrf_rows["transcript"].is_in(cohorts["tcga"]["transcript"].implode())
validation["mmrf_target_by_transcript"] = {
    "in_tcga_training": correlations(mmrf_rows.filter(seen)),
    "not_in_tcga_training": correlations(mmrf_rows.filter(~seen)),
}
for group, result in validation["mmrf_target_by_transcript"].items():
    logger.info(
        "mmrf_target, %s, %d rows of %d genes: model %s, rules %s, gain %s",
        TRANSCRIPT_GROUPS[group],
        result["n"],
        result["genes"],
        estimate(result["model"], result["model_ci"]),
        estimate(result["rules"], result["rules_ci"]),
        estimate(result["gain"], result["gain_ci"]),
    )
(args.out_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
logger.info(
    "gtex, %d variant-tissue rows: NMD-triggering mean %.2f, from %.2f to %.2f between the %d tissues with at"
    " least %d rows",
    len(gtex),
    mean_nmd_triggering(gtex),
    min(tissue_means.values()),
    max(tissue_means.values()),
    len(tissue_means),
    MIN_TISSUE_ROWS,
)

# %% [markdown]
# The scatters show the measured NMD efficiency of each row over the prediction of the model, per cohort. The gray
# line marks equal values.

# %%
scatter = facet_frame(
    {
        f"{COHORT_NAMES[cohort]}, {len(rows):,} rows\nSpearman {validation[cohort]['model']:.2f}": rows
        for cohort, rows in cohorts.items()
    }
)
plot = (
    ggplot(scatter, aes("prediction", "NMD_efficiency"))
    + geom_point(color=BLUE, alpha=0.3, size=1.5, stroke=0)
    + geom_abline(color=AXIS)
    + facet_wrap("panel", scales="free")
    + labs(x="Predicted NMD efficiency", y="Measured NMD efficiency")
)
save_figure(plot, args.out_dir, "measured_vs_predicted", 8, 7.5)

# %% [markdown]
# Per cohort, the dots show the mean measured over the mean predicted NMD efficiency in each quantile bin of the
# prediction, with the 95% CI of the measured mean. The gray line marks equal values, and the blue line is the
# least-squares fit of measured on predicted over all rows, the linear recalibration that the main README suggests.

# %%
calibrations = {}
for cohort, rows in cohorts.items():
    bins, slope, intercept = calibration(rows)
    title = f"{COHORT_NAMES[cohort]}, {len(rows):,} rows\nFit: slope {slope:.2f}, intercept {intercept:.2f}"
    calibrations[title] = bins.with_columns(slope=slope, intercept=intercept)
plot = calibration_plot(facet_frame(calibrations)) + facet_wrap("panel", scales="free")
save_figure(plot, args.out_dir, "calibration", 8, 9)

# %% [markdown]
# The boxes show the measured NMD efficiency of the rows by `nmd_escape`, per cohort. A box spans the quartiles, and
# its whiskers reach the furthest row within 1.5 times the interquartile range.

# %%
boxes = facet_frame(
    {
        COHORT_NAMES[cohort]: rows.with_columns(
            pl.col("nmd_escape").replace_strict(
                {escape: f"{escape}\n{count:,} rows" for escape, count in rows["nmd_escape"].value_counts().iter_rows()}
            )
        )
        for cohort, rows in cohorts.items()
    }
)
plot = (
    ggplot(boxes, aes("nmd_escape", "NMD_efficiency"))
    + geom_boxplot(color=BLUE, fill=BLUE, alpha=0.1, width=0.5, outlier_size=1, outlier_alpha=0.3, outlier_stroke=0)
    + facet_wrap("panel", nrow=1, scales="free")
    + labs(x="nmd_escape", y="Measured NMD efficiency")
    + theme(panel_grid_major_x=element_blank())
)
save_figure(plot, args.out_dir, "efficiency_by_nmd_escape", 9.5, 3.5)

# %% [markdown]
# The dots show the Spearman correlation of the model and of the canonical rules with the measured NMD efficiency, per
# cohort, with the 95% bootstrap CIs. These are the numbers of the table in the section "Model validation" of the main
# README.

# %%
plot = interval_plot(
    intervals([validation[cohort] for cohort in cohorts], MODEL_AND_RULES),
    [f"{COHORT_NAMES[cohort]}\n{validation[cohort]['n']:,} rows" for cohort in cohorts],
) + labs(title="Spearman correlation with the measured NMD efficiency, 95% CIs", x="Spearman correlation")
save_figure(plot, args.out_dir, "spearman_per_cohort", 6.5, 3.5)

# %% [markdown]
# The dots show the mean measured NMD efficiency of the GTEx variant-tissue rows with `nmd_escape` False, per tissue
# with at least 100 variant-tissue rows, with the 95% CI of the mean. The gray line marks the mean over all
# variant-tissue rows with `nmd_escape` False.

# %%
tissue_rows = validation["gtex_variant_tissue_rows"]
# The tissues from the lowest mean to the highest, which coord_flip draws from the bottom to the top
tissues = list(tissue_rows["mean_nmd_triggering_per_tissue"])
nmd_triggering = gtex.filter(~pl.col("nmd_escape"), pl.col("tissue").is_in(tissues)).cast({"tissue": pl.Enum(tissues)})
plot = (
    ggplot(nmd_triggering, aes("tissue", "NMD_efficiency"))
    + geom_hline(yintercept=tissue_rows["mean_nmd_triggering"], color=AXIS, linetype="dashed")
    + annotate(
        "text",
        x=tissues[-1],
        y=tissue_rows["mean_nmd_triggering"],
        label=f"All tissues {tissue_rows['mean_nmd_triggering']:.2f}",
        ha="left",
        nudge_y=0.01,
        color=SECONDARY_INK,
        size=9,
    )
    + stat_summary(fun_data="mean_cl_normal", color=BLUE)
    + coord_flip()
    + labs(
        title=f"GTEx, NMD-triggering variant-tissue rows of {len(tissues)} tissues",
        x="",
        y="Mean measured NMD efficiency, 95% CI",
    )
    + theme(panel_grid_major_y=element_blank())
)
save_figure(plot, args.out_dir, "gtex_tissues", 6, 9)

# %% [markdown]
# The dots show the Spearman correlations and the gain of the model over the rules in MMRF-TARGET, for the rows whose
# transcript is in the TCGA training rows of --train-dir and for the other rows, with 95% bootstrap CIs. The numbers
# are `mmrf_target_by_transcript` in validation.json.

# %%
transcript_split = validation["mmrf_target_by_transcript"]
results = list(transcript_split.values())
gain_panel = "Gain of the model over the rules"
transcript_intervals = facet_frame(
    {"Spearman correlation": intervals(results, MODEL_AND_RULES), gain_panel: intervals(results, GAIN)}
)
zero = pl.DataFrame({"panel": [gain_panel], "x": [0.0]}).cast({"panel": transcript_intervals["panel"].dtype})
plot = (
    interval_plot(
        transcript_intervals,
        [f"{TRANSCRIPT_GROUPS[group]}\n{result['n']:,} rows" for group, result in transcript_split.items()],
    )
    + geom_vline(aes(xintercept="x"), data=zero, color=AXIS)
    + facet_wrap("panel", scales="free_x")
    + labs(title="MMRF-TARGET by transcript, 95% CIs", x="")
)
save_figure(plot, args.out_dir, "mmrf_target_by_transcript", 9, 3)
