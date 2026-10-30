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

# %% [markdown]
# # Model validation
#
# This notebook draws the figures of the training run of scripts/train_model.py, and it computes the numbers of the
# section "Model validation" of the main README: how well nmd_efficiency_rf.onnx and nmd_escape alone rank the measured
# NMD efficiency in the benchmarks of nmd_scanner.benchmark. It writes validation.json and figures/, with each figure
# as PNG and PDF, to OUT_DIR. scripts/README.md shows how to run it.
#
# Cohorts: the benchmarks of nmd_scanner.benchmark.DATASETS. Each module of nmd_scanner.benchmark.datasets describes
# its study and its NMD efficiency. TCGA has the out-of-fold predictions in oof_predictions.parquet of TRAIN_DIR.
#
# Steps of the validation:
#
# 1. Read the benchmark tables of BENCHMARK_DIR, except TCGA.
# 2. Run nmd_scanner.annotate_variants on their unique variants. It adds the prediction of the model as
#    nmd_efficiency_prediction. Join each measurement to the NMD-Scanner result of its variant-transcript pair.
# 3. Keep the measurements with nmd_efficiency_prediction_status "ok". Log the measurements that this drops, and why.
# 4. A benchmark with the column tissue has one row per variant and tissue. Its cohort has one row per variant, with
#    the mean NMD efficiency over the tissues.
# 5. Per cohort, compute the Spearman correlation of the predictions and of (not nmd_escape) with the measured NMD
#    efficiency. The 95% CIs are percentiles of N_RESAMPLES bootstrap resamples of the genes. A resample draws genes
#    with replacement and takes all measurements of each drawn gene. Measurements of one gene share a transcript
#    structure, so they are not independent. The gene of a measurement is the gene of its transcript in the GFF3. The
#    gain, model minus rules, has a paired CI from the same resamples. Its P value is from a two-sided z-test, with the
#    standard deviation of the resampled gains as the standard error. Each cohort draws from its own generator, seeded
#    with SEED.
# 6. Per cohort, compute the mean measured NMD efficiency of the measurements with nmd_escape False. For GTEx, also
#    compute it over the variant-tissue pairs, in total and per tissue. Per tissue, also compute the Spearman
#    correlations of step 5, without CIs. The per-tissue numbers leave out the EXCLUDED_TISSUES, which have too few
#    samples.
# 7. For MMRF-TARGET, compute the numbers of step 5 also for the variant-sample pairs on transcripts of the TCGA
#    training variant-sample pairs, and for the other variant-sample pairs.
# 8. Per cohort, split the r2_score of scikit-learn for the predictions into the squared Pearson correlation, minus a
#    term for the offset of the means, minus a term for the mismatch of the spreads.

# %% [markdown]
# The paths, relative to the repo root:

# %%
from pathlib import Path

GFF3 = Path("data/gencode.v42.annotation.gff3.gz")
FASTA = Path("data/GRCh38.fa")
# The --out-dir of scripts/build_benchmarks.py, with the same GFF3 and FASTA
BENCHMARK_DIR = Path("data/benchmarks")
# The --out-dir of the scripts/train_model.py run that trained nmd_efficiency_rf.onnx
TRAIN_DIR = Path("data/train")
OUT_DIR = Path("data/validate")

# %%
import json
import logging

import numpy as np
import plotnine as pn
import polars as pl
from pyfaidx import Fasta
from scipy.stats import norm, rankdata
from sklearn.metrics import r2_score

from nmd_scanner import annotate_variants, read_annotation
from nmd_scanner.benchmark import DATASETS, KEY, VARIANT_COLUMNS, join_features
from nmd_scanner.benchmark.figures import (
    AXIS,
    BLUE,
    MUTED_INK,
    ORANGE,
    SECONDARY_INK,
    SURFACE,
    calibration,
    calibration_plot,
    save_figure,
)
from nmd_scanner.model import MODEL_PATH

logger = logging.getLogger("validate_model")

# %%
SEED = 0
N_RESAMPLES = 2000
# The tissues to drop, because GTEx v8 has too few samples of them (see the AbExp paper)
EXCLUDED_TISSUES = ["Bladder", "Cervix - Ectocervix", "Cervix - Endocervix", "Fallopian Tube", "Kidney - Medulla"]
# The scores of the Spearman bars and of the MMRF-TARGET intervals: the key in the dict of correlations, the direct
# label, the color, and the offset from the position of the group
MODEL_AND_RULES = [("model", "Model", BLUE, 0.15), ("rules", "Canonical rules", ORANGE, -0.15)]
GAIN = [("gain", "", BLUE, 0.0)]
# The keys of the MMRF-TARGET variant-sample pairs by transcript in validation.json, and their names in the log and the
# figure
TRANSCRIPT_GROUPS = {
    "in_tcga_training": "Transcript in the TCGA training variant-sample pairs",
    "not_in_tcga_training": "Transcript not in the TCGA training variant-sample pairs",
}
# The name of the measurements of each cohort in the figures and the log. A measurement of the GTEx cohort is the mean
# of one variant over the tissues.
MEASUREMENTS = {
    "tcga": "variant-sample pairs",
    "geuvadis": "variants",
    "gtex": "variants",
    "mmrf_target": "variant-sample pairs",
}
# The terms of r2_terms: the key, the sign in the sum, and the name in the figure
R2_TERMS = [
    ("r_squared", 1, "Squared Pearson r"),
    ("mean_offset", -1, "Minus mean offset"),
    ("spread_mismatch", -1, "Minus spread mismatch"),
    ("r2_score", 1, "= R2"),
]


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def percentiles(values: np.ndarray) -> list[float]:
    return [float(v) for v in np.percentile(values, [2.5, 97.5])]


def correlations(cohort_df: pl.DataFrame) -> dict:
    """
    Spearman correlation of the prediction and of the rules with NMD_efficiency, and their difference, with bootstrap
    CIs. The rules score is (not nmd_escape), so that both scores rank high NMD efficiency first. Each resample draws
    the genes of the column gene with replacement, and it takes all measurements of each drawn gene.
    """
    y = cohort_df["NMD_efficiency"].to_numpy()
    model = cohort_df["prediction"].to_numpy()
    rules = (~cohort_df["nmd_escape"]).cast(pl.Float64).to_numpy()
    # The positions of the measurements of each gene, in the order of the genes
    genes_df = cohort_df.with_row_index().group_by("gene").agg("index").sort("gene")
    gene_positions = [positions.to_numpy() for positions in genes_df["index"]]
    rng = np.random.default_rng(SEED)
    model_resampled = np.empty(N_RESAMPLES)
    rules_resampled = np.empty(N_RESAMPLES)
    for resample in range(N_RESAMPLES):
        i = np.concatenate([gene_positions[gene] for gene in rng.integers(0, len(gene_positions), len(gene_positions))])
        model_resampled[resample] = spearman(model[i], y[i])
        rules_resampled[resample] = spearman(rules[i], y[i])
    model_value = spearman(model, y)
    rules_value = spearman(rules, y)
    gain = model_value - rules_value
    gain_resampled = model_resampled - rules_resampled
    return {
        "n": len(y),
        "genes": len(gene_positions),
        "model": model_value,
        "model_ci": percentiles(model_resampled),
        "rules": rules_value,
        "rules_ci": percentiles(rules_resampled),
        "gain": gain,
        "gain_ci": percentiles(gain_resampled),
        "gain_p": float(2 * norm.sf(abs(gain) / gain_resampled.std(ddof=1))),
    }


def mean_nmd_triggering(cohort_df: pl.DataFrame) -> float:
    """The mean NMD_efficiency of the measurements with nmd_escape False."""
    return float(cohort_df.filter(~pl.col("nmd_escape"))["NMD_efficiency"].mean())


def r2_terms(cohort_df: pl.DataFrame) -> dict[str, float]:
    """r2_score = r_squared - mean_offset - spread_mismatch, as the R2 figure explains. ddof=0 makes the sum exact."""
    r, offset, sd_prediction, sd_measured = cohort_df.select(
        r=pl.corr("prediction", "NMD_efficiency"),
        offset=pl.col("prediction").mean() - pl.col("NMD_efficiency").mean(),
        sd_prediction=pl.col("prediction").std(ddof=0),
        sd_measured=pl.col("NMD_efficiency").std(ddof=0),
    ).row(0)
    terms = {
        "r_squared": r**2,
        "mean_offset": offset**2 / sd_measured**2,
        "spread_mismatch": (sd_prediction - r * sd_measured) ** 2 / sd_measured**2,
        "r2_score": float(r2_score(cohort_df["NMD_efficiency"].to_numpy(), cohort_df["prediction"].to_numpy())),
    }
    difference = terms["r_squared"] - terms["mean_offset"] - terms["spread_mismatch"] - terms["r2_score"]
    assert abs(difference) < 1e-9, terms
    return terms


def estimate(value: float, ci: list[float]) -> str:
    return f"{value:.2f} ({ci[0]:.2f} to {ci[1]:.2f})"


def p_value(p: float) -> str:
    """P with 2 significant digits, and below 0.001 as a power of 10, such as P = 1.6×10⁻⁷."""
    if p >= 0.001:
        return f"P = {p:.2g}"
    mantissa, exponent = f"{p:.1e}".split("e")
    return f"P = {mantissa}×10{str(int(exponent)).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))}"


def facet_frame(panels: dict[str, pl.DataFrame]) -> pl.DataFrame:
    """
    Concatenate the frames of panels for facet_wrap, with the columns of all frames. The column panel has their keys as
    categories, in order.
    """
    frames = [frame_df.with_columns(panel=pl.lit(title)) for title, frame_df in panels.items()]
    return pl.concat(frames, how="diagonal_relaxed").cast({"panel": pl.Enum(list(panels))})


def intervals(correlations_by_group: dict[str, dict], scores: list[tuple[str, str, str, float]]) -> pl.DataFrame:
    """
    One row per group of correlations_by_group and score, with the value of the score and its 95% CI as low and high.
    Group i has the position -i plus the offset of the score. The rows of the last group, at the bottom, have the direct
    labels.
    """
    return pl.DataFrame(
        [
            {
                "position": -i + offset,
                "value": group_correlations[score],
                "low": group_correlations[f"{score}_ci"][0],
                "high": group_correlations[f"{score}_ci"][1],
                "color": color,
                "label": label if i == len(correlations_by_group) - 1 else "",
            }
            for i, group_correlations in enumerate(correlations_by_group.values())
            for score, label, color, offset in scores
        ]
    )


# %%
logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s")
for name in ("validate_model", "nmd_scanner.benchmark"):
    logging.getLogger(name).setLevel(logging.INFO)
OUT_DIR.mkdir(parents=True, exist_ok=True)
trained = TRAIN_DIR / "nmd_efficiency_rf.onnx"
if trained.read_bytes() != MODEL_PATH.read_bytes():
    raise ValueError(f"{trained} differs from {MODEL_PATH}, so its out-of-fold predictions are of another model")

# %% [markdown]
# ## Training

# %% [markdown]
# The histogram shows the target `NMD_efficiency` of the TCGA variant-sample pairs. 0 means no NMD, and 1 means that the
# variant allele fraction in the RNA is half of that in the DNA.

# %%
oof_df = pl.read_parquet(TRAIN_DIR / "oof_predictions.parquet")
plot = (
    pn.ggplot(oof_df, pn.aes("NMD_efficiency"))
    + pn.geom_histogram(bins=60, fill=BLUE, color=SURFACE, size=0.25)
    + pn.labs(
        title=f"NMD_efficiency of {len(oof_df):,} TCGA variant-sample pairs",
        x="NMD_efficiency",
        y="Variant-sample pairs",
    )
    + pn.theme(panel_grid_major_x=pn.element_blank())
)
save_figure(plot, OUT_DIR, "target_histogram", 6, 3.5)

# %% [markdown]
# The scatter shows the measured NMD efficiency of each TCGA variant-sample pair over its out-of-fold prediction, with
# Spearman and R2 over all variant-sample pairs. The gray line marks equal values.

# %%
oof_measured = oof_df["NMD_efficiency"].to_numpy()
oof_spearman = spearman(oof_df["prediction"].to_numpy(), oof_measured)
oof_r2 = r2_score(oof_measured, oof_df["prediction"].to_numpy())
plot = (
    pn.ggplot(oof_df, pn.aes("prediction", "NMD_efficiency"))
    + pn.geom_point(color=BLUE, alpha=0.3, size=1.5, stroke=0)
    + pn.geom_abline(color=AXIS)
    + pn.labs(
        title=f"TCGA, out of fold, {len(oof_df):,} variant-sample pairs\nSpearman {oof_spearman:.2f}, R2 {oof_r2:.2f}",
        x="Predicted NMD efficiency",
        y="Measured NMD efficiency",
    )
)
save_figure(plot, OUT_DIR, "tcga_measured_vs_predicted", 5, 4.5)

# %% [markdown]
# The dots show the mean measured over the mean predicted NMD efficiency in each quantile bin of the out-of-fold
# prediction, with the 95% CI of the measured mean. The blue line is the least-squares fit of measured on predicted over
# all variant-sample pairs, and the gray line marks equal values.

# %%
bins_df, slope, intercept = calibration(oof_df)
plot = calibration_plot(bins_df.with_columns(slope=slope, intercept=intercept)) + pn.labs(
    title=f"TCGA, out of fold, {len(oof_df):,} variant-sample pairs\nFit: slope {slope:.2f}, intercept {intercept:.2f}"
)
save_figure(plot, OUT_DIR, "calibration_oof", 5, 6)

# %% [markdown]
# The bars show the impurity importance and the permutation importance of the model inputs in the final random forest.
# The permutation importance is the loss in R2 on the training variant-sample pairs when one input is permuted, with the
# standard deviation over the permutations.

# %%
importance_df = pl.read_parquet(TRAIN_DIR / "feature_importance.parquet").sort("value", "input")
# Each panel orders its bars by importance, so a bar is one kind and input. Bars of equal importance are ordered
# by input name. The axis shows only the input.
bars = importance_df["kind"] + "|" + importance_df["input"]
importance_df = importance_df.with_columns(bar=bars.cast(pl.Enum(bars)))
plot = (
    pn.ggplot(importance_df, pn.aes("bar", "value"))
    + pn.geom_col(fill=BLUE, width=0.6)
    + pn.geom_linerange(
        pn.aes(ymin="value - sd", ymax="value + sd"), data=importance_df.drop_nulls(), color=SECONDARY_INK
    )
    + pn.scale_x_discrete(labels=lambda bars: [bar.split("|")[1] for bar in bars])
    + pn.coord_flip()
    + pn.facet_wrap("kind", scales="free")
    + pn.labs(x="", y="Importance")
    + pn.theme(panel_grid_major_y=pn.element_blank())
)
save_figure(plot, OUT_DIR, "feature_importance", 9, 5)

# %% [markdown]
# The learning curve shows the R2 of the final hyperparameters on the training variant-sample pairs and on the held-out
# fold, for growing numbers of training variant-sample pairs, with bands of one standard deviation over the folds. Each
# chromosome fold trains on the same numbers of variant-sample pairs, up to all training variant-sample pairs of the
# fold with the fewest.

# %%
learning_df = pl.read_parquet(TRAIN_DIR / "learning_curve.parquet").sort("training_variant_sample_pairs", "fold")
train_sizes = learning_df["training_variant_sample_pairs"].unique().sort().to_numpy()
# The scores as learning_curve returns them: one row per training size, and one column per fold
train_scores, test_scores = (
    learning_df[column].to_numpy().reshape(len(train_sizes), -1) for column in ("training_r2", "held_out_r2")
)
curves_df = pl.concat(
    [
        pl.DataFrame(
            {
                "training_variant_sample_pairs": train_sizes,
                "mean": scores.mean(axis=1),
                "sd": scores.std(axis=1),
                "label": label,
                "color": color,
            }
        )
        for scores, label, color in (
            (train_scores, "Training\nvariant-sample pairs", MUTED_INK),
            (test_scores, "Held-out fold", BLUE),
        )
    ]
)
training_df, held_out_df = (curves_df.filter(pl.col("color") == color) for color in (MUTED_INK, BLUE))
# The label of the training curve goes on two lines centered above its fifth point, clear of the band there and at the
# neighbouring points
band_top = training_df[3:6].select((pl.col("mean") + pl.col("sd")).max()).item()
training_label_df = training_df[4].with_columns(y=band_top + 0.003)
plot = (
    pn.ggplot(curves_df, pn.aes("training_variant_sample_pairs", "mean"))
    + pn.geom_ribbon(pn.aes(ymin="mean - sd", ymax="mean + sd", fill="color"), alpha=0.1)
    + pn.geom_line(pn.aes(color="color"))
    + pn.geom_point(pn.aes(color="color"), size=2.5)
    + pn.geom_text(
        pn.aes(y="y", label="label"), data=training_label_df, ha="center", va="bottom", color=SECONDARY_INK, size=9
    )
    + pn.geom_text(
        pn.aes(label="label"),
        data=held_out_df.tail(1),
        ha="left",
        nudge_x=0.03 * train_sizes[-1],
        color=SECONDARY_INK,
        size=9,
    )
    + pn.scale_color_identity()
    + pn.scale_fill_identity()
    + pn.scale_x_continuous(expand=(0.05, 0, 0.2, 0))
    + pn.labs(title="Learning curve of the final hyperparameters", x="Training variant-sample pairs", y="R2")
)
save_figure(plot, OUT_DIR, "learning_curve", 6, 4)

# %% [markdown]
# The curve shows the out-of-bag error of the final hyperparameters, fitted on all TCGA variant-sample pairs with a
# growing number of trees. The out-of-bag prediction of a variant-sample pair is the mean prediction of the trees that
# did not draw it, and the out-of-bag error is 1 - R2 of these predictions.

# %%
oob_df = pl.read_parquet(TRAIN_DIR / "oob_error.parquet").select("trees", error=1 - pl.col("oob_r2"))
final_trees = json.loads((TRAIN_DIR / "cv_results.json").read_text())["params"]["n_estimators"]
plot = (
    pn.ggplot(oob_df, pn.aes("trees", "error"))
    + pn.geom_vline(xintercept=final_trees, color=AXIS)
    + pn.annotate(
        "text",
        x=final_trees,
        y=np.inf,
        label=f"Final model, {final_trees} trees",
        ha="left",
        va="top",
        nudge_x=5,
        color=SECONDARY_INK,
        size=9,
    )
    + pn.geom_line(color=BLUE)
    + pn.geom_point(color=BLUE, size=2.5)
    + pn.labs(title="Out-of-bag error of the final hyperparameters", x="Trees", y="1 - out-of-bag R2")
)
save_figure(plot, OUT_DIR, "oob_error", 6, 4)

# %% [markdown]
# ## Validation

# %%
benchmarks = {name: pl.read_parquet(BENCHMARK_DIR / f"{name}.parquet") for name in DATASETS if name != "tcga"}
variants_df = pl.concat([table_df.select(VARIANT_COLUMNS) for table_df in benchmarks.values()]).unique()
fasta = Fasta(str(FASTA))
results_df = annotate_variants(variants_df.to_pandas(), read_annotation(GFF3, fasta), fasta, sequences=False)
# The bootstrap of correlations and the sums of the means depend on the order of the rows
joined = {
    name: join_features(table_df, results_df, DATASETS[name].NAME)
    .select(*table_df.columns, "nmd_escape", prediction="nmd_efficiency_prediction")
    .sort("source_row")
    for name, table_df in benchmarks.items()
}
# The cohorts and their names in the figures
cohorts = {"tcga": oof_df}
cohort_names = {"tcga": f"{DATASETS['tcga'].NAME}, out of fold"}
for name, cohort_df in joined.items():
    cohort_names[name] = DATASETS[name].NAME
    if "tissue" in cohort_df.columns:
        cohort_df = (
            cohort_df.group_by(KEY)
            .agg(
                pl.col("NMD_efficiency").mean(),
                pl.col("gene", "prediction", "nmd_escape").first(),
                pl.col("source_row").min(),
            )
            .sort("source_row")
        )
        cohort_names[name] += ", mean per variant"
    cohorts[name] = cohort_df
# The variant-tissue pairs of GTEx
gtex_df = joined["gtex"]

# %%
validation = {}
for cohort, cohort_df in cohorts.items():
    cohort_correlations = correlations(cohort_df)
    validation[cohort] = {
        **cohort_correlations,
        "mean_nmd_triggering": mean_nmd_triggering(cohort_df),
        "r2_decomposition": r2_terms(cohort_df),
    }
    logger.info(
        "%s, %d %s of %d genes: model %s, rules %s, gain %s, %s, NMD-triggering mean %.2f",
        cohort,
        cohort_correlations["n"],
        MEASUREMENTS[cohort],
        cohort_correlations["genes"],
        estimate(cohort_correlations["model"], cohort_correlations["model_ci"]),
        estimate(cohort_correlations["rules"], cohort_correlations["rules_ci"]),
        estimate(cohort_correlations["gain"], cohort_correlations["gain_ci"]),
        p_value(cohort_correlations["gain_p"]),
        validation[cohort]["mean_nmd_triggering"],
    )
tissues_df = (
    gtex_df.filter(~pl.col("tissue").is_in(EXCLUDED_TISSUES))
    .group_by("tissue")
    .agg(
        n=pl.len(),
        mean=pl.col("NMD_efficiency").filter(~pl.col("nmd_escape")).mean(),
        model=pl.corr("prediction", "NMD_efficiency", method="spearman"),
        rules=pl.corr(~pl.col("nmd_escape"), "NMD_efficiency", method="spearman"),
    )
)
tissue_means = dict(tissues_df.sort("mean", "tissue").select("tissue", "mean").iter_rows())
tissue_spearman_df = tissues_df.drop("mean").sort("model", "tissue")
validation["gtex_variant_tissue_pairs"] = {
    "n": len(gtex_df),
    "tissues": gtex_df["tissue"].n_unique(),
    "mean_nmd_triggering": mean_nmd_triggering(gtex_df),
    "excluded_tissues": EXCLUDED_TISSUES,
    "mean_nmd_triggering_per_tissue": tissue_means,
}
mmrf_df = cohorts["mmrf_target"]
seen = mmrf_df["transcript"].is_in(cohorts["tcga"]["transcript"].implode())
validation["mmrf_target_by_transcript"] = {
    "in_tcga_training": correlations(mmrf_df.filter(seen)),
    "not_in_tcga_training": correlations(mmrf_df.filter(~seen)),
}
for group, group_correlations in validation["mmrf_target_by_transcript"].items():
    logger.info(
        "mmrf_target, %s, %d variant-sample pairs of %d genes: model %s, rules %s, gain %s, %s",
        TRANSCRIPT_GROUPS[group],
        group_correlations["n"],
        group_correlations["genes"],
        estimate(group_correlations["model"], group_correlations["model_ci"]),
        estimate(group_correlations["rules"], group_correlations["rules_ci"]),
        estimate(group_correlations["gain"], group_correlations["gain_ci"]),
        p_value(group_correlations["gain_p"]),
    )
(OUT_DIR / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
logger.info(
    "gtex, %d variant-tissue pairs: NMD-triggering mean %.2f, from %.2f to %.2f between the %d tissues without"
    " the %d with too few samples",
    len(gtex_df),
    mean_nmd_triggering(gtex_df),
    min(tissue_means.values()),
    max(tissue_means.values()),
    len(tissue_means),
    len(EXCLUDED_TISSUES),
)

# %% [markdown]
# The scatters show the measured NMD efficiency of each measurement over the prediction of the model, per cohort. The
# gray line marks equal values.

# %%
scatter_df = facet_frame(
    {
        f"{cohort_names[cohort]}, {len(cohort_df):,} {MEASUREMENTS[cohort]}\n"
        f"Spearman {validation[cohort]['model']:.2f}": cohort_df
        for cohort, cohort_df in cohorts.items()
    }
)
plot = (
    pn.ggplot(scatter_df, pn.aes("prediction", "NMD_efficiency"))
    + pn.geom_point(color=BLUE, alpha=0.3, size=1.5, stroke=0)
    + pn.geom_abline(color=AXIS)
    + pn.facet_wrap("panel", scales="free")
    + pn.labs(x="Predicted NMD efficiency", y="Measured NMD efficiency")
)
save_figure(plot, OUT_DIR, "measured_vs_predicted", 8, 7.5)

# %% [markdown]
# Per cohort, the dots show the mean measured over the mean predicted NMD efficiency in each quantile bin of the
# prediction, with the 95% CI of the measured mean. The gray line marks equal values, and the blue line is the
# least-squares fit of measured on predicted over all measurements, the linear recalibration that the main README
# suggests.

# %%
calibrations = {}
for cohort, cohort_df in cohorts.items():
    bins_df, slope, intercept = calibration(cohort_df)
    title = (
        f"{cohort_names[cohort]}, {len(cohort_df):,} {MEASUREMENTS[cohort]}\n"
        f"Fit: slope {slope:.2f}, intercept {intercept:.2f}"
    )
    calibrations[title] = bins_df.with_columns(slope=slope, intercept=intercept)
plot = calibration_plot(facet_frame(calibrations)) + pn.facet_wrap("panel", scales="free")
save_figure(plot, OUT_DIR, "calibration", 8, 9)

# %% [markdown]
# The boxes show the measured NMD efficiency of the measurements by `nmd_escape`, per cohort. A box spans the quartiles,
# and its whiskers reach the furthest measurement within 1.5 times the interquartile range. The number below a box
# counts its measurements.

# %%
boxes_df = facet_frame(
    {
        f"{cohort_names[cohort]}\n{len(cohort_df):,} {MEASUREMENTS[cohort]}": cohort_df.with_columns(
            pl.col("nmd_escape").replace_strict(
                {escape: f"{escape}\n{count:,}" for escape, count in cohort_df["nmd_escape"].value_counts().iter_rows()}
            )
        )
        for cohort, cohort_df in cohorts.items()
    }
)
plot = (
    pn.ggplot(boxes_df, pn.aes("nmd_escape", "NMD_efficiency"))
    + pn.geom_boxplot(color=BLUE, fill=BLUE, alpha=0.1, width=0.5, outlier_size=1, outlier_alpha=0.3, outlier_stroke=0)
    + pn.facet_wrap("panel", nrow=1, scales="free")
    + pn.labs(x="nmd_escape", y="Measured NMD efficiency")
    + pn.theme(panel_grid_major_x=pn.element_blank())
)
save_figure(plot, OUT_DIR, "efficiency_by_nmd_escape", 9.5, 3.5)

# %% [markdown]
# The bars show the Spearman correlation of the model and of the canonical rules with the measured NMD efficiency, per
# cohort, with the 95% bootstrap CIs. These are the numbers of the table in the section "Model validation" of the main
# README. The bracket of a cohort gives the P value of the gain of the model over the rules.

# %%
bars_df = intervals({cohort: validation[cohort] for cohort in cohorts}, MODEL_AND_RULES)
# The bracket of each cohort spans the centers of its two bars, just right of their CIs
brackets_df = pl.DataFrame(
    {
        "position": -np.arange(len(cohorts), dtype=float),
        "rho": [max(validation[cohort]["model_ci"][1], validation[cohort]["rules_ci"][1]) + 0.03 for cohort in cohorts],
        "label": [p_value(validation[cohort]["gain_p"]) for cohort in cohorts],
    }
)
top, bottom = (offset for *_, offset in MODEL_AND_RULES)
# The line of each bracket, and a short tick at each end that points to the bars
bracket_lines_df = pl.concat(
    brackets_df.select(x=pl.col("position") + start, xend=pl.col("position") + end, y=pl.col("rho") - tick, yend="rho")
    for start, end, tick in [(top, bottom, 0.0), (top, top, 0.015), (bottom, bottom, 0.015)]
)
plot = (
    pn.ggplot(bars_df, pn.aes("position", "value"))
    + pn.geom_col(pn.aes(fill="color"), width=0.27)
    + pn.geom_errorbar(pn.aes(ymin="low", ymax="high"), width=0.12, color=SECONDARY_INK)
    + pn.geom_segment(pn.aes("x", "y", xend="xend", yend="yend"), data=bracket_lines_df, color=SECONDARY_INK)
    + pn.geom_text(
        pn.aes(y="rho", label="label"), data=brackets_df, ha="left", nudge_y=0.012, color=SECONDARY_INK, size=9
    )
    + pn.scale_fill_identity(
        guide="legend",
        breaks=[color for _, _, color, _ in MODEL_AND_RULES],
        labels=[label for _, label, _, _ in MODEL_AND_RULES],
    )
    + pn.scale_x_continuous(
        breaks=-np.arange(len(cohorts)),
        labels=[f"{cohort_names[cohort]}\n{validation[cohort]['n']:,} {MEASUREMENTS[cohort]}" for cohort in cohorts],
        expand=(0, 0.4),
    )
    + pn.scale_y_continuous(limits=(0, 1), breaks=np.linspace(0, 1, 6), expand=(0, 0))
    + pn.coord_flip()
    + pn.labs(x="", y="Spearman ρ with measured NMD efficiency (95% CI)")
    + pn.theme(legend_title=pn.element_blank(), panel_grid_major_y=pn.element_blank())
)
save_figure(plot, OUT_DIR, "spearman_per_cohort", 7, 3.5)

# %% [markdown]
# Per cohort, the bars split the R2 of the model into three terms: the squared Pearson correlation r of the prediction
# with the measured NMD efficiency, minus the squared difference of their means, minus the squared difference of the
# standard deviation of the prediction and r times that of the measured NMD efficiency. Each term is a share of the
# variance of the measured NMD efficiency, and R2 is the r2_score of scikit-learn. The linear recalibration that the main
# README suggests removes both error terms, so that R2 becomes the squared Pearson correlation.

# %%
r2_panels = {}
for cohort, cohort_df in cohorts.items():
    values = [sign * validation[cohort]["r2_decomposition"][term] for term, sign, _ in R2_TERMS]
    r2_panels[f"{cohort_names[cohort]}\n{len(cohort_df):,} {MEASUREMENTS[cohort]}"] = pl.DataFrame(
        {"term": [name for _, _, name in R2_TERMS], "value": values, "label": [f"{value:+.2f}" for value in values]}
    )
positive = pl.col("value") >= 0
# coord_flip draws the first term at the bottom, so the terms go in reverse order.
r2_df = (
    facet_frame(r2_panels)
    .cast({"term": pl.Enum([name for _, _, name in reversed(R2_TERMS)])})
    .with_columns(
        label_position=pl.col("value") + pl.when(positive).then(0.02).otherwise(-0.02),
        ha=pl.when(positive).then(pl.lit("left")).otherwise(pl.lit("right")),
    )
)
plot = (
    pn.ggplot(r2_df, pn.aes("term", "value"))
    + pn.geom_hline(yintercept=0, color=AXIS)
    + pn.geom_col(fill=BLUE, width=0.6)
    + pn.geom_text(pn.aes(y="label_position", label="label", ha="ha"), color=SECONDARY_INK, size=9)
    + pn.scale_y_continuous(breaks=[-0.4, 0, 0.4], expand=(0, 0.3))
    + pn.coord_flip()
    + pn.facet_wrap("panel", nrow=1)
    + pn.labs(x="", y="Share of the variance of the measured NMD efficiency")
    + pn.theme(panel_grid_major_y=pn.element_blank())
)
save_figure(plot, OUT_DIR, "r2_decomposition", 10, 3)

# %% [markdown]
# The dots show the mean measured NMD efficiency of the GTEx variant-tissue pairs with `nmd_escape` False, per tissue
# without the 5 tissues with too few samples, with the 95% CI of the mean. The gray line marks the mean over all
# variant-tissue pairs with `nmd_escape` False.

# %%
pairs_result = validation["gtex_variant_tissue_pairs"]
# The tissues from the lowest mean to the highest, which coord_flip draws from the bottom to the top
tissues = list(pairs_result["mean_nmd_triggering_per_tissue"])
nmd_triggering_df = gtex_df.filter(~pl.col("nmd_escape"), pl.col("tissue").is_in(tissues)).cast(
    {"tissue": pl.Enum(tissues)}
)
plot = (
    pn.ggplot(nmd_triggering_df, pn.aes("tissue", "NMD_efficiency"))
    + pn.geom_hline(yintercept=pairs_result["mean_nmd_triggering"], color=AXIS, linetype="dashed")
    + pn.annotate(
        "text",
        x=tissues[-1],
        y=pairs_result["mean_nmd_triggering"],
        label=f"All tissues {pairs_result['mean_nmd_triggering']:.2f}",
        ha="left",
        nudge_y=0.01,
        color=SECONDARY_INK,
        size=9,
    )
    + pn.stat_summary(fun_data="mean_cl_normal", color=BLUE)
    + pn.coord_flip()
    + pn.labs(
        title=f"GTEx, NMD-triggering variant-tissue pairs of {len(tissues)} tissues",
        x="",
        y="Mean measured NMD efficiency, 95% CI",
    )
    + pn.theme(panel_grid_major_y=pn.element_blank())
)
save_figure(plot, OUT_DIR, "gtex_tissues", 8, 9)

# %% [markdown]
# The dots show the Spearman correlation of the model and of the canonical rules with the measured NMD efficiency, per
# GTEx tissue without the 5 tissues with too few samples. The dashed line marks equal values. The labels name the 3
# tissues with the highest and the 3 with the lowest correlation of the model.

# %%
# The label of the lowest tissue hangs below its point, so that it clears the label of the next one.
tissue_labels_df = pl.concat(
    [
        tissue_spearman_df.tail(3).with_columns(x=pl.col("rules") - 0.008, ha=pl.lit("right"), va=pl.lit("center")),
        tissue_spearman_df.head(3).with_columns(
            x=pl.col("rules") + 0.008, ha=pl.lit("left"), va=pl.Series(["top", "center", "center"])
        ),
    ]
).with_columns(label=pl.format("{}, {} variants", "tissue", "n"))
plot = (
    pn.ggplot(tissue_spearman_df, pn.aes("rules", "model"))
    + pn.geom_abline(color=AXIS, linetype="dashed")
    + pn.geom_point(color=BLUE, size=2.5)
    # Both axes span the same range, so that the square panel gives them the same scale
    + pn.geom_blank(pn.aes(x="model", y="rules"))
    + pn.geom_text(pn.aes(x="x", label="label", ha="ha", va="va"), data=tissue_labels_df, color=SECONDARY_INK, size=9)
    + pn.labs(
        title=f"GTEx, Spearman ρ with the measured NMD efficiency in {len(tissue_spearman_df)} tissues",
        x="Spearman ρ, canonical rules",
        y="Spearman ρ, model",
    )
    + pn.theme(aspect_ratio=1)
)
save_figure(plot, OUT_DIR, "gtex_tissue_spearman", 6.5, 6.5)

# %% [markdown]
# The dots show the Spearman correlations and the gain of the model over the rules in MMRF-TARGET, for the
# variant-sample pairs whose transcript is in the TCGA training variant-sample pairs of TRAIN_DIR and for the other
# variant-sample pairs, with 95% bootstrap CIs. The numbers are `mmrf_target_by_transcript` in validation.json.

# %%
transcript_split = validation["mmrf_target_by_transcript"]
gain_panel = "Gain of the model over the rules"
transcript_intervals_df = facet_frame(
    {
        "Spearman correlation": intervals(transcript_split, MODEL_AND_RULES),
        gain_panel: intervals(transcript_split, GAIN),
    }
)
zero_df = pl.DataFrame({"panel": [gain_panel], "x": [0.0]}).cast({"panel": transcript_intervals_df["panel"].dtype})
plot = (
    pn.ggplot(transcript_intervals_df, pn.aes(y="position", color="color"))
    + pn.geom_segment(pn.aes(x="low", xend="high", yend="position"), size=1)
    + pn.geom_point(pn.aes(x="value"), size=2.5)
    + pn.geom_text(pn.aes(x="high", label="label"), ha="left", nudge_x=0.01, color=SECONDARY_INK, size=9)
    + pn.geom_vline(pn.aes(xintercept="x"), data=zero_df, color=AXIS)
    + pn.scale_color_identity()
    + pn.scale_y_continuous(
        breaks=-np.arange(len(transcript_split)),
        labels=[
            f"{TRANSCRIPT_GROUPS[group]}\n{group_correlations['n']:,} variant-sample pairs"
            for group, group_correlations in transcript_split.items()
        ],
        expand=(0, 0.4),
    )
    + pn.facet_wrap("panel", scales="free_x")
    + pn.labs(title="MMRF-TARGET by transcript, 95% CIs", x="")
    + pn.theme(axis_title_y=pn.element_blank(), panel_grid_major_y=pn.element_blank())
)
save_figure(plot, OUT_DIR, "mmrf_target_by_transcript", 9, 3)
