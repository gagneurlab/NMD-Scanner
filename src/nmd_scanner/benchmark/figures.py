"""
The plotnine theme, the calibration plot and save_figure of scripts/validate_model.py.py.
calibration computes the statistics of the calibration plot, and save_figure saves each figure as PNG and PDF. This
module needs the extra scripts of nmd-scanner.

The theme and the colors follow the style of the AbExp paper figures: the theme_cowplot of cowplot with its background
grid, at FONT_SIZE and LINE_SIZE. BLUE marks the model and the only series of a figure, and LIGHT_GRAY marks the rules.
GRAY marks a series that gives context, and the reference lines. Under simulated protanopia, deuteranopia and
tritanopia, BLUE stays apart from both grays. BLUE and GRAY have about the same lightness, so a figure tells them apart
by direct labels, not only by their color.
"""

from pathlib import Path

import matplotlib
import numpy as np
import plotnine as pn
import polars as pl
from scipy.stats import t

BLUE = "#55A2CF"
LIGHT_GRAY = "#BDBDBD"
GRAY = "#999999"
INK = "black"
GRID = "#D9D9D9"
STRIP = "#CCCCCC"
SURFACE = "#ffffff"
DPI = 600
# The base font size in points, and the width of the axis lines in points
FONT_SIZE = 12
LINE_SIZE = 1
# The calibration plots bin the measurements into this many quantile bins of the prediction
CALIBRATION_BINS = 10

# theme_cowplot sets the axis text, the legend text and the strip text at 12/14 of the base size, and the caption at
# 11/14 of it. background_grid adds the major grid lines.
pn.theme_set(
    pn.theme_classic(base_size=FONT_SIZE)
    + pn.theme(
        text=pn.element_text(color=INK),
        plot_background=pn.element_rect(fill=SURFACE, color=SURFACE),
        plot_title=pn.element_text(size=FONT_SIZE, ha="left"),
        plot_caption=pn.element_text(size=FONT_SIZE * 11 / 14, ha="left"),
        strip_background=pn.element_rect(fill=STRIP, color=STRIP),
        strip_text=pn.element_text(size=FONT_SIZE * 12 / 14),
        axis_title=pn.element_text(size=FONT_SIZE),
        axis_text=pn.element_text(size=FONT_SIZE * 12 / 14, color=INK),
        axis_text_y=pn.element_text(ma="right"),
        axis_line=pn.element_line(color=INK, size=LINE_SIZE),
        axis_ticks=pn.element_line(color=INK, size=LINE_SIZE),
        axis_ticks_length=FONT_SIZE / 4,
        axis_ticks_minor=pn.element_blank(),
        panel_grid_major=pn.element_line(color=GRID, size=0.7 * LINE_SIZE),
        panel_grid_minor=pn.element_blank(),
        panel_spacing_y=0.04,
        legend_title=pn.element_blank(),
        legend_text=pn.element_text(size=FONT_SIZE * 12 / 14),
        legend_key_size=FONT_SIZE,
    )
)
# Embed TrueType fonts, so that the text in the PDF stays text
matplotlib.rcParams["pdf.fonttype"] = 42


def calibration(cohort_df: pl.DataFrame) -> tuple[pl.DataFrame, float, float]:
    """
    One row per bin of CALIBRATION_BINS quantile bins of the column prediction of cohort_df, from the lowest bin to the
    highest. A bin has the mean prediction, the mean NMD_efficiency, and the 95% CI of the mean NMD_efficiency from the
    t distribution as low and high. Also the slope and the intercept of the least-squares line of NMD_efficiency on
    prediction over all measurements.
    """
    bins_df = (
        cohort_df.group_by(bin=pl.col("prediction").qcut(CALIBRATION_BINS, allow_duplicates=True))
        .agg(
            pl.col("prediction", "NMD_efficiency").mean(),
            sem=pl.col("NMD_efficiency").std() / pl.len().sqrt(),
            n=pl.len(),
        )
        .sort("prediction")
    )
    half_width = pl.Series(t.ppf(0.975, bins_df["n"].to_numpy() - 1)) * bins_df["sem"]
    measured = pl.col("NMD_efficiency")
    bins_df = bins_df.select("prediction", measured, low=measured - half_width, high=measured + half_width)
    slope, intercept = np.polyfit(cohort_df["prediction"].to_numpy(), cohort_df["NMD_efficiency"].to_numpy(), 1)
    return bins_df, float(slope), float(intercept)


def calibration_plot(bins_df: pl.DataFrame) -> pn.ggplot:
    """
    The calibration plot of bins_df, the bins of calibration with its slope and intercept as columns. With facets,
    bins_df has the bins and the fit of each panel. Both axes of a panel span the same range, so that the square panel gives
    them the same scale.
    """
    return (
        pn.ggplot(bins_df, pn.aes("prediction", "NMD_efficiency"))
        + pn.geom_abline(color=GRAY)
        # One line per fit
        + pn.geom_abline(
            pn.aes(slope="slope", intercept="intercept"),
            data=bins_df.unique(["slope", "intercept"]),
            color=BLUE,
        )
        + pn.geom_linerange(pn.aes(ymin="low", ymax="high"), color=BLUE)
        + pn.geom_point(color=BLUE, size=2.5)
        + pn.geom_blank(pn.aes(x="low", y="prediction"))
        + pn.geom_blank(pn.aes(x="high", y="prediction"))
        + pn.labs(
            x="Mean predicted NMD efficiency",
            y="Mean measured NMD efficiency",
            caption=f"Dots: mean in each of {CALIBRATION_BINS} quantile bins of the prediction,\n"
            "with the 95% CI of the measured mean.\n"
            "Blue line: least-squares fit over all measurements.\n"
            "Gray line: measured = predicted.",
        )
        + pn.theme(aspect_ratio=1)
    )


def save_figure(plot: pn.ggplot, out_dir: Path, name: str, width: float, height: float) -> pn.ggplot:
    """
    Save plot as out_dir/figures/name.png and name.pdf, width x height inches at DPI, and return it with this size. As
    the last expression of a notebook cell, the plot is the cell result, and Jupyter shows it. A plain python run
    ignores it.
    """
    plot = plot + pn.theme(figure_size=(width, height))
    figures = out_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        pn.ggsave(plot, figures / f"{name}.{suffix}", dpi=DPI, verbose=False)
    return plot
