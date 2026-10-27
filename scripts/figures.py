"""
The plotnine theme, the calibration plot and save_figure of train_model.py and validate_model.py. calibration computes
the statistics of the calibration plot, and save_figure saves each figure as PNG and PDF.

The colors come from a categorical palette that is checked for color vision deficiency: BLUE marks the model and every
single series, ORANGE marks the rules, and MUTED_INK marks a series that gives context.
"""

import sys
from pathlib import Path

import matplotlib
import numpy as np
import polars as pl
from plotnine import (
    aes,
    element_blank,
    element_line,
    element_rect,
    element_text,
    geom_abline,
    geom_blank,
    geom_linerange,
    geom_point,
    ggplot,
    ggsave,
    labs,
    theme,
    theme_minimal,
    theme_set,
)
from scipy.stats import t

BLUE = "#2a78d6"
ORANGE = "#eb6834"
INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED_INK = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#ffffff"
DPI = 600
# The calibration plots bin the rows into this many quantile bins of the prediction
CALIBRATION_BINS = 10

theme_set(
    theme_minimal(base_size=9)
    + theme(
        plot_background=element_rect(fill=SURFACE, color=SURFACE),
        plot_title=element_text(size=10, color=INK, ha="left"),
        plot_caption=element_text(color=SECONDARY_INK, ha="left"),
        strip_text=element_text(size=9, color=INK, ha="left", va="top", ma="left"),
        axis_title=element_text(color=SECONDARY_INK),
        axis_text=element_text(size=9, color=MUTED_INK),
        axis_text_y=element_text(ma="right"),
        axis_line=element_line(color=AXIS, size=0.75),
        axis_ticks=element_line(color=AXIS, size=0.75),
        axis_ticks_minor=element_blank(),
        panel_grid_major=element_line(color=GRID, size=0.75),
        panel_grid_minor=element_blank(),
        panel_spacing_y=0.04,
    )
)
# Embed TrueType fonts, so that the text in the PDF stays text
matplotlib.rcParams["pdf.fonttype"] = 42


def in_notebook() -> bool:
    """Whether the code runs in a Jupyter kernel."""
    return "ipykernel" in sys.modules


def calibration(rows: pl.DataFrame) -> tuple[pl.DataFrame, float, float]:
    """
    One row per bin of CALIBRATION_BINS quantile bins of the column prediction of rows, from the lowest bin to the
    highest. A row has the mean prediction, the mean NMD_efficiency, and the 95% CI of the mean NMD_efficiency from the
    t distribution as low and high. Also the slope and the intercept of the least-squares line of NMD_efficiency on prediction over all rows.
    """
    bins = (
        rows.group_by(bin=pl.col("prediction").qcut(CALIBRATION_BINS, allow_duplicates=True))
        .agg(
            pl.col("prediction", "NMD_efficiency").mean(),
            sem=pl.col("NMD_efficiency").std() / pl.len().sqrt(),
            n=pl.len(),
        )
        .sort("prediction")
    )
    half_width = pl.Series(t.ppf(0.975, bins["n"].to_numpy() - 1)) * bins["sem"]
    measured = pl.col("NMD_efficiency")
    bins = bins.select("prediction", measured, low=measured - half_width, high=measured + half_width)
    slope, intercept = np.polyfit(rows["prediction"].to_numpy(), rows["NMD_efficiency"].to_numpy(), 1)
    return bins, float(slope), float(intercept)


def calibration_plot(bins: pl.DataFrame) -> ggplot:
    """
    The calibration plot of bins, the bins of calibration with its slope and intercept as columns. With facets, bins
    has the bins and the fit of each panel. Both axes of a panel span the same range, so that the square panel gives
    them the same scale.
    """
    return (
        ggplot(bins, aes("prediction", "NMD_efficiency"))
        + geom_abline(color=AXIS)
        # One line per fit
        + geom_abline(
            aes(slope="slope", intercept="intercept"),
            data=bins.unique(["slope", "intercept"]),
            color=BLUE,
        )
        + geom_linerange(aes(ymin="low", ymax="high"), color=BLUE)
        + geom_point(color=BLUE, size=2.5)
        + geom_blank(aes(x="low", y="prediction"))
        + geom_blank(aes(x="high", y="prediction"))
        + labs(
            x="Mean predicted NMD efficiency",
            y="Mean measured NMD efficiency",
            caption=f"Dots: mean in each of {CALIBRATION_BINS} quantile bins of the prediction,\n"
            "with the 95% CI of the measured mean.\n"
            "Blue line: least-squares fit over all rows. Gray line: measured = predicted.",
        )
        + theme(aspect_ratio=1)
    )


def save_figure(plot: ggplot, out_dir: Path, name: str, width: float, height: float) -> None:
    """Save plot as out_dir/figures/name.png and name.pdf, width x height inches at DPI. In a notebook, show it."""
    plot = plot + theme(figure_size=(width, height))
    figures = out_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        ggsave(plot, figures / f"{name}.{suffix}", dpi=DPI, verbose=False)
    if in_notebook():
        plot.show()
