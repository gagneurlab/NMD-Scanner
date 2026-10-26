"""
The figure style of train_model.py and validate_model.py, and save_figure, which saves each figure as PNG and PDF.

The colors come from a categorical palette that is checked for color vision deficiency: BLUE marks the model and every
single series, ORANGE marks the rules, and MUTED_INK marks a series that gives context.
"""

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.stats import sem, t

BLUE = "#2a78d6"
ORANGE = "#eb6834"
INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED_INK = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#ffffff"
DPI = 600
# A notebook shows each figure at this many pixels per inch
NOTEBOOK_PIXELS_PER_INCH = 100
# The calibration plots bin the rows into this many quantile bins of the prediction
CALIBRATION_BINS = 10
# The dots of every figure, and the whiskers or lines that join them
MARKER_STYLE = {"markersize": 6, "markeredgecolor": SURFACE, "markeredgewidth": 1, "linewidth": 1.5}

matplotlib.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titlelocation": "left",
        "axes.titlecolor": INK,
        "axes.labelcolor": SECONDARY_INK,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.75,
        "xtick.color": AXIS,
        "ytick.color": AXIS,
        "xtick.labelcolor": MUTED_INK,
        "ytick.labelcolor": MUTED_INK,
        "legend.frameon": False,
        "legend.labelcolor": SECONDARY_INK,
        # Embed TrueType fonts, so that the text in the PDF stays text
        "pdf.fonttype": 42,
    }
)


def in_notebook() -> bool:
    """Whether the code runs in a Jupyter kernel."""
    return "ipykernel" in sys.modules


def new_figure(width: float, height: float, nrows: int = 1, ncols: int = 1) -> tuple[Figure, np.ndarray]:
    """A figure of width x height inches with a grid of nrows x ncols axes, returned as a flat array."""
    fig = Figure(figsize=(width, height), layout="constrained")
    return fig, fig.subplots(nrows, ncols, squeeze=False).ravel()


def mean_ci_half_width(values: pd.Series) -> float:
    """The half width of the 95% CI of the mean of values, from the t distribution."""
    return float(t.ppf(0.975, len(values) - 1) * sem(values))


def plot_intervals(
    ax: Axes,
    values: np.ndarray,
    lows: np.ndarray,
    highs: np.ndarray,
    positions: np.ndarray,
    color: str,
) -> None:
    """A dot at each value on the x axis, at its position on the y axis, with a whisker from its low to its high."""
    ax.errorbar(values, positions, xerr=[values - lows, highs - values], fmt="o", color=color, **MARKER_STYLE)


def plot_model_and_rules(ax: Axes, results: list[dict], labels: list[str]) -> None:
    """
    For each result, the Spearman correlation of the model and of the rules with their 95% CIs, in rows from top to
    bottom. A result has the keys model, model_ci, rules and rules_ci. The first row has direct labels.
    """
    positions = -np.arange(len(results), dtype=float)
    for offset, score, label, color in ((0.15, "model", "Model", BLUE), (-0.15, "rules", "Canonical rules", ORANGE)):
        values = np.array([result[score] for result in results])
        cis = np.array([result[f"{score}_ci"] for result in results])
        plot_intervals(ax, values, cis[:, 0], cis[:, 1], positions + offset, color)
        ax.annotate(
            label,
            (cis[0, 1], positions[0] + offset),
            xytext=(4, 0),
            textcoords="offset points",
            va="center",
            color=SECONDARY_INK,
        )
    ax.set_yticks(positions, labels)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(right=ax.get_xlim()[1] + 0.08)


def plot_measured_vs_predicted(ax: Axes, prediction: np.ndarray, measured: np.ndarray, title: str) -> None:
    """A scatter of the measured over the predicted NMD efficiency, with a line where both are equal."""
    ax.scatter(prediction, measured, s=36, color=BLUE, alpha=0.3, edgecolors=SURFACE, linewidths=0.75)
    ax.axline((0, 0), slope=1, color=AXIS, linewidth=1)
    ax.set_title(title)
    ax.set_xlabel("Predicted NMD efficiency")
    ax.set_ylabel("Measured NMD efficiency")


def plot_calibration(ax: Axes, prediction: np.ndarray, measured: np.ndarray, title: str) -> None:
    """
    A calibration plot: per bin of CALIBRATION_BINS quantile bins of the prediction, the mean measured over the mean
    predicted NMD efficiency, with the 95% CI of the measured mean from the t distribution. Add the least-squares line
    of measured on predicted over all rows, with its slope and intercept in the title, and the line where both are
    equal. Both axes have the same scale.
    """
    rows = pd.DataFrame({"prediction": prediction, "measured": measured})
    bins = rows.groupby(pd.qcut(rows["prediction"], CALIBRATION_BINS, labels=False, duplicates="drop"))
    means = bins.mean()
    half_width = bins["measured"].agg(mean_ci_half_width)
    slope, intercept = np.polyfit(prediction, measured, 1)
    ax.axline((0, 0), slope=1, color=AXIS, linewidth=1, label="Measured = predicted")
    ax.axline((0, intercept), slope=slope, color=BLUE, linewidth=1.5, label="Least-squares fit over all rows")
    ax.errorbar(
        means["prediction"],
        means["measured"],
        yerr=half_width,
        fmt="o",
        color=BLUE,
        label=f"Mean in each of {CALIBRATION_BINS} quantile bins of the prediction, 95% CI",
        **MARKER_STYLE,
    )
    low = min(means["prediction"].min(), (means["measured"] - half_width).min())
    high = max(means["prediction"].max(), (means["measured"] + half_width).max())
    margin = 0.05 * (high - low)
    ax.set_xlim(low - margin, high + margin)
    ax.set_ylim(low - margin, high + margin)
    ax.set_aspect("equal")
    ax.set_title(f"{title}\nFit: slope {slope:.2f}, intercept {intercept:.2f}")
    ax.set_xlabel("Mean predicted NMD efficiency")
    ax.set_ylabel("Mean measured NMD efficiency")


def save_figure(fig: Figure, out_dir: Path, name: str) -> None:
    """Save fig as out_dir/figures/name.png and name.pdf at DPI. In a notebook, show the PNG."""
    figures = out_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(figures / f"{name}.{suffix}", dpi=DPI)
    if in_notebook():
        from IPython.display import Image, display

        display(
            Image(filename=str(figures / f"{name}.png"), width=round(fig.get_figwidth() * NOTEBOOK_PIXELS_PER_INCH))
        )
