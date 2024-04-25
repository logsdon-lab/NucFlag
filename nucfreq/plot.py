import numpy as np
import polars as pl
import portion as pt

import matplotlib.pyplot as plt

from .constants import PLOT_WIDTH, PLOT_HEIGHT, PLOT_YLIM
from .misassembly import Misassembly

from typing import Any


def plot_coverage(
    df: pl.DataFrame,
    misassemblies: dict[Misassembly, set[pt.Interval]],
    contig_name: str,
) -> tuple[plt.Figure, Any]:
    fig, ax = plt.subplots(figsize=(PLOT_WIDTH, PLOT_HEIGHT))

    (_,) = ax.plot(
        df["position"],
        df["first"],
        "o",
        color="black",
        markeredgewidth=0.0,
        markersize=2,
        label="Most Frequent Base",
    )
    (_,) = ax.plot(
        df["position"],
        df["second"],
        "o",
        color="red",
        markeredgewidth=0.0,
        markersize=2,
        label="Second Most Frequent Base",
    )

    # Add misassembly rect patches to highlight region.
    for misasm, regions in misassemblies.items():
        color = misasm.as_color()
        for region in regions:
            plt.axvspan(
                region.lower, region.upper, color=color, alpha=0.4, label=misasm
            )

    # Add legend. Deduplicate multiple labels.
    # https://stackoverflow.com/a/36189073
    handles, labels = plt.gca().get_legend_handles_labels()
    labels, ids = np.unique(labels, return_index=True)
    handles = [handles[i] for i in ids]
    plt.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.2),
        ncols=len(labels),
        borderaxespad=0,
        fancybox=True,
    )

    maxval = df["position"].max()
    minval = df["position"].min()
    subval = 0

    title = "{}:{}-{}\n".format(contig_name, minval, maxval)
    ax.set_title(title, fontweight="bold")

    if maxval < 1_000_000:
        xlabels = [format((label - subval), ",.0f") for label in ax.get_xticks()]
        lab = "bp"
    elif maxval < 10_000_000:
        xlabels = [format((label - subval) / 1000, ",.1f") for label in ax.get_xticks()]
        lab = "kbp"
    else:
        xlabels = [format((label - subval) / 1000, ",.1f") for label in ax.get_xticks()]
        lab = "kbp"

    ax.set_ylim(0, PLOT_YLIM)
    ax.set_xlabel("Assembly position ({})".format(lab), fontweight="bold")
    ax.set_ylabel("Sequence read depth", fontweight="bold")
    ax.set_xticklabels(xlabels)

    # Hide the right and top spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    fig.tight_layout()

    return fig, ax
