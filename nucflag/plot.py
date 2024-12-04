import warnings
from typing import Any, DefaultDict

import matplotlib
import matplotlib.axes
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
from intervaltree import Interval, IntervalTree

from .misassembly import Misassembly
from .region import ActionOpt, Region

PLOT_FONT_SIZE = 16
PLOT_HEIGHT = 5
PLOT_WIDTH = 24
PLOT_YLIM = 100
PLOT_HEIGHT_SCALE_FACTOR = 2

# No margins.
matplotlib.use("agg")
warnings.filterwarnings("ignore")

# Set text size
matplotlib.rcParams.update({"font.size": PLOT_FONT_SIZE})
plt.rcParams["axes.xmargin"] = 0


def minimalize_ax(ax: matplotlib.axes.Axes) -> None:
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_frame_on(False)


def plot_coverage(
    df: pl.DataFrame,
    misassemblies: dict[Misassembly, IntervalTree],
    contig_name: str,
    overlay_regions: DefaultDict[int, set[Region]] | None,
) -> tuple[plt.Figure, Any]:
    region_bounds = Interval(df["position"].min(), df["position"].max())

    subplot_handles = []
    subplot_labels = []

    number_of_overlap_beds = len(overlay_regions.keys()) if overlay_regions else 0
    if overlay_regions:
        fig, axs = plt.subplots(
            # | ---------- |
            # | \/\__/\||/ |
            # | cov_legend |
            # | bed_legend |
            # nrows (n for bed, 1 for the coverage plot, 1 for coverage plot legend, and n for bed legends)
            len(overlay_regions.keys()) + 1 + 1 + len(overlay_regions.keys()),
            # ncols
            1,
            figsize=(
                PLOT_WIDTH,
                PLOT_HEIGHT + (PLOT_HEIGHT_SCALE_FACTOR * number_of_overlap_beds),
            ),
            sharex=True,
            gridspec_kw={
                "height_ratios": [
                    *[0.5 for _ in range(number_of_overlap_beds)],
                    3,
                    0.2,
                    *[0.3 for _ in range(number_of_overlap_beds)],
                ]
            },
            layout="constrained",
        )
        plt.subplots_adjust(hspace=0.6)

        # Last axis with coverage plt.
        ax: matplotlib.axes.Axes = axs[number_of_overlap_beds]
        # Add bed regions.
        for i, regions in overlay_regions.items():
            # Make axis as minimal as possible.
            bed_axs: matplotlib.axes.Axes = axs[i]
            minimalize_ax(bed_axs)

            # Map uniq types to new color if none given.
            uniq_types = sorted({r.desc for r in regions if r.desc})
            cmap = dict(
                zip(
                    uniq_types,
                    (
                        c
                        for c in plt.colormaps.get_cmap("hsv")(
                            np.linspace(0, 1, len(uniq_types))
                        )
                    ),
                )
            )
            for row in regions:
                # Skip rows not within bounds of df.
                if not region_bounds.overlaps(row.region):
                    continue

                if not row.action or (row.action and row.action.opt != ActionOpt.PLOT):
                    continue
                width = row.region.length()
                # Use color provided. Default to random generated ones otherwise.
                if row.action.desc:
                    color = row.action.desc
                elif row.desc:
                    color = cmap[row.desc]
                else:
                    raise ValueError(f"Region {row} has no description.")

                rect = ptch.Rectangle(
                    (row.region.begin, 0),
                    width,
                    0.5,
                    linewidth=1,
                    edgecolor=None,
                    facecolor=color,
                    alpha=0.75,
                    label=row.desc,
                )
                bed_axs.add_patch(rect)

            # Get legend elements.
            handles, labels = bed_axs.get_legend_handles_labels()
            labels, ids = np.unique(labels, return_index=True)
            handles = [handles[i] for i in ids]
            subplot_handles.append(handles)
            subplot_labels.append(labels)
    else:
        fig, axs = plt.subplots(
            2,
            1,
            figsize=(PLOT_WIDTH, PLOT_HEIGHT),
            layout="constrained",
            gridspec_kw={"height_ratios": [3, 0.2]},
        )
        ax = axs[0]

    (_,) = ax.plot(
        df["position"],
        df["second"],
        "o",
        color="red",
        markeredgewidth=0.0,
        markersize=2,
        label="Second Most Frequent Base",
    )
    (_,) = ax.plot(
        df["position"],
        df["first"],
        "o",
        color="black",
        markeredgewidth=0.0,
        markersize=2,
        label="Most Frequent Base",
    )
    # Add misassembly rect patches to highlight region.
    for misasm, misasm_regions in misassemblies.items():
        color = misasm.as_color()
        for misasm_region in misasm_regions:
            ax.axvspan(
                misasm_region.begin,
                misasm_region.end,
                color=color,
                alpha=0.4,
                label=misasm,
            )

    # Add legend for coverage plot in separate axis. Deduplicate multiple labels.
    # https://stackoverflow.com/a/36189073
    handles, labels = ax.get_legend_handles_labels()
    labels, ids = np.unique(labels, return_index=True)
    handles = [handles[i] for i in ids]
    legend_cov_ax: matplotlib.axes.Axes = axs[number_of_overlap_beds + 1]
    minimalize_ax(legend_cov_ax)
    legend_cov_ax.legend(
        handles,
        labels,
        loc="center",
        ncols=len(labels),
        borderaxespad=0,
        fancybox=True,
    )
    # Add legends for each overlapped bedfile.
    for i, (sp_handles, sp_labels) in enumerate(zip(subplot_handles, subplot_labels)):
        # Remove plot elements from legend ax.
        # Offset by 2 for coverage plot and its legend.
        legend_ax: matplotlib.axes.Axes = axs[i + number_of_overlap_beds + 2]
        minimalize_ax(legend_ax)

        legend_ax.legend(
            sp_handles,
            sp_labels,
            loc="center",
            # Must have at least one col.
            ncols=max(len(sp_handles) // 3, 1),
            borderaxespad=0,
            fancybox=True,
        )

    maxval = df["position"].max()
    minval = df["position"].min()
    subval = 0

    title = "{}:{}-{}\n".format(contig_name, minval, maxval)
    plt.suptitle(title, fontweight="bold")

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

    return fig, axs
