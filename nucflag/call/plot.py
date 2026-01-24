import random
import logging
import warnings
import matplotlib
import matplotlib.axes
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.patches as ptch

from typing import Any
from collections import OrderedDict
from intervaltree import Interval  # type: ignore[import-untyped]
from matplotlib.collections import PatchCollection
from functools import lru_cache

from .region import ActionOpt, Region
from ..common import minimalize_ax

PLOT_FONT_SIZE = 16
PLOT_HEIGHT = 5
PLOT_WIDTH = 24
PLOT_HEIGHT_SCALE_FACTOR = 2

# No margins.
matplotlib.use("agg")
warnings.filterwarnings("ignore")

# Set text size
matplotlib.rcParams.update({"font.size": PLOT_FONT_SIZE})
plt.rcParams["axes.xmargin"] = 0


@lru_cache
def get_random_color(name: str) -> str:
    # Seed the RNG with the name
    random.seed(name)
    # Then generate rgb and hexcode
    rgb = tuple(random.randint(0, 255) for _ in range(3))
    assert len(rgb) == 3
    # https://stackoverflow.com/a/3380739
    return "#%02x%02x%02x" % rgb


def plot_coverage(
    itv: Interval,
    df_pileup: pl.DataFrame,
    tracks: OrderedDict[str, set[Region]],
    ovl_tracks: OrderedDict[str, set[Region]],
    plot_ylim: float | int = 100,
) -> tuple[plt.Figure, Any]:
    """
    :param itv: Interval to plot
    :type itv: Interval
    :param df_pileup: Pileup dataframe
    :type df_pileup: pl.DataFrame
    :param tracks: Tracks to plot
    :type tracks: OrderedDict[str, set[Region]]
    :param ovl_tracks: Tracks to plot on top of coverage plot.
    :type ovl_tracks: OrderedDict[str, set[Region]]
    :param plot_ylim: Plot y-limit.
    :type plot_ylim: float | int
    :return: Figure and its axes.
    :rtype: tuple[Figure, Any]
    """

    color: str | None
    subplot_patches: dict[str, list[ptch.Rectangle]] = {}
    number_of_overlap_beds = len(tracks.keys()) if tracks else 0
    if tracks:
        fig, axs = plt.subplots(
            # | ---------- |
            # | \/\__/\||/ |
            # | bed_legend |
            # | cov_legend |
            # nrows (n for bed, 1 for the coverage plot, n for bed legends, and 1 for coverage plot legend)
            len(tracks.keys()) + 1 + len(tracks.keys()) + 1,
            # ncols
            1,
            figsize=(
                PLOT_WIDTH,
                PLOT_HEIGHT + (PLOT_HEIGHT_SCALE_FACTOR * number_of_overlap_beds),
            ),
            sharex=True,
            gridspec_kw={
                "height_ratios": [
                    *[0.3 for _ in range(number_of_overlap_beds)],
                    3,
                    *[0.15 for _ in range(number_of_overlap_beds)],
                    0.15,
                ]
            },
            layout="constrained",
        )

        # Last axis with coverage plt.
        ax: matplotlib.axes.Axes = axs[number_of_overlap_beds]
        # Add bed regions.
        for i, (name, regions) in enumerate(tracks.items()):
            # Make axis as minimal as possible.
            bed_axs: matplotlib.axes.Axes = axs[i]
            minimalize_ax(bed_axs, remove_ticks=True)

            bed_axs.set_ylabel(
                name,
                rotation="horizontal",
                ma="right",
                ha="right",
                va="center",
            )

            # Map uniq types to new color if none given.
            patches: list[ptch.Rectangle] = []
            ylim = bed_axs.get_ylim()[1]
            for row in regions:
                # Skip rows not within bounds of df.
                if not itv.overlaps(row.region):
                    continue

                if not row.action or (row.action and row.action.opt != ActionOpt.PLOT):
                    continue

                region_st = max(itv.begin, row.region.begin)
                region_end = min(itv.end, row.region.end)
                width = region_end - region_st
                # Use color provided. Default to random generated ones otherwise.
                if row.action.desc:
                    color = row.action.desc
                elif row.desc:
                    color = get_random_color(row.desc)
                else:
                    raise ValueError(f"Region {row} has no description.")

                rect = ptch.Rectangle(
                    (region_st, 0),
                    width,
                    ylim,
                    linewidth=1,
                    edgecolor=None,
                    facecolor=color,
                    alpha=0.75,
                    label=row.desc,
                )
                patches.append(rect)

            bed_axs.add_collection(PatchCollection(patches, match_original=True))
            subplot_patches[name] = patches
    else:
        fig, axs = plt.subplots(
            2,
            1,
            figsize=(PLOT_WIDTH, PLOT_HEIGHT),
            layout="constrained",
            gridspec_kw={"height_ratios": [3, 0.15]},
        )
        ax = axs[0]

    for typ, markercolor, color, label in (
        ("insertion", "#800080", "#800080", "Insertions"),
        ("deletion", "#f8f4ff", "#000000", "Deletions"),
        ("mismatch", "#ff0000", "#ff0000", "Mismatches"),
        ("cov", "#000000", "#000000", "Coverage"),
    ):
        ax.plot(
            df_pileup["pos"],
            df_pileup[typ],
            marker="o",
            linestyle="None",
            color=color,
            markersize=2,
            markeredgewidth=0.2,
            markerfacecolor=markercolor,
            label=label,
        )

    for track_name, regions in ovl_tracks.items():
        for region in regions:
            if region.name == "correct" and track_name == "Calls":
                continue
            if region.action:
                color = region.action.desc
            elif region.desc:
                color = get_random_color(region.desc)
            else:
                raise ValueError(f"Region {region} has no description.")

            if not region.region.overlaps(itv):
                continue

            region_st = max(itv.begin, region.region.begin)
            region_end = min(itv.end, region.region.end)
            ax.axvspan(
                region_st,
                region_end,
                color=color,
                alpha=0.4,
                label=region.desc,
            )

    # Add legend for coverage plot in separate axis. Deduplicate multiple labels.
    # https://stackoverflow.com/a/36189073
    handles, labels = ax.get_legend_handles_labels()
    labels_handles = dict(zip(labels, handles))
    legend_cov_ax: matplotlib.axes.Axes = axs[-1]
    minimalize_ax(legend_cov_ax, remove_ticks=True)
    legend_cov_ax.legend(
        labels_handles.values(),
        labels_handles.keys(),
        loc="center",
        alignment="left",
        ncols=len(labels),
        borderaxespad=0,
        handlelength=1.0,
        handleheight=1.0,
        fancybox=False,
        frameon=False,
        prop={"size": 12},
        title="Pileup",
    )
    # Add legends for each overlapped bedfile.
    for i, (name, sp_patches) in enumerate(subplot_patches.items()):
        # Remove plot elements from legend ax.
        # Offset by 1 for coverage plot.
        legend_ax: matplotlib.axes.Axes = axs[i + number_of_overlap_beds + 1]
        minimalize_ax(legend_ax, remove_ticks=True)

        # Filter rectangle patches.
        sp_patch_labels = set()
        sp_filtered_patches: list[ptch.Rectangle] = []
        for patch in sp_patches:
            patch_lbl = patch.get_label()
            if patch_lbl in sp_patch_labels:
                continue
            sp_filtered_patches.append(patch)
            sp_patch_labels.add(patch_lbl)

        # Sort patches by label.
        # Special exception for MAPQ since it's a range and cannot sort lexographically. ex. 5-10 and 50-55
        if name == "MAPQ":
            sp_filtered_patches.sort(key=lambda p: int(p.get_label().split("-")[0]))
        elif name == "Self-identity":
            # 99.25%-99.50%
            sp_filtered_patches.sort(
                key=lambda p: float(p.get_label().split("-")[0][0:-1])
            )
        else:
            sp_filtered_patches.sort(key=lambda p: str(p.get_label()))

        legend_ax.legend(
            handles=sp_filtered_patches,
            loc="center",
            alignment="left",
            title=name,
            # At least 1-15 columns in legend.
            ncols=min(max(len(sp_filtered_patches), 1), 15),
            handlelength=1.0,
            handleheight=1.0,
            borderaxespad=0,
            fancybox=False,
            frameon=False,
            prop={"size": 12},
        )

    fig.suptitle("{}:{}-{}\n".format(itv.data, itv.begin, itv.end), fontweight="bold")

    if itv.end < 1_000_000:
        xlabels = [format(label, ",.0f") for label in ax.get_xticks()]
        lab = "bp"
    elif itv.end < 10_000_000:
        xlabels = [format(label / 1000, ",.1f") for label in ax.get_xticks()]
        lab = "kbp"
    else:
        xlabels = [format(label / 1000, ",.1f") for label in ax.get_xticks()]
        lab = "kbp"

    mean_cov = df_pileup["cov"].mean()
    if isinstance(plot_ylim, float):
        assert isinstance(mean_cov, (float, int)), "Invalid mean coverage type."
        plot_ylim = mean_cov * plot_ylim
    elif isinstance(plot_ylim, int):
        plot_ylim = plot_ylim
    else:
        logging.error(
            f"Invalid {plot_ylim} of type {type(plot_ylim)}. Defaulting to 100."
        )
        plot_ylim = 100

    ax.set_ylim(0, plot_ylim)
    ax.set_xlabel("Assembly position ({})".format(lab), fontweight="bold")
    ax.set_ylabel(
        "Sequence\nread\ndepth",
        fontweight="bold",
        rotation="horizontal",
        ma="center",
        ha="right",
        va="center",
    )
    ax.set_xticklabels(xlabels)

    # Hide the right and top spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")

    return fig, axs
