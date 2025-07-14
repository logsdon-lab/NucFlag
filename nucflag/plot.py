import warnings
import matplotlib
import matplotlib.axes
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.patches as ptch

from typing import Any
from collections import OrderedDict
from intervaltree import Interval
from matplotlib.collections import PatchCollection

from .region import ActionOpt, Region

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


def minimalize_ax(ax: matplotlib.axes.Axes) -> None:
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_frame_on(False)


def plot_coverage(
    itv: Interval,
    df_pileup: pl.DataFrame,
    df_misasm: pl.DataFrame,
    overlay_regions: OrderedDict[str, set[Region]] | None,
) -> tuple[plt.Figure, Any]:
    subplot_patches: dict[str, list[ptch.Rectangle]] = {}

    # TODO: Move misassembly types to separate row.
    # TODO: Reorder tracks so in order of appearance.
    # TODO: Add label to mapq and bin track.
    # TODO: Rename bin to sequence composition/identity and add percentage.

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

        # Last axis with coverage plt.
        ax: matplotlib.axes.Axes = axs[number_of_overlap_beds]
        # Add bed regions.
        for i, (name, regions) in enumerate(overlay_regions.items()):
            # Make axis as minimal as possible.
            bed_axs: matplotlib.axes.Axes = axs[i]
            minimalize_ax(bed_axs)

            # Map uniq types to new color if none given.
            uniq_types = sorted({r.desc for r in regions if r.desc})
            cmap = dict(
                zip(
                    uniq_types,
                    (np.random.rand(3) for _ in range(len(uniq_types))),
                )
            )
            patches: list[ptch.Rectangle] = []
            for row in regions:
                # Skip rows not within bounds of df.
                if not itv.overlaps(row.region):
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
                patches.append(rect)

            bed_axs.add_collection(PatchCollection(patches, match_original=True))
            subplot_patches[name] = patches
    else:
        fig, axs = plt.subplots(
            2,
            1,
            figsize=(PLOT_WIDTH, PLOT_HEIGHT),
            layout="constrained",
            gridspec_kw={"height_ratios": [3, 0.2]},
        )
        ax = axs[0]

    ax.plot(
        df_pileup["pos"],
        df_pileup["indel"],
        marker="o",
        linestyle="None",
        markersize=2,
        color="purple",
        label="Indels",
    )
    if "mismatch" in df_pileup.columns:
        ax.plot(
            df_pileup["pos"],
            df_pileup["mismatch"],
            marker="o",
            markersize=2,
            linestyle="None",
            color="red",
            label="Mismatches",
        )
    ax.plot(
        df_pileup["pos"],
        df_pileup["cov"],
        marker="o",
        linestyle="None",
        markersize=2,
        color="black",
        label="Coverage",
    )
    # Add misassembly rect patches to highlight region.
    for region in df_misasm.iter_rows(named=True):
        if region["name"] == "good":
            continue
        try:
            rgb: str = region["itemRgb"]
            rgb_color = [int(c) / 255 for c in rgb.split(",")]
        except ValueError:
            continue
        st, end = region["chromStart"], region["chromEnd"]
        ax.axvspan(
            st,
            end,
            color=rgb_color,
            alpha=0.4,
            label=region["name"],
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
        fancybox=False,
        frameon=False,
        prop={"size": 12},
    )
    # Add legends for each overlapped bedfile.
    for i, (name, sp_patches) in enumerate(subplot_patches.items()):
        # Remove plot elements from legend ax.
        # Offset by 2 for coverage plot and its legend.
        legend_ax: matplotlib.axes.Axes = axs[i + number_of_overlap_beds + 2]
        minimalize_ax(legend_ax)

        # Filter rectangle patches.
        sp_patch_labels = set()
        sp_filtered_patches = []
        for patch in sp_patches:
            patch_lbl = patch.get_label()
            if patch_lbl in sp_patch_labels:
                continue
            sp_filtered_patches.append(patch)
            sp_patch_labels.add(patch_lbl)

        # Sort patches by label.
        sp_filtered_patches.sort(key=lambda p: p.get_label())

        # Don't allow title to be just index.
        try:
            _ = int(name)
            title = None
        except ValueError:
            title = name

        legend_ax.legend(
            handles=sp_filtered_patches,
            loc="center",
            alignment="left",
            title=title,
            # At least 1-15 columns in legend.
            ncols=min(max(len(sp_filtered_patches), 1), 15),
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

    ax.set_ylim(0, df_pileup["cov"].mean() * 3)
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
