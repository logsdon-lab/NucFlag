import logging
import argparse
import polars as pl
import pyideogram as pyid
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from typing import Any
from matplotlib.axes import Axes
from matplotlib import cbook
from matplotlib.colors import rgb2hex
from matplotlib.patches import FancyBboxPatch, Patch
from pyideogram.matplotlib_extension import SideRound

from ..common import BED9_COLS, minimalize_ax

BAND_COLORS = pyid.BANDCOL | {"none": (1.0, 1.0, 1.0)}
LBL_KWARGS = dict(rotation=0, ha="right", va="center")


logger = logging.getLogger(__name__)


def draw_ideogramh_no_cytobands(
    ax: Axes, chrom: str, length: int, textkwargs: dict[str, Any] = {}, **ideokwargs
):
    ideokwargs = cbook.normalize_kwargs(ideokwargs, Patch)
    ideokwargs = {"edgecolor": "k"} | ideokwargs
    barpatches = ax.barh(
        [chrom],
        length,
        1,
        0,
        color=BAND_COLORS["gneg"],
        **ideokwargs,
    )
    cornerref = (0, 3, 1, 2)
    for val in [
        "edgecolor",
        "linewidth",
        "hatch",
        "xerr",
        "yerr",
        "error_kw",
        "ecolor",
        "capsize",
        "orientation",
    ]:
        ideokwargs.pop(val, None)

    fs = 5
    textkwargs = {
        "path_effects": [
            pe.withStroke(linewidth=fs / 3, foreground="white", alpha=0.8)
        ],
        "fontsize": fs,
        "va": "center_baseline",
        "ha": "center",
        "weight": "bold",
        "clip_on": True,
    } | textkwargs

    figW, figH = ax.get_figure().get_size_inches()
    _, _, w, h = ax.get_position().bounds
    disp_ratio = (figW * w) / (figH * h)
    ratio = ax.get_data_ratio() * disp_ratio
    # Get patch
    patch = barpatches[0]
    bb = patch.get_bbox()
    patch.remove()
    Box = SideRound(
        xround=min([((0.55 * 0.9) / ratio) / bb.width, 2]),
        yround=0.9,
        corners=cornerref,
    )
    fp = FancyBboxPatch(
        (bb.xmin, bb.ymin),
        abs(bb.width),
        abs(bb.height),
        boxstyle=Box,
        ec=patch.get_edgecolor(),
        fc=patch.get_facecolor(),
        linewidth=patch.get_linewidth(),
        hatch=patch.get_hatch(),
        label=patch.get_label(),
    )
    fp._internal_update(ideokwargs)
    ax.add_patch(fp)


def create_ideogram(args: argparse.Namespace) -> int:
    df_calls = pl.read_csv(
        args.infile,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        schema=dict(BED9_COLS),
    ).with_columns(length=pl.col("chromEnd") - pl.col("chromStart"))

    df_fai = df_calls.group_by(["#chrom"]).agg(
        length=pl.col("chromEnd").max() - pl.col("chromStart").min()
    )
    fai_map = dict(df_fai.iter_rows())
    if args.cytobands:
        cytobands = pyid.dataloader.load_cytobands(args.cytobands)
    else:
        cytobands = None

    color_key = {
        name: rgb2hex([int(e) / 255.0 for e in itemRgb.split(",")])
        if not itemRgb.startswith("#")
        else itemRgb
        for name, itemRgb in df_calls.select("name", "itemRgb").unique().iter_rows()
    }
    max_length = df_fai["length"].max()
    chrom_names = (
        df_fai.filter(pl.col("length") > args.filter_length)
        .sort(by="length", descending=True)["#chrom"]
        .unique(maintain_order=True)
    )
    logger.info(f"Plotting {len(chrom_names)} chromosomes.")
    if chrom_names.is_empty():
        return 1

    # chrom, calls, spacer
    base_height_ratios = [1.0, 0.66]
    width_ratios = [0.1, 0.9]
    num_tracks = len(base_height_ratios)
    height_ratios = base_height_ratios * len(chrom_names)
    # Add extra for legend
    height_ratios.append(3.0)

    fig, axes = plt.subplots(
        ncols=len(width_ratios),
        # 1 additional track for legend
        nrows=len(chrom_names) * num_tracks + 1,
        figsize=(20, len(chrom_names) * args.track_height),
        height_ratios=height_ratios,
        width_ratios=width_ratios,
    )
    fig.subplots_adjust(wspace=0.02)

    ax_col_idx_stats = 0
    ax_col_idx_chrom = 1

    for chrom_name_idx, chrom_name in enumerate(chrom_names):
        df_chrom_calls = df_calls.filter(pl.col("#chrom") == chrom_name)
        chrom_length = fai_map[chrom_name]

        logger.info(
            f"On contig #{chrom_name_idx + 1} {chrom_name} ({chrom_length // 1_000_000} Mbp) ..."
        )
        ax_row_idx_chrom = chrom_name_idx * num_tracks
        ax_row_indices_tracks = range(
            ax_row_idx_chrom + 1, ax_row_idx_chrom + num_tracks
        )
        ax_chrom: Axes = axes[ax_row_idx_chrom, ax_col_idx_chrom]
        ax_chrom_stats: Axes = axes[ax_row_idx_chrom, ax_col_idx_stats]

        ax_chrom.xaxis.set_tick_params(which="both", length=0, labelleft=False)
        ax_chrom.yaxis.set_tick_params(which="both", length=0)
        # pyideogram removes xyticks
        minimalize_ax(ax_chrom)
        minimalize_ax(ax_chrom_stats, remove_ticks=True)

        for ax_row_idx in ax_row_indices_tracks:
            ax_track: Axes = axes[ax_row_idx, ax_col_idx_chrom]
            ax_stats: Axes = axes[ax_row_idx, ax_col_idx_stats]

            ax_track.set_xlim(0, max_length)
            ax_track.set_ylim(0, 1)
            minimalize_ax(ax_track, remove_ticks=True)
            minimalize_ax(ax_stats, remove_ticks=True)

            # Write stats
            type_counts = dict(
                df_chrom_calls.group_by(["name"])
                .agg(length=pl.col("length").sum())
                .iter_rows()
            )
            ax_stats.pie(
                type_counts.values(),
                colors=[color_key[typ] for typ in type_counts.keys()],
                radius=1,
            )

            # Write regions
            for row in df_chrom_calls.iter_rows(named=True):
                color = color_key[row["name"]]
                ax_track.axvspan(
                    xmin=row["chromStart"], xmax=row["chromEnd"], color=color
                )

        ax_chrom.set_xlim(0, max_length)
        if cytobands:
            pyid.ideogramh(
                chrom=chrom_name,
                bands=cytobands,
                ax=ax_chrom,
                color=BAND_COLORS,
                label="",
            )
        else:
            draw_ideogramh_no_cytobands(
                ax=ax_chrom, chrom=chrom_name, length=chrom_length
            )

        ax_chrom.set_yticks([], [])
        ax_chrom.set_ylabel(chrom_name, **LBL_KWARGS)

    # Add legend.
    ax_legend: Axes = axes[len(chrom_names) * num_tracks, ax_col_idx_chrom]
    minimalize_ax(ax_legend, remove_ticks=True)
    minimalize_ax(
        axes[len(chrom_names) * num_tracks, ax_col_idx_stats], remove_ticks=True
    )
    ax_legend.legend(
        handles=[Patch(facecolor=color, label=lbl) for lbl, color in color_key.items()],
        bbox_to_anchor=(0.5, 0.5),
        loc="center",
        ncol=6,
        frameon=False,
        edgecolor="black",
        handlelength=0.7,
        handleheight=0.7,
    )
    # Reduce white space between haps
    logger.info(f"Saving to {args.output_prefix}.(pdf|png)")
    fig.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight", dpi=600)
    fig.savefig(f"{args.output_prefix}.png", bbox_inches="tight", dpi=600)

    return 0
