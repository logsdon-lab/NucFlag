import logging
import argparse
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from matplotlib.axes import Axes
from matplotlib.colors import rgb2hex
from matplotlib.patches import Patch


from ..common import BED9_COLS, STATUSES, minimalize_ax


logger = logging.getLogger(__name__)


def create_breakdown_plot(args: argparse.Namespace) -> int:
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
    color_key = {
        name: rgb2hex([int(e) / 255.0 for e in itemRgb.split(",")])
        if not itemRgb.startswith("#")
        else itemRgb
        for name, itemRgb in df_calls.select("name", "itemRgb").unique().iter_rows()
    }
    chrom_names = (
        df_fai.filter(pl.col("length") > args.filter_length)
        .sort(by="length", descending=True)["#chrom"]
        .unique(maintain_order=True)
    )
    if chrom_names.is_empty():
        return 1

    df_fai = df_fai.filter(pl.col("#chrom").is_in(chrom_names)).cast(
        {"#chrom": pl.Enum(chrom_names)}
    )
    logger.info(f"Plotting {len(chrom_names)} chromosomes.")

    df_grouped_calls = (
        df_calls.filter(pl.col("#chrom").is_in(chrom_names))
        .group_by(["#chrom", "name"])
        .agg(pl.col("length").sum())
        .pivot(
            on="name",
            index="#chrom",
            values="length",
            maintain_order=True,
        )
        .with_columns(
            # Ensure column exists.
            # https://github.com/pola-rs/polars/issues/18372#issuecomment-2390371173
            **{status: pl.coalesce(pl.col(f"^{status}$"), 0.0) for status in STATUSES},
        )
        # TODO: There's got to be a better way to do this.
        .unpivot(index="#chrom", value_name="length", variable_name="name")
        .cast({"name": pl.Enum(STATUSES), "#chrom": pl.Enum(chrom_names)})
        .sort(by=["#chrom", "name"])
    )
    # Scale fig width by number of chroms
    fig, axes = plt.subplots(
        nrows=2,
        height_ratios=[0.9, 0.1],
        figsize=(len(chrom_names), 20),
        layout="constrained",
    )
    ax: Axes = axes[0]
    bottom = np.zeros(len(chrom_names))

    if args.type == "percent":
        ax.set_ylabel("Percent (%)")
        df_grouped_calls = (
            df_grouped_calls.join(df_fai, on="#chrom", how="left")
            .with_columns(length=(pl.col("length") / pl.col("length_right")) * 100.0)
            .drop("length_right")
        )
    else:
        ax.set_ylabel("Length (Mbp)")
        df_grouped_calls = df_grouped_calls.with_columns(pl.col("length") / 1_000_000)

    # https://matplotlib.org/stable/gallery/lines_bars_and_markers/bar_label_demo.html#bar-chart-with-labels
    for group, df_group in df_grouped_calls.group_by(["name"], maintain_order=True):
        group = group[0]
        values = df_group["length"].to_numpy()
        color = color_key.get(group)
        p = ax.bar(chrom_names, values, 0.6, label=group, bottom=bottom, color=color)
        bottom += values
        text_color = "#000000"
        ax.bar_label(
            p,
            label_type="center",
            fmt="{:,.1f}",
            color=text_color,
            path_effects=[pe.withStroke(linewidth=4, foreground="white")],
        )

    ax.set_title("Misassemblies by contig", loc="left", pad=20)
    ax.set_xticklabels(
        ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor"
    )

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    # Add legend.
    legend_ax: Axes = axes[1]
    minimalize_ax(legend_ax, remove_ticks=True)
    legend_ax.legend(
        handles=[Patch(facecolor=color, label=lbl) for lbl, color in color_key.items()],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.5),
        ncol=len(color_key),
        frameon=False,
        edgecolor="black",
        handlelength=0.7,
        handleheight=0.7,
    )
    logger.info(f"Saving to {args.output_prefix}.(pdf|png)")
    fig.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{args.output_prefix}.png", bbox_inches="tight", dpi=300)

    return 0
