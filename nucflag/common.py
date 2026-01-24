import polars as pl
from matplotlib.axes import Axes

BED9_COLS = [
    ("#chrom", pl.String),
    ("chromStart", pl.UInt64),
    ("chromEnd", pl.UInt64),
    ("name", pl.String),
    ("score", pl.UInt64),
    ("strand", pl.String),
    ("thickStart", pl.UInt64),
    ("thickEnd", pl.UInt64),
    ("itemRgb", pl.String),
]

STATUSES = (
    "correct",
    "insertion",
    "deletion",
    "softclip",
    "het_mismap",
    "collapse",
    "misjoin",
    "mismatch",
    "false_dup",
    "homopolymer",
    "dinucleotide",
    "simple_repeat",
    "other_repeat",
    "scaffold",
)
PRESETS = ("ont_r9", "ont_r10", "hifi")


def minimalize_ax(ax: Axes, *, remove_ticks: bool = False) -> None:
    for spine in ["left", "right", "bottom", "top"]:
        ax.spines[spine].set_visible(False)
    if remove_ticks:
        ax.tick_params(
            axis="both",
            left=False,
            top=False,
            right=False,
            bottom=False,
            labelleft=False,
            labeltop=False,
            labelright=False,
            labelbottom=False,
        )


def add_group_columns(df_region: pl.DataFrame) -> pl.DataFrame:
    return (
        df_region.sort(by=["#chrom", "chromStart"])
        .with_columns(
            group=(pl.col("chromEnd") != pl.col("chromStart").shift(-1))
            .fill_null(False)
            .rle_id()
            .over("#chrom"),
        )
        .with_columns(
            group=pl.when(pl.col("group") % 2 != 0)
            .then(pl.col("group") - 1)
            .otherwise(pl.col("group"))
            .over("#chrom")
        )
        .with_columns(
            length=pl.col("chromEnd") - pl.col("chromStart"),
            minStart=pl.col("chromStart").min().over(["#chrom", "group"]),
            maxEnd=pl.col("chromEnd").max().over(["#chrom", "group"]),
        )
    )
