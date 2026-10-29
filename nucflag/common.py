import polars as pl
from matplotlib.axes import Axes

BED9P_COLS = [
    ("#chrom", pl.String),
    ("chromStart", pl.UInt64),
    ("chromEnd", pl.UInt64),
    ("name", pl.String),
    ("score", pl.UInt64),
    ("strand", pl.String),
    ("thickStart", pl.UInt64),
    ("thickEnd", pl.UInt64),
    ("itemRgb", pl.String),
    ("zscore", pl.Float32),
    ("af", pl.Float32),
]

STATUSES = (
    "correct",
    "insertion",
    "deletion",
    "softclip",
    "low_quality",
    "het_or_mismap",
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
NON_ERROR_STATUSES = ("scaffold", "het_or_mismap")
PRESETS = ("ont_r9", "ont_r10", "hifi")
DEFAULT_WG_WINDOW = 10_000_000
CORRECT_ITEM_RGB = "206,206,206"


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


def group_dataframe_by_contiguous_itvs(df_region: pl.DataFrame) -> pl.DataFrame:
    """
    Group contiguous intervals. Any break creates a new group id.
    """

    rows = list(
        df_region.sort(by=["#chrom", "chromStart"])
        .unique(maintain_order=True)
        .iter_rows(named=True)
    )
    group_n = 0
    grouped_rows = []
    for i, row in enumerate(rows):
        try:
            next_row = rows[i + 1]
        except IndexError:
            row["group"] = group_n
            grouped_rows.append(row)
            continue

        row["group"] = group_n
        if row["chromEnd"] != next_row["chromStart"]:
            group_n += 1

        grouped_rows.append(row)

    schema = df_region.schema | {"group": pl.UInt64}
    return (
        pl.DataFrame(grouped_rows, orient="row", schema=schema)
        .with_columns(
            length=pl.col("chromEnd") - pl.col("chromStart"),
            minStart=pl.col("chromStart").min().over(["#chrom", "group"]),
            maxEnd=pl.col("chromEnd").max().over(["#chrom", "group"]),
        )
        .with_columns(group_length=pl.col("maxEnd") - pl.col("minStart"))
    )
