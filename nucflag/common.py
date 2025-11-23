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
    "indel",
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


def minimalize_ax(ax: Axes, *, remove_ticks: bool = False):
    for spine in ["left", "right", "bottom", "top"]:
        ax.spines[spine].set_visible(False)
    if remove_ticks:
        ax.set_xticks([], [])
        ax.set_yticks([], [])
