from enum import StrEnum, auto
from typing import Generator, NamedTuple

import polars as pl
import matplotlib.colors

from bisect import bisect
from intervaltree import Interval  # type: ignore[import-untyped]


class ActionOpt(StrEnum):
    IGNORE = auto()
    PLOT = auto()
    NOOP = auto()


class Action(NamedTuple):
    opt: ActionOpt
    desc: str | None

    def as_str(self) -> str:
        action = f"{self.opt}"
        if self.desc:
            action = f"{action}:{self.desc}"
        return action


class Region(NamedTuple):
    name: str
    region: Interval
    desc: str | None
    action: Action | None

    def as_tsv(self) -> str:
        action_str = self.action.as_str() if self.action else "None"
        return f"{self.name}\t{self.region.begin}\t{self.region.end}\t{self.desc}\t{action_str}"


DF_MAPQ_COLORS = pl.DataFrame(
    {
        "mapq_rng": [
            "0",
            "1-5",
            "5-10",
            "10-15",
            "15-20",
            "20-25",
            "25-30",
            "30-35",
            "35-40",
            "40-45",
            "45-50",
            "50-55",
        ],
        "mapq_color": [
            "#666666",
            "#8f59a7",
            "#5954a8",
            "#01aef3",
            "#04b99e",
            "#8bc83b",
            "#cdde25",
            "#fff600",
            "#ffc309",
            "#fa931a",
            "#f8631f",
            "#f21821",
        ],
    }
)


def add_mapq_overlay_region(
    name: str, df: pl.DataFrame
) -> Generator[Region, None, None]:
    regions = (
        df.with_columns(mapq_grp=pl.col("mapq").rle_id())
        .group_by(["mapq_grp"])
        .agg(
            st=pl.col("pos").min(), end=pl.col("pos").max(), mapq=pl.col("mapq").first()
        )
        .with_columns(mapq_lower_bound=(pl.col("mapq") // 5) * 5)
        .with_columns(
            mapq_rng=pl.when(pl.col("mapq_lower_bound") == 60)
            .then(pl.lit("55-60"))
            # 0
            .when(pl.col("mapq") == 0)
            .then(pl.lit("0"))
            # 0-5 to 1-5
            .when(pl.col("mapq_lower_bound") == 0)
            .then(pl.lit("1-5"))
            # 5-60
            .otherwise(
                pl.col("mapq_lower_bound").cast(pl.String)
                + "-"
                + (pl.col("mapq_lower_bound") + 5).clip(0, 60).cast(pl.String)
            )
        )
        .join(DF_MAPQ_COLORS, on=["mapq_rng"], how="left")
        .with_columns(pl.col("mapq_color").fill_null(pl.lit("#FF8DA1")))
        .select("st", "end", "mapq_rng", "mapq_color")
    )
    for st, end, mapq_rng, mapq_color in regions.iter_rows():
        yield Region(
            name,
            Interval(st, end),
            desc=mapq_rng,
            action=Action(ActionOpt.PLOT, mapq_color),
        )


def add_ident_overlay_region(
    name: str, df: pl.DataFrame, breakpoint_colors: tuple[list[float], list[str]]
) -> Generator[Region, None, None]:
    breakpoints, colors = breakpoint_colors
    regions = (
        df.with_columns(bin_grp=pl.col("bin").rle_id())
        .group_by(["bin_grp"])
        .agg(
            st=pl.col("pos").min(),
            end=pl.col("pos").max(),
            bin_ident=pl.col("bin_ident").mean(),
        )
        .select("st", "end", "bin_ident")
    )
    for st, end, bin_ident in regions.iter_rows():
        if bin_ident == 0.0:
            continue
        try:
            idx_end = bisect(breakpoints, bin_ident)
            ident_end = breakpoints[idx_end]
            if idx_end == 0:
                ident_st = 0.0
            else:
                ident_st = breakpoints[idx_end - 1]
                color = colors[idx_end - 1]
            desc = f"{ident_st}%-{ident_end}%"
        except IndexError:
            ident_end = breakpoints[-1]
            color = colors[-1]
            desc = str(ident_end)

        yield Region(
            name,
            Interval(st, end),
            desc=desc,
            action=Action(ActionOpt.PLOT, color),
        )


def add_misassemblies_overlay_region(df: pl.DataFrame) -> Generator[Region, None, None]:
    df = df.with_columns(
        itemRgb=pl.when(pl.col("itemRgb").str.contains(","))
        .then(
            pl.col("itemRgb")
            .str.split(",")
            .list.eval(pl.element().cast(pl.Float32) / 255.0)
            .map_elements(lambda x: matplotlib.colors.to_hex(x), return_dtype=pl.String)
        )
        .otherwise(pl.col("itemRgb"))
    ).select("chromStart", "chromEnd", "name", "itemRgb")
    for row in df.iter_rows(named=True):
        yield Region(
            row["name"],
            Interval(row["chromStart"], row["chromEnd"]),
            desc=row["name"],
            action=Action(ActionOpt.PLOT, row["itemRgb"]),
        )
