from enum import StrEnum, auto
from typing import Generator, NamedTuple

import polars as pl
from intervaltree import Interval


class IgnoreOpt(StrEnum):
    ABSOLUTE = auto()
    RELATIVE = auto()


class ActionOpt(StrEnum):
    IGNORE = auto()
    PLOT = auto()
    NOOP = auto()


class Action(NamedTuple):
    opt: ActionOpt
    desc: IgnoreOpt | str | None

    def as_str(self):
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


def add_bin_overlay_region(
    name: str, df: pl.DataFrame
) -> Generator[Region, None, None]:
    regions = (
        df.with_columns(bin_grp=pl.col("bin").rle_id())
        .group_by(["bin_grp"])
        .agg(st=pl.col("pos").min(), end=pl.col("pos").max(), bin=pl.col("bin").first())
        .select("st", "end", "bin")
    )
    for st, end, bin_num in regions.iter_rows():
        yield Region(
            name,
            Interval(st, end),
            desc=f"b{bin_num}",
            action=Action(ActionOpt.PLOT, None),
        )
