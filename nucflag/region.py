from enum import StrEnum, auto
from typing import Generator, NamedTuple

import polars as pl
from intervaltree import Interval


class RegionStatus(StrEnum):
    MISASSEMBLED = auto()
    GOOD = auto()


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


class Region(NamedTuple):
    name: str
    region: Interval
    desc: str | None
    action: Action | None


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
        .with_columns(
            mapq_color=pl.when(pl.col("mapq_rng") == "0")
            .then(pl.lit("#666666"))
            .when(pl.col("mapq_rng") == "1-5")
            .then(pl.lit("#8f59a7"))
            .when(pl.col("mapq_rng") == "5-10")
            .then(pl.lit("#5954a8"))
            .when(pl.col("mapq_rng") == "10-15")
            .then(pl.lit("#01aef3"))
            .when(pl.col("mapq_rng") == "15-20")
            .then(pl.lit("#04b99e"))
            .when(pl.col("mapq_rng") == "20-25")
            .then(pl.lit("#8bc83b"))
            .when(pl.col("mapq_rng") == "25-30")
            .then(pl.lit("#cdde25"))
            .when(pl.col("mapq_rng") == "30-35")
            .then(pl.lit("#fff600"))
            .when(pl.col("mapq_rng") == "35-40")
            .then(pl.lit("#ffc309"))
            .when(pl.col("mapq_rng") == "40-45")
            .then(pl.lit("#fa931a"))
            .when(pl.col("mapq_rng") == "45-50")
            .then(pl.lit("#f8631f"))
            .when(pl.col("mapq_rng") == "50-55")
            .then(pl.lit("#f21821"))
            .otherwise(pl.lit("#FF8DA1"))
        )
        .select("st", "end", "mapq_rng", "mapq_color")
    )
    for st, end, mapq_rng, mapq_color in regions.iter_rows():
        yield Region(
            name,
            Interval(st, end),
            desc=mapq_rng,
            action=Action(ActionOpt.PLOT, mapq_color),
        )
