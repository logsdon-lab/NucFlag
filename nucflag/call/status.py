import polars as pl

from typing import Literal
from collections import Counter, defaultdict
from intervaltree import Interval, IntervalTree

from ..common import STATUSES, group_dataframe_by_contiguous_itvs


def group_dataframe_by_region(
    df_region: pl.DataFrame,
    bed_group_by_regions: str,
    groupby: Literal["region", "name"],
) -> pl.DataFrame:
    df_bed_regions = pl.read_csv(
        bed_group_by_regions,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        new_columns=["#chrom", "chromStart", "chromEnd"],
    ).with_columns(length=pl.col("chromEnd") - pl.col("chromStart"))

    itrees: defaultdict[str, IntervalTree] = defaultdict(IntervalTree)
    if groupby == "name":
        df_bed_regions = (
            df_bed_regions.rename({"column_4": "name"})
            .unique(subset=["#chrom", "chromStart", "chromEnd", "name"])
            .with_columns(group_length=pl.col("length").sum().over("name"))
        )
    else:
        df_bed_regions = df_bed_regions.unique(
            subset=["#chrom", "chromStart", "chromEnd"]
        ).with_columns(group_length=pl.col("length"))

    group_lengths: Counter[str | int] = Counter()
    for i, row in enumerate(df_bed_regions.iter_rows(named=True)):
        group_id = row.get("name", i)
        group_length = row["chromEnd"] - row["chromStart"]
        group_lengths[group_id] += group_length
        itrees[row["#chrom"]].add(
            Interval(row["chromStart"], row["chromEnd"], group_id)
        )

    new_rows = []
    for row in df_region.iter_rows(named=True):
        itree_grp = itrees.get(row["#chrom"])
        if not itree_grp:
            row["group"] = None
            new_rows.append(row)
            continue
        ovl = itree_grp.overlap(row["chromStart"], row["chromEnd"])
        if not ovl:
            row["group"] = None
            new_rows.append(row)
            continue

        # Overlaps. Trim to boundaries of interval
        for ovl_itv in sorted(ovl):
            st = max(row["chromStart"], ovl_itv.begin)
            end = min(row["chromEnd"], ovl_itv.end)
            new_row = row.copy()
            new_row["chromStart"] = st
            new_row["chromEnd"] = end
            new_row["thickStart"] = st
            new_row["thickEnd"] = end
            new_row["group"] = ovl_itv.data
            new_rows.append(new_row)

    df_rows = (
        pl.DataFrame(new_rows, orient="row")
        .with_columns(
            length=pl.col("chromEnd") - pl.col("chromStart"),
            minStart=pl.col("chromStart").min().over(["#chrom", "group"]),
            maxEnd=pl.col("chromEnd").max().over(["#chrom", "group"]),
        )
        .join(
            pl.DataFrame(
                list(group_lengths.items()),
                schema=["group", "group_length"],
                orient="row",
            ),
            on="group",
            how="left",
        )
    )
    return df_rows


def generate_status_from_regions(
    df_region: pl.DataFrame,
    groupby: Literal["region", "name"],
    bed_group_by_regions: str | None = None,
) -> pl.DataFrame:
    if bed_group_by_regions:
        df_region_grp = group_dataframe_by_region(
            df_region, bed_group_by_regions, groupby=groupby
        )
    else:
        # Regions don't have coordinates so need to group by break in contiguity of adjacent intervals.
        df_region_grp = group_dataframe_by_contiguous_itvs(df_region).with_columns(
            group_length=pl.col("maxEnd") - pl.col("minStart")
        )

    if groupby == "region":
        cols_agg = {
            "chromStart": pl.col("minStart").first(),
            "chromEnd": pl.col("maxEnd").first(),
            "perc": (pl.col("length").sum() / pl.col("group_length").first()) * 100.0,
        }
        cols_pivot_index = ["#chrom", "chromStart", "chromEnd"]
        cols_select = [
            pl.col("#chrom"),
            pl.col("chromStart"),
            pl.col("chromEnd"),
            pl.col("status"),
            *[pl.col(status) for status in STATUSES],
        ]
    else:
        cols_agg = {
            "perc": (pl.col("length").sum() / pl.col("group_length").first()) * 100.0,
        }
        cols_pivot_index = ["#chrom", "group", "group_length"]
        cols_select = [
            pl.col("#chrom"),
            pl.col("group"),
            pl.col("group_length"),
            pl.col("status"),
            *[pl.col(status) for status in STATUSES],
        ]

    df_final = (
        df_region_grp.group_by(["#chrom", "name", "group", "group_length"])
        .agg(**cols_agg)
        .pivot(
            on="name",
            index=cols_pivot_index,
            values="perc",
            maintain_order=True,
        )
        # Ensure column exists.
        # https://github.com/pola-rs/polars/issues/18372#issuecomment-2390371173
        .with_columns(
            **{status: pl.coalesce(pl.col(f"^{status}$"), 0.0) for status in STATUSES},
        )
        .with_columns(
            status=pl.when(pl.col("correct") == 100.0)
            .then(pl.lit("correct"))
            .otherwise(pl.lit("misassembled")),
        )
        .select(cols_select)
        .fill_null(0.0)
    )
    if groupby == "region":
        df_final = (
            df_final
            # chromStart and chromEnd get cast to float for some reason.
            .cast({"chromStart": pl.UInt64, "chromEnd": pl.UInt64}).sort(
                by=["#chrom", "chromStart"]
            )
        )
    else:
        df_final = df_final.cast({"group_length": pl.UInt64})
    return df_final
