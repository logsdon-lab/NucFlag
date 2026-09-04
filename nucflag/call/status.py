import logging
import polars as pl
import polars.selectors as cs

from typing import Any, Literal, Collection, TextIO
from collections import Counter, defaultdict
from intervaltree import Interval, IntervalTree

from ..qv.run import expr_qv
from ..common import STATUSES, group_dataframe_by_contiguous_itvs


logger = logging.getLogger(__name__)


def group_dataframe_by_region(
    df_region: pl.DataFrame,
    bed_group_by_regions: TextIO,
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
        if "column_4" not in df_bed_regions.columns:
            raise ValueError(
                f"Missing name (4th) column in {bed_group_by_regions.name}"
            )

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
    group_rows: Counter[str | int] = Counter()
    for i, row in enumerate(df_bed_regions.iter_rows(named=True)):
        group_id = row.get("name", i)
        group_length = row["chromEnd"] - row["chromStart"]
        group_lengths[group_id] += group_length
        group_rows[group_id] += 1
        itrees[row["#chrom"]].add(
            Interval(row["chromStart"], row["chromEnd"], group_id)
        )

    new_rows_by_chrom_group: defaultdict[
        tuple[str, int], list[dict[str, Any]]
    ] = defaultdict(list)
    for row in df_region.iter_rows(named=True):
        itree_grp = itrees.get(row["#chrom"])
        if not itree_grp:
            continue
        ovl = itree_grp.overlap(row["chromStart"], row["chromEnd"])
        if not ovl:
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
            new_rows_by_chrom_group[(row["#chrom"], ovl_itv.data)].append(new_row)

    # Merge overlaps.
    # In case where input regions have overlap, merge rows based on group/name/overlap
    new_rows = []
    for _, itvs in new_rows_by_chrom_group.items():
        itvs_sorted = sorted(itvs, key=lambda x: x["chromStart"])
        itvs_merged = [itvs_sorted[0]]
        for itv in itvs_sorted[1:]:
            prev_itv = itvs_merged[-1]
            dst = itv["chromStart"] - prev_itv["chromEnd"]
            # Overlaps (Not book-ended) and same name
            if dst < 0 and prev_itv["name"] == itv["name"]:
                prev_itv["chromEnd"] = itv["chromEnd"]
            else:
                itvs_merged.append(itv)
        new_rows.extend(itvs_merged)

    df_rows = (
        pl.DataFrame(new_rows, orient="row", infer_schema_length=None)
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
        .join(
            pl.DataFrame(
                list(group_rows.items()),
                schema=["group", "group_rows"],
                orient="row",
            ),
            on="group",
            how="left",
        )
    )
    return df_rows


def generate_status_from_regions(
    df_region: pl.DataFrame,
    groupby: Literal["region", "name", "all"],
    thr_qv: int,
    ignore_calls_qv: Collection[str],
    metric: Literal[
        "count",
        "length",
    ] = "length",
    bed_group_by_regions: TextIO | None = None,
) -> pl.DataFrame:
    if groupby == "all":
        # Group rows and get length of contigs
        # Overwrite groups to merge everything
        df_region_grp = group_dataframe_by_contiguous_itvs(df_region)
        df_region_grp = df_region_grp.with_columns(
            group=pl.lit("all"),
            group_length=pl.col("group_length").unique().sum(),
            group_rows=pl.col("#chrom").count(),
        )
    elif bed_group_by_regions:
        logger.info(f"Grouping regions by {bed_group_by_regions.name}.")
        df_region_grp = group_dataframe_by_region(
            df_region, bed_group_by_regions, groupby=groupby
        )
    else:
        logger.info("Grouping contiguous regions.")
        # Regions don't have coordinates so need to group by break in contiguity of adjacent intervals.
        df_region_grp = group_dataframe_by_contiguous_itvs(df_region)

    if metric == "count":
        logger.info("Tallying by count. QV and status columns will be omitted.")
        expr_metric = pl.col("name").count()
    else:
        logger.info("Tallying by length.")
        expr_metric = (pl.col("length").sum() / pl.col("group_length").first()) * 100.0

    if groupby == "region":
        cols_groupby = ["#chrom", "name", "group", "group_length"]
        cols_agg = {
            "chromStart": pl.col("minStart").first(),
            "chromEnd": pl.col("maxEnd").first(),
            "metric": expr_metric,
        }
        cols_pivot_index = ["#chrom", "chromStart", "chromEnd"]
        cols_select = {
            "#chrom": pl.col("#chrom"),
            "chromStart": pl.col("chromStart"),
            "chromEnd": pl.col("chromEnd"),
            "status": pl.col("status"),
            **{status: pl.col(status) for status in STATUSES},
            "QV": pl.col("QV"),
        }
    else:
        cols_groupby = ["name", "group", "group_length", "group_rows"]
        cols_agg = {"metric": expr_metric}
        cols_pivot_index = ["group", "group_length", "group_rows"]
        cols_select = {
            "group": pl.col("group"),
            "group_length": pl.col("group_length"),
            "group_rows": pl.col("group_rows"),
            "status": pl.col("status"),
            **{status: pl.col(status) for status in STATUSES},
            "QV": pl.col("QV"),
        }
    error_statuses = set(STATUSES).difference((*ignore_calls_qv, "correct"))
    if metric != "count":
        logger.info(
            f"Treating {error_statuses} as error statuses and regions with QV >= {thr_qv} as correct."
        )
    df_final = (
        df_region_grp.group_by(cols_groupby)
        .agg(**cols_agg)
        .pivot(
            on="name",
            index=cols_pivot_index,
            values="metric",
            maintain_order=True,
        )
        # Ensure column exists.
        # https://github.com/pola-rs/polars/issues/18372#issuecomment-2390371173
        .with_columns(
            **{status: pl.coalesce(pl.col(f"^{status}$"), 0.0) for status in STATUSES},
        )
    )

    if metric == "length":
        df_final = (
            df_final.with_columns(
                # We're adding percents here but doesn't matter since proportion used for QV
                bpError=pl.sum_horizontal(cs.contains(error_statuses)),  # pyright:ignore
                bpTotal=pl.sum_horizontal(cs.contains(*ignore_calls_qv, "correct")),
            )
            .with_columns(QV=expr_qv())
            .with_columns(
                status=pl.when(pl.col("QV") >= pl.lit(thr_qv))
                .then(pl.lit("correct"))
                .otherwise(pl.lit("misassembled")),
            )
            .select(cols_select.values())
        )
    else:
        cols_select.pop("QV")
        cols_select.pop("status")
        df_final = df_final.select(cols_select.values())

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
