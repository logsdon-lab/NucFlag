import polars as pl
from collections import defaultdict
from intervaltree import Interval, IntervalTree

from .common import consecutive
from ..misassembly import Misassembly


def get_secondary_allele_coords(
    df_cov: pl.DataFrame,
    *,
    second_thr: int,
    group_distance: int,
    min_group_size: int,
    thr_het_ratio: float,
) -> IntervalTree:
    df_second_outliers = (
        df_cov.filter(pl.col("second") > second_thr)
        .with_columns(het_ratio=pl.col("second") / (pl.col("first") + pl.col("second")))
        .with_columns(
            het_cls=pl.when(pl.col("het_ratio") > thr_het_ratio)
            .then(pl.lit("Error"))
            .otherwise(pl.lit("Het"))
        )
    )

    het_coords = IntervalTree()
    for grp in consecutive(
        df_second_outliers.filter(pl.col("het_cls") == "Het").get_column("position"),
        stepsize=group_distance,
    ):
        if len(grp) < min_group_size:
            continue
        start, end = grp[0], grp[-1]
        ht = (
            df_second_outliers.filter(pl.col("position").is_between(start, end))
            .get_column("second")
            .median()
        )
        het_coords.add(Interval(start, end, (Misassembly.HET, ht)))

    err_coords: IntervalTree = IntervalTree()
    for grp in consecutive(
        df_second_outliers.filter(pl.col("het_cls") == "Error").get_column("position"),
        stepsize=group_distance,
    ):
        if len(grp) < min_group_size:
            continue

        start, end = grp[0], grp[-1]
        ht = (
            df_second_outliers.filter(pl.col("position").is_between(start, end))
            .get_column("second")
            .median()
        )
        err_interval = Interval(start, end, (Misassembly.COLLAPSE_OTHER, ht))
        # Remove hets.
        het_coords.chop(err_interval.begin, err_interval.end)

        err_coords.add(err_interval)

    return err_coords.union(het_coords)


def identify_collapses(
    peaks: IntervalTree,
    second_outliers_coords: IntervalTree,
    classified_second_outliers: set[Interval],
    misassemblies: defaultdict[Misassembly, IntervalTree],
    *,
    collapse_height_thr: int,
) -> None:
    # Intersect intervals and classify collapses.
    for peak in peaks.iter():
        height = peak.data
        # If height of suspected collapsed region is greater than thr, is a collapse.
        if height < collapse_height_thr:
            continue

        overlaps: set[Interval] = second_outliers_coords.overlap(peak)
        for overlap in overlaps:
            new_overlap_interval = Interval(
                min(peak.begin, overlap.begin), max(peak.end, overlap.end)
            )
            misassemblies[Misassembly.COLLAPSE_VAR].add(new_overlap_interval)
            classified_second_outliers.add(overlap)

        if not overlaps:
            misassemblies[Misassembly.COLLAPSE].add(peak)
