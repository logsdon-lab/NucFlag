import polars as pl
from collections import defaultdict
from intervaltree import Interval, IntervalTree

from nucflag.classifier.common import consecutive

from ..misassembly import Misassembly


def identify_zero_cov_regions(
    df_cov: pl.DataFrame,
    misassemblies: defaultdict[Misassembly, IntervalTree],
) -> None:
    # Classify regions with no coverage.
    df_no_cov = df_cov.filter(pl.col("first") == 0)
    # Group consecutive positions.
    for region in consecutive(df_no_cov["position"], stepsize=1):
        if region.size == 0:
            continue
        # In case of 1 len interval
        if region[0] == region[-1]:
            region_itv = Interval(region[0], region[0] + 1)
        else:
            region_itv = Interval(region[0], region[-1])

        misassemblies[Misassembly.MISJOIN].add(region_itv)


def identify_misjoins(
    valleys: IntervalTree,
    second_outliers_coords: IntervalTree,
    classified_second_outliers: set[Interval],
    misassemblies: defaultdict[Misassembly, IntervalTree],
    *,
    misjoin_height_thr: float,
) -> None:
    # Classify misjoins.
    for valley in valleys.iter():
        second_overlaps = second_outliers_coords.overlap(valley)

        # Filter on relative height of valley.
        # Only allow valleys where max change in y doesn't account for half of the valley's depth
        depth = valley.data
        valley_below_thr = depth > misjoin_height_thr

        # Check for overlaps of valleys with regions with high secondary base support.
        for overlap in second_overlaps:
            het_cls, second_overlap_ht = overlap.data

            # Adjust misjoin_height_thr by how much overlap height is.
            # \/    ||
            #    =  \/
            # /\
            adj_valley_below_thr = (depth + second_overlap_ht) > misjoin_height_thr
            # Merge intervals.
            new_overlap_interval = Interval(
                min(valley.begin, overlap.begin), max(valley.end, overlap.end)
            )
            valid_overlap_intervals: list[tuple[Misassembly, Interval]] = []
            if adj_valley_below_thr:
                valid_overlap_intervals.append(
                    (Misassembly.MISJOIN, new_overlap_interval)
                )
                classified_second_outliers.add(overlap)
            else:
                valid_overlap_intervals.append((het_cls, new_overlap_interval))
                classified_second_outliers.add(overlap)

            for misassembly, misassembly_interval in valid_overlap_intervals:
                # Trim back collapses
                misassemblies[Misassembly.COLLAPSE].chop(
                    misassembly_interval.begin, misassembly_interval.end
                )
                misassemblies[Misassembly.COLLAPSE_VAR].chop(
                    misassembly_interval.begin, misassembly_interval.end
                )
                misassemblies[misassembly].add(misassembly_interval)

        if not second_overlaps and valley_below_thr:
            # Trim back collapses
            misassemblies[Misassembly.COLLAPSE].chop(valley.begin, valley.end)
            misassemblies[Misassembly.COLLAPSE_VAR].chop(valley.begin, valley.end)
            # Add misjoin.
            misassemblies[Misassembly.MISJOIN].add(valley)
