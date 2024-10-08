import polars as pl
from collections import defaultdict
from intervaltree import Interval, IntervalTree

from .common import calculate_het_ratio
from ..misassembly import Misassembly


def identify_misjoins(
    df_cov: pl.DataFrame,
    valleys: IntervalTree,
    second_outliers_coords: IntervalTree,
    classified_second_outliers: set[Interval],
    misassemblies: defaultdict[Misassembly, IntervalTree],
    *,
    misjoin_height_thr: int,
    second_thr: int,
    het_ratio_thr: int,
) -> None:
    # Classify misjoins.
    # Threshold for max allowed drop between two base positions in coverage.
    thr_max_dy = misjoin_height_thr * 0.5

    for valley in valleys.iter():
        second_overlaps = second_outliers_coords.overlap(valley)
        overlaps_gap = misassemblies[Misassembly.GAP].overlaps(valley)
        overlaps_collapse = misassemblies[Misassembly.COLLAPSE].overlaps(
            valley
        ) or misassemblies[Misassembly.COLLAPSE_VAR].overlaps(valley)

        # Ignore if overlaps existing gap, collapse, or collapse with variant.
        if overlaps_gap or overlaps_collapse:
            continue

        # Filter on relative height of valley.
        # Only allow valleys where max change in y doesn't account for half of the valley's depth
        depth, max_dy = valley.data
        valley_below_thr = depth > misjoin_height_thr and max_dy < thr_max_dy

        for overlap in second_overlaps:
            # Merge intervals.
            new_overlap_interval = Interval(
                min(valley.begin, overlap.begin), max(valley.end, overlap.end)
            )
            if valley_below_thr:
                misassemblies[Misassembly.MISJOIN].add(new_overlap_interval)
                classified_second_outliers.add(overlap)
            else:
                het_ratio = calculate_het_ratio(df_cov, overlap, second_thr)
                if het_ratio >= het_ratio_thr:
                    misassemblies[Misassembly.ERROR].add(new_overlap_interval)
                    classified_second_outliers.add(overlap)

        if not second_overlaps and valley_below_thr:
            misassemblies[Misassembly.MISJOIN].add(valley)
