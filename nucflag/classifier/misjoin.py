from collections import defaultdict
from intervaltree import Interval, IntervalTree

from .common import calculate_het_ratio
from ..misassembly import Misassembly


def identify_misjoins(
    valleys: IntervalTree,
    second_outliers_coords: IntervalTree,
    classified_second_outliers: set[Interval],
    misassemblies: defaultdict[Misassembly, IntervalTree],
    *,
    misjoin_height_thr: float,
    het_ratio_thr: float,
) -> None:
    # Classify misjoins.
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
        depth = valley.data
        valley_below_thr = depth > misjoin_height_thr

        # Check for overlaps of valleys with regions with high secondary base support.
        for overlap in second_overlaps:
            _, second_overlap_ht = overlap.data

            # Adjust misjoin_height_thr by how much overlap height is.
            # \/    ||
            #    =  \/
            # /\
            adj_valley_below_thr = (depth + second_overlap_ht) > misjoin_height_thr
            # Merge intervals.
            new_overlap_interval = Interval(
                min(valley.begin, overlap.begin), max(valley.end, overlap.end)
            )
            if adj_valley_below_thr:
                misassemblies[Misassembly.MISJOIN].add(new_overlap_interval)
                classified_second_outliers.add(overlap)
            else:
                het_ratio = calculate_het_ratio(overlap)
                if het_ratio >= het_ratio_thr:
                    misassemblies[Misassembly.ERROR].add(new_overlap_interval)
                    classified_second_outliers.add(overlap)

        if not second_overlaps and valley_below_thr:
            misassemblies[Misassembly.MISJOIN].add(valley)
