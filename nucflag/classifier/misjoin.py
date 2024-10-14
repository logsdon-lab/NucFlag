from collections import defaultdict
from intervaltree import Interval, IntervalTree

from .common import calculate_het_ratio, subtract_interval
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
        overlaps_gap = misassemblies[Misassembly.GAP].overlap(valley)
        overlaps_collapse = misassemblies[Misassembly.COLLAPSE].overlaps(
            valley
        ) or misassemblies[Misassembly.COLLAPSE_VAR].overlaps(valley)

        # Ignore if overlaps existing gap, collapse, or collapse with variant.
        if overlaps_collapse:
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
            valid_overlap_intervals: list[tuple[Misassembly, Interval]] = []
            if adj_valley_below_thr:
                valid_overlap_intervals.append(
                    (Misassembly.MISJOIN, new_overlap_interval)
                )
                classified_second_outliers.add(overlap)
            else:
                het_ratio = calculate_het_ratio(overlap)
                if het_ratio >= het_ratio_thr:
                    valid_overlap_intervals.append(
                        (Misassembly.ERROR, new_overlap_interval)
                    )
                    classified_second_outliers.add(overlap)

            # Subtract interval by gaps.
            if overlaps_gap and valid_overlap_intervals:
                (
                    valid_overlap_misassembly,
                    valid_overlap_interval,
                ) = valid_overlap_intervals.pop()
                for interval in subtract_interval(
                    valid_overlap_interval, by=overlaps_gap
                ):
                    valid_overlap_intervals.append(
                        (valid_overlap_misassembly, interval)
                    )

            for misassembly, misassembly_interval in valid_overlap_intervals:
                misassemblies[misassembly].add(misassembly_interval)

        if not second_overlaps and valley_below_thr:
            valid_overlap_intervals = []
            # Subtract interval by gaps.
            if overlaps_gap:
                for interval in subtract_interval(valley, by=overlaps_gap):
                    valid_overlap_intervals.append(interval)
            else:
                valid_overlap_intervals.append(valley)

            for interval in valid_overlap_intervals:
                misassemblies[Misassembly.MISJOIN].add(interval)
