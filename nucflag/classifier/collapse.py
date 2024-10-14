from collections import defaultdict
from intervaltree import Interval, IntervalTree
from ..misassembly import Misassembly


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
