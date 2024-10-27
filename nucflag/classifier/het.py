from collections import defaultdict

from intervaltree import Interval, IntervalTree

from ..misassembly import Misassembly


def identify_hets(
    second_outliers_coords: IntervalTree,
    classified_second_outliers: set[Interval],
    misassemblies: defaultdict[Misassembly, IntervalTree],
):
    # Assign heterozygous site for remaining secondary regions not categorized.
    for second_outlier in second_outliers_coords:
        if second_outlier in classified_second_outliers:
            continue
        het_cls, _ = second_outlier.data
        misassemblies[het_cls].add(second_outlier)
