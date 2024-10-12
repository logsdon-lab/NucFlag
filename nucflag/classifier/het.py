from collections import defaultdict

from intervaltree import Interval, IntervalTree

from .common import calculate_het_ratio
from ..misassembly import Misassembly


def identify_hets(
    second_outliers_coords: IntervalTree,
    classified_second_outliers: set[Interval],
    misassemblies: defaultdict[Misassembly, IntervalTree],
    *,
    first_mean: float,
    het_ratio_thr: float,
):
    # Assign heterozygous site for remaining secondary regions not categorized.
    for second_outlier in second_outliers_coords:
        if second_outlier in classified_second_outliers:
            continue

        # If misflagging, may be due to second_thr.
        # 2nd order poly fit will give greater weight to smaller values lowering peak of fn.
        # Increase second_thr fixes these cases.

        # Calculate het ratio of region and using mean of entire contig.
        # Prevents situations where clear misassembly and enrichment of second bases but is missed because first most base is also high.
        het_ratio_region = calculate_het_ratio(second_outlier)
        het_ratio_all = second_outlier.data[1] / (first_mean + second_outlier.data[1])
        # If under het ratios, treat as het. Otherwise, some error.
        if het_ratio_region >= het_ratio_thr or het_ratio_all >= het_ratio_thr:
            misassemblies[Misassembly.ERROR].add(second_outlier)
        else:
            misassemblies[Misassembly.HET].add(second_outlier)
