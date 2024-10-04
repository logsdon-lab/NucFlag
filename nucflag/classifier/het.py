from collections import defaultdict

import polars as pl
from intervaltree import Interval, IntervalTree

from .common import calculate_het_ratio
from ..misassembly import Misassembly


def identify_hets(
    df_cov: pl.DataFrame,
    second_outliers_coords: IntervalTree,
    classified_second_outliers: set[Interval],
    misassemblies: defaultdict[Misassembly, IntervalTree],
    *,
    second_thr: int,
    het_ratio_thr: int,
):
    # Assign heterozygous site for remaining secondary regions not categorized.
    for second_outlier in second_outliers_coords:
        if second_outlier in classified_second_outliers:
            continue

        het_ratio = calculate_het_ratio(df_cov, second_outlier, second_thr)

        # If under het ratio, treat as het. Otherwise, some error.
        if het_ratio >= het_ratio_thr:
            misassemblies[Misassembly.ERROR].add(second_outlier)
        else:
            misassemblies[Misassembly.HET].add(second_outlier)
