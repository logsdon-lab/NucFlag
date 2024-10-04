from collections import defaultdict


import polars as pl
from intervaltree import Interval, IntervalTree

from .common import consecutive
from ..misassembly import Misassembly


def identify_gaps(
    df_cov: pl.DataFrame,
    misassemblies: defaultdict[Misassembly, IntervalTree],
    *,
    max_allowed_gap_size: int,
) -> None:
    # Classify gaps.
    df_gaps = df_cov.filter(pl.col("first") == 0)
    # Group consecutive positions.
    for grp in consecutive(df_gaps["position"], stepsize=1):
        if grp.size == 0:
            continue
        # In case of 1 len interval
        if grp[0] == grp[-1]:
            gap_len = 1
            gap = Interval(grp[0], grp[0] + 1)
        else:
            gap_len = grp[-1] - grp[0]
            gap = Interval(grp[0], grp[-1])

        if gap_len < max_allowed_gap_size:
            continue

        misassemblies[Misassembly.GAP].add(gap)
