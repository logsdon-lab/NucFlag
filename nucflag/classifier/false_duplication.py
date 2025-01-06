import polars as pl
from collections import defaultdict, deque
from typing import Iterable, Callable
from intervaltree import Interval, IntervalTree
from .collapse import Misassembly


def fn_cmp_def(_itv_1: Interval, _itv_2: Interval) -> bool:
    return True


def fn_merge_itv_def(itv_1: Interval, itv_2: Interval) -> Interval:
    return Interval(begin=itv_1.begin, end=itv_2.end, data=None)


def merge_itvs(
    itvs: Iterable[Interval],
    dst: int = 1,
    fn_cmp: Callable[[Interval, Interval], bool] | None = None,
    fn_merge_itv: Callable[[Interval, Interval], Interval] | None = None,
):
    if not fn_cmp:
        fn_cmp = fn_cmp_def
    if not fn_merge_itv:
        fn_merge_itv = fn_merge_itv_def

    final_itvs = []
    sorted_itvs = deque(sorted(itvs))
    while sorted_itvs:
        try:
            itv_1 = sorted_itvs.popleft()
        except IndexError:
            break
        try:
            itv_2 = sorted_itvs.popleft()
        except IndexError:
            final_itvs.append(itv_1)
            break
        dst_between = itv_2.begin - itv_1.end
        passes_cmp = fn_cmp(itv_1, itv_2)
        if dst_between <= dst and passes_cmp:
            sorted_itvs.appendleft(fn_merge_itv(itv_1, itv_2))
        else:
            final_itvs.append(itv_1)
            sorted_itvs.appendleft(itv_2)

    return final_itvs


def identify_false_duplication(
    df_cov: pl.DataFrame,
    valleys: IntervalTree,
    misassemblies: defaultdict[Misassembly, IntervalTree],
    false_dupe_height_thr: float,
    misjoin_abs_height_thr: float,
    bp_merge: int,
    stdev_merge: int,
) -> None:
    # Generate subset of valleys that are potentially false dupes and apply a second merging.
    valid_itvs = []
    for valley in valleys.iter():
        depth, _ = valley.data

        # If pass valley threshold.
        if depth < false_dupe_height_thr:
            continue

        valid_itvs.append(valley)

    # Calculate stdev of entire included region.
    stdev_all = df_cov.filter(pl.col("include"))["first"].std()

    def check_in_between(itv_1: Interval, itv_2: Interval) -> bool:
        # Then get region in between potential false dupes.
        df_rgn = df_cov.filter(
            pl.col("include") & pl.col("position").is_between(itv_1.end, itv_2.begin)
        )
        mean_rgn = df_rgn["first"].mean()
        if not mean_rgn or not stdev_all:
            return False
        # Do not merge if the mean of the region is greater than the required false dupe depth + 2 stdev of the entire region.
        return mean_rgn < (false_dupe_height_thr + stdev_all * stdev_merge)

    # We merge here because the signal may not be noisy and interrupted by a small peak in coverage.
    # This would truncate the region as we require a misjoin to consider a valley a false dupe.
    valleys_potential_dupes = IntervalTree(
        merge_itvs(
            valid_itvs,
            dst=bp_merge,
            fn_cmp=check_in_between,
            fn_merge_itv=lambda i1, i2: Interval(
                i1.begin,
                i2.end,
                (
                    min(i1.data[0], i2.data[0]),
                    max(i1.data[1], i2.data[1]),
                ),
            ),
        )
    )

    for valley in valleys_potential_dupes.iter():
        misjoin_overlap = misassemblies[Misassembly.MISJOIN].overlap(valley)
        if not misjoin_overlap:
            continue
        df_below_misjoin_thr = df_cov.filter(
            pl.col("position").is_between(valley.begin, valley.end)
        ).filter(pl.col("first") < misjoin_abs_height_thr)
        prop_below_misjoin_thr = df_below_misjoin_thr.shape[0] / valley.length()
        # If majority of valley passes misjoin threshold, treat as misjoin.
        if prop_below_misjoin_thr > 0.5:
            continue

        # Remove overlapped misjoins.
        misassemblies[Misassembly.MISJOIN].remove_overlap(valley.begin, valley.end)
        new_begin, new_end = valley.begin, valley.end
        # Extend false dupe interval based on overlap.
        for itv_ovl in misjoin_overlap:
            new_begin = min(itv_ovl.begin, new_begin)
            new_end = max(itv_ovl.end, new_end)
        # Add false duplication.
        misassemblies[Misassembly.FALSE_DUP].add(Interval(new_begin, new_end))
