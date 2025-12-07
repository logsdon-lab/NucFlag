import argparse
import statistics

import matplotlib.colors

import polars as pl
import matplotlib.pyplot as plt

from typing import TextIO, Callable
from collections import defaultdict

from intervaltree import Interval, IntervalTree

from .constants import GOOD_REGIONS
from ..common import BED9_COLS


def get_consensus_calls(args: argparse.Namespace) -> int:
    files: list[TextIO] = args.infiles

    filter_calls = set(GOOD_REGIONS)
    if args.filter_calls:
        filter_calls.update(args.filter_calls)

    assert len(files) > 1, f"Only one file ({[file.name for file in files]}) provided."

    # {chrom: {file: IntervalTree[Interval[start, end, name]]}}
    itree_calls: defaultdict[str, defaultdict[int, IntervalTree]] = defaultdict(
        lambda: defaultdict(IntervalTree)
    )
    for i, file in enumerate(files):
        df_call = pl.read_csv(
            file,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            columns=list(range(0, 4)),
            schema=dict(BED9_COLS[0:4]),
            truncate_ragged_lines=True,
        )
        df_call = df_call.filter(~pl.col("name").is_in(filter_calls))
        for call in df_call.iter_rows(named=True):
            itree_calls[call["#chrom"]][i].add(
                Interval(call["chromStart"], call["chromEnd"], call["name"])
            )

    if args.perc_ovl:
        assert args.perc_ovl < 1.0, f"Percent overlap ({args.perc_ovl}) exceeds 100%"

    fn_ovl_check: Callable[[float, float], tuple[bool, float]] | None
    if args.perc_ovl and args.type_ovl == "A":

        def fn_ovl_check(a_b_ovl_frac, b_a_ovl_frac):
            return (b_a_ovl_frac > args.perc_ovl, b_a_ovl_frac)

    elif args.perc_ovl and args.type_ovl == "B":

        def fn_ovl_check(a_b_ovl_frac, b_a_ovl_frac):
            return (
                a_b_ovl_frac > args.perc_ovl,
                a_b_ovl_frac,
            )

    elif args.perc_ovl and args.type_ovl == "A_B":

        def fn_ovl_check(a_b_ovl_frac, b_a_ovl_frac):
            return (
                a_b_ovl_frac > args.perc_ovl and b_a_ovl_frac > args.perc_ovl,
                (a_b_ovl_frac + b_a_ovl_frac) / 2.0,
            )
    else:
        fn_ovl_check = None

    # https://stackoverflow.com/a/43009963
    cmap = plt.cm.get_cmap("rainbow")
    if args.perc_ovl:
        norm = matplotlib.colors.Normalize(vmin=0, vmax=100)
    else:
        norm = matplotlib.colors.Normalize(vmin=1, vmax=10)

    # Determine order based on number of calls or chosen index.
    for chrom, file_itrees in itree_calls.items():
        if args.relative_to:
            idx_file_w_max_calls = args.relative_to
        else:
            idx_file_w_max_calls = max(
                file_itrees.items(), key=lambda itree: len(itree[1])
            )[0]

        if not file_itrees.get(idx_file_w_max_calls):
            raise ValueError(
                f"Invalid index ({idx_file_w_max_calls}) for number of files: {len(file_itrees)}"
            )

        idx_other_calls = set(file_itrees.keys()).difference([idx_file_w_max_calls])

        # Store name + idx of callset in interval
        all_other_calls = IntervalTree(
            Interval(itv.begin, itv.end, (itv.data, idx))
            for idx in idx_other_calls
            for itv in file_itrees[idx]
        )
        itv: Interval
        ovl: set[Interval]
        for itv in file_itrees[idx_file_w_max_calls].iter():
            # a - itv
            ovl = all_other_calls.overlap(itv)

            # Check overlaps to ensure meets overlap criteria.
            if ovl and args.perc_ovl and fn_ovl_check:
                itv_len = itv.length()
                new_ovl = set()
                for ovl_itv in ovl:
                    # b - itv_ovl
                    ovl_itv_len = ovl_itv.length()
                    a_b_ovl_frac = itv_len / ovl_itv_len
                    b_a_ovl_frac = ovl_itv_len / itv_len
                    is_valid, ovl_frac = fn_ovl_check(a_b_ovl_frac, b_a_ovl_frac)
                    if is_valid:
                        # (name, file_index, ovl_frac)
                        mdata: tuple[str, int, str] = (*ovl_itv.data, ovl_frac)
                        valid_itv = Interval(ovl_itv.begin, ovl_itv.end, mdata)
                        new_ovl.add(valid_itv)

                ovl = new_ovl

            valid_n = len(ovl) + 1 >= args.number_ovl
            if not valid_n:
                continue

            ovl_names = f"{itv.data},{','.join(oitv.data[0] for oitv in ovl)}"
            # Color by average overlap percent
            # Or number of overlaps
            if args.perc_ovl:
                ovl_avg = int(statistics.mean([oitv.data[2] for oitv in ovl]) * 100)
            else:
                ovl_avg = len(ovl) + 1

            color = cmap(norm(ovl_avg))
            row = (
                chrom,
                str(itv.begin),
                str(itv.end),
                ovl_names,
                ".",
                str(ovl_avg),
                str(itv.begin),
                str(itv.end),
                ",".join(str(int(c * 255)) for c in color[0:3]),
            )
            print("\t".join(row), file=args.outfile)

    return 0
