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

    ignore_calls = set(GOOD_REGIONS)
    if args.ignore_calls:
        ignore_calls.update(args.ignore_calls)

    assert len(files) > 1, f"Only one file ({[file.name for file in files]}) provided."

    # {chrom: IntervalTree[Interval[start, end, name]]}
    itree_calls: defaultdict[str, IntervalTree] = defaultdict(IntervalTree)
    for file in files:
        df_call = pl.read_csv(
            file,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            columns=list(range(0, 4)),
            schema=dict(BED9_COLS[0:4]),
            truncate_ragged_lines=True,
        )
        df_call = df_call.filter(~pl.col("name").is_in(ignore_calls))
        for call in df_call.iter_rows(named=True):
            itree_calls[call["#chrom"]].add(
                Interval(
                    max(call["chromStart"] - args.dst, 0),
                    call["chromEnd"] + args.dst,
                    # Add filename so unique if name is same.
                    (call["name"], file.name),
                )
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
        norm = matplotlib.colors.Normalize(vmin=1, vmax=len(files))

    # Reduce overlapping intervals
    def data_reducer(
        curr: tuple[set[str], int, str], new: tuple[set[str], int, str]
    ) -> tuple[set[str], int, str]:
        return (curr[0].union(new[0]), max(curr[1], new[1]), curr[2])

    for chrom, itrees in itree_calls.items():
        itv: Interval
        ovl: set[Interval]
        all_itvs = IntervalTree()

        for itv in itrees.iter():
            # a - itv
            ovl = itrees.overlap(itv)

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
                        # (name, filename, ovl_frac)
                        mdata: tuple[str, str, str] = (*ovl_itv.data, ovl_frac)
                        valid_itv = Interval(ovl_itv.begin, ovl_itv.end, mdata)
                        new_ovl.add(valid_itv)

                ovl = new_ovl

            # Count by total file names so 2 <= valid_n <= n_files
            valid_ovls = set(oitv.data[1] for oitv in ovl)
            valid_n = len(valid_ovls) >= args.total_number_ovl
            if not valid_n:
                continue

            # Color by average overlap percent
            # Or number of overlaps
            if args.perc_ovl:
                ovl_avg = int(statistics.mean([oitv.data[2] for oitv in ovl]) * 100)
            else:
                ovl_avg = len(valid_ovls)

            ovl_names = set(oitv.data[0] for oitv in ovl)
            color = cmap(norm(ovl_avg))
            rgb = ",".join(str(int(c * 255)) for c in color[0:3])
            all_itvs.add(Interval(itv.begin, itv.end, (ovl_names, ovl_avg, rgb)))

        # Merge overlaps
        all_itvs.merge_overlaps(
            data_reducer=data_reducer,
        )
        for itv in sorted(all_itvs):
            ovl_names, ovl_avg, rgb = itv.data
            # Don't modify zero start if dst non-zero.
            st = str(0 if itv.begin == 0 else itv.begin + args.dst)
            end = str(itv.end - args.dst)
            row = (
                chrom,
                st,
                end,
                ",".join(sorted(ovl_names)),
                ".",
                str(ovl_avg),
                st,
                end,
                ",".join(str(int(c * 255)) for c in color[0:3]),
            )
            print("\t".join(row), file=args.outfile)

    return 0
