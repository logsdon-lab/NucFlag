#!/usr/bin/env python
import io
import os
import sys
import pprint
import tomllib
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import polars as pl
import matplotlib.pyplot as plt

from intervaltree import Interval

from .plot import plot_coverage
from .io import (
    read_ignored_regions,
    read_overlay_regions,
)
from .region import Region, add_mapq_overlay_region

from rs_nucflag import run_nucflag


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use per-base read coverage to classify/plot misassemblies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--infile",
        required=True,
        help="Input bam file or per-base coverage tsv file with 3-columns (position, first, second). If a bam file is provided, it must be indexed.",
    )
    parser.add_argument(
        "-b",
        "--input_regions",
        default=None,
        type=argparse.FileType("rt"),
        help="Bed file with regions to check.",
    )
    parser.add_argument(
        "-d",
        "--output_plot_dir",
        default=None,
        help="Output plot dir.",
    )
    parser.add_argument(
        "--output_cov_dir",
        default=None,
        help="Output coverage dir. Generates gzipped coverage bed files per region.",
    )
    parser.add_argument(
        "-o",
        "--output_misasm",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output bed file with misassembled regions.",
    )
    parser.add_argument(
        "-s",
        "--output_status",
        default=None,
        type=argparse.FileType("wt"),
        help="Bed file with status of contigs. With format: contig,start,end,misassembled|good",
    )
    parser.add_argument(
        "-t", "--threads", default=4, type=int, help="Threads for reading bam file."
    )
    parser.add_argument(
        "-p",
        "--processes",
        default=4,
        type=int,
        help="Processes for classifying/plotting.",
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        type=str,
        help="Additional threshold/params as toml file.",
    )
    parser.add_argument(
        "--ignore_regions",
        default=None,
        type=argparse.FileType("rt"),
        help="Bed file with regions to ignore. With format: contig|all,start,end,label,ignore:absolute|relative",
    )
    parser.add_argument(
        "--overlay_regions",
        nargs="*",
        type=argparse.FileType("rt"),
        help="Overlay additional regions as 4-column bedfile alongside coverage plot.",
    )
    # parser.add_argument("-v", "--version", action="version", version=version("nucflag"))
    return parser.parse_args()


def plot_misassemblies(
    itv: Interval,
    df_cov: pl.DataFrame,
    df_misasm: pl.DataFrame,
    overlay_regions: defaultdict[int, set[Region]],
    output_dir: str,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    sys.stderr.write(f"Adding mapq track for {itv.data}.\n")
    overlay_regions[0] = set(
        add_mapq_overlay_region(itv.data, df_cov.select("pos", "mapq"))
    )
    _ = plot_coverage(
        itv=itv, df_cov=df_cov, df_misasm=df_misasm, overlay_regions=overlay_regions
    )

    sys.stderr.write(f"Plotting {itv.data}.\n")

    output_plot = os.path.join(output_dir, f"{itv.data}.png")
    plt.savefig(output_plot, dpi=600, bbox_inches="tight")
    return df_cov, df_misasm


def main() -> int:
    args = parse_args()
    if args.output_plot_dir:
        os.makedirs(args.output_plot_dir, exist_ok=True)
    if args.output_cov_dir:
        os.makedirs(args.output_cov_dir, exist_ok=True)

    if isinstance(args.config, io.IOBase):
        config = tomllib.load(args.config)
        sys.stderr.write(f"Using config:\n{pprint.pformat(config)}\n")
    else:
        sys.stderr.write("Using default config.\n")
        config = args.config

    # Load ignored regions.
    if args.ignore_regions:
        ignored_regions: defaultdict[str, set[Region]] = read_ignored_regions(
            args.ignore_regions
        )
        total_ignored_positions = sum(len(v) for _, v in ignored_regions.items())
        if args.input_regions and "all" in ignored_regions:
            total_ignored_regions = len(args.input_regions.readlines())
        else:
            total_ignored_regions = len(ignored_regions)

        sys.stderr.write(
            f"Ignoring {total_ignored_positions} position(s) from {total_ignored_regions} region(s).\n"
        )
    else:
        total_ignored_positions = 0
        ignored_regions = defaultdict(set)

    # Load additional regions to overlay.
    if args.overlay_regions:
        overlay_regions: defaultdict[str, defaultdict[int, set[Region]]] = defaultdict(
            lambda: defaultdict(set)
        )
        # Pass reference of overlay regions to update.
        overlay_regions = read_overlay_regions(
            args.overlay_regions, ignored_regions=ignored_regions
        )
        added_ignored_positions = (
            sum(len(v) for _, v in ignored_regions.items()) - total_ignored_positions
        )
        sys.stderr.write(f"Overlapping {len(args.overlay_regions)} bedfile(s).\n")
        sys.stderr.write(
            f"Added {added_ignored_positions} ignored position(s) from bedfile(s).\n"
        )
    else:
        overlay_regions = defaultdict(lambda: defaultdict(set))

    results = run_nucflag(
        args.infile,
        args.input_regions.name if args.input_regions else None,
        args.threads,
        args.config,
    )
    if not args.output_plot_dir:
        return 0

    with ProcessPoolExecutor(max_workers=args.processes, max_tasks_per_child=1) as pool:
        plot_args = [
            (
                Interval(res.st, res.end, res.ctg),
                res.cov,
                res.regions,
                overlay_regions.get(res.ctg, defaultdict()),
                args.output_plot_dir,
            )
            for res in results
        ]
        futures = [
            (args[0].data, pool.submit(plot_misassemblies, *args)) for args in plot_args
        ]
        dfs: list[tuple[pl.DataFrame, pl.DataFrame]] = []
        for ctg, future in futures:
            if future.exception():
                print(f"Failed to plot {ctg} ({future.exception()})")
                continue
            dfs.append(future.result())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
