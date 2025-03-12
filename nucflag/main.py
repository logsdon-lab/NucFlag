#!/usr/bin/env python
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
    write_output,
    BED9_COLS,
    BED_STATUS_COLS,
)
from .region import Region, add_mapq_overlay_region

from py_nucflag import run_nucflag


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use per-base read coverage to classify/plot misassemblies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--infile",
        required=True,
        help="Indexed BAM file aligned to a genome assembly.",
    )
    parser.add_argument(
        "-b",
        "--input_regions",
        default=None,
        type=argparse.FileType("rt"),
        help="BED file with regions to check.",
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
        "--output_regions",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help=f"Output bed file with checked regions. With format: {BED9_COLS}",
    )
    parser.add_argument(
        "-s",
        "--output_status",
        default=None,
        type=argparse.FileType("wt"),
        help=f"Bed file with status of contigs. With format: {BED_STATUS_COLS}",
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=4,
        type=int,
        help="Threads for nucflag classification.",
    )
    parser.add_argument(
        "-p",
        "--processes",
        default=4,
        type=int,
        help="Processes for plotting.",
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        type=str,
        help="Threshold/params as toml file.",
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
    df_cov: pl.DataFrame | None,
    df_misasm: pl.DataFrame,
    overlay_regions: defaultdict[int, set[Region]],
    plot_dir: str | None,
    cov_dir: str | None,
) -> pl.DataFrame:
    # Plot contig.
    if plot_dir and isinstance(df_cov, pl.DataFrame):
        sys.stderr.write(f"Adding mapq track for {itv.data}.\n")
        overlay_regions[0] = set(
            add_mapq_overlay_region(itv.data, df_cov.select("pos", "mapq"))
        )
        sys.stderr.write(f"Plotting {itv.data}.\n")
        _ = plot_coverage(
            itv=itv, df_cov=df_cov, df_misasm=df_misasm, overlay_regions=overlay_regions
        )

        output_plot = os.path.join(plot_dir, f"{itv.data}.png")
        sys.stderr.write(f"Saving plot to {output_plot}.\n")
        plt.savefig(output_plot, dpi=600, bbox_inches="tight")

    # Output coverage.
    if cov_dir and isinstance(df_cov, pl.DataFrame):
        sys.stderr.write(f"Saving coverage data for {itv.data}.\n")
        output_cov = os.path.join(cov_dir, f"{itv.data}.tsv.gz")
        df_cov.write_csv(output_cov, include_header=True)

    return df_misasm


def main() -> int:
    args = parse_args()
    if args.output_plot_dir:
        os.makedirs(args.output_plot_dir, exist_ok=True)
    if args.output_cov_dir:
        os.makedirs(args.output_cov_dir, exist_ok=True)

    if isinstance(args.config, str):
        with open(args.config, "rb") as fh:
            config_str = pprint.pformat(tomllib.load(fh))
        sys.stderr.write(f"Using config:\n{config_str}\n")
    else:
        sys.stderr.write("Using default config.\n")

    # Load ignored regions.
    if args.ignore_regions:
        raise NotImplementedError
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

    sys.stderr.write(f"Running nucflag with {args.threads} threads.\n")
    results = run_nucflag(
        args.infile,
        args.input_regions.name if args.input_regions else None,
        args.threads,
        args.config,
    )

    dfs_regions: list[pl.DataFrame] = []
    sys.stderr.write(f"Generating outputs with {args.processes} processes.\n")
    with ProcessPoolExecutor(max_workers=args.processes, max_tasks_per_child=1) as pool:
        plot_args = [
            (
                Interval(res.st, res.end, res.ctg),
                res.cov,
                res.regions,
                overlay_regions.get(res.ctg, defaultdict()),
                args.output_plot_dir,
                args.output_cov_dir,
            )
            for res in results
        ]
        futures = [
            (args[0].data, pool.submit(plot_misassemblies, *args)) for args in plot_args
        ]
        for ctg, future in futures:
            if future.exception():
                sys.stderr.write(f"Failed to plot {ctg} ({future.exception()})\n")
                continue
            dfs_regions.append(future.result())

    sys.stderr.write(f"Writing region BED file to {args.output_regions.name}\n")
    write_output(
        dfs_regions,
        args.output_regions,
        args.output_status,
    )
    sys.stderr.write("Done!\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
