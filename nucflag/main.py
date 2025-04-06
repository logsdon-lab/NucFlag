#!/usr/bin/env python
import os
import sys
import pprint
import tomllib
import tempfile
import argparse
from collections import OrderedDict, defaultdict
from concurrent.futures import ProcessPoolExecutor

import polars as pl
import matplotlib.pyplot as plt

from intervaltree import Interval

from .plot import plot_coverage
from .io import (
    read_overlay_regions,
    write_output,
    BED9_COLS,
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
        help="Bed file with status of contigs and percentage breakdown of each misassembly type.",
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
    parser.add_argument(
        "--add_mapq",
        action="store_true",
        help="Add mapq as an overlay track.",
    )
    # parser.add_argument("-v", "--version", action="version", version=version("nucflag"))
    return parser.parse_args()


def plot_misassemblies(
    itv: Interval,
    df_cov: pl.DataFrame | None,
    df_misasm: pl.DataFrame,
    overlay_regions: OrderedDict[str, set[Region]],
    plot_dir: str | None,
    cov_dir: str | None,
    add_mapq: bool,
) -> pl.DataFrame:
    # Plot contig.
    if plot_dir and isinstance(df_cov, pl.DataFrame):
        if add_mapq:
            sys.stderr.write(f"Adding mapq track for {itv.data}.\n")
            overlay_regions["mapq"] = set(
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
    ignored_regions_files: list[str] = []
    if args.ignore_regions:
        sys.stderr.write(f"Ignoring region(s) from {args.ignore_regions}.\n")
        ignored_regions_files.append(args.ignore_regions)

    # Load additional regions to overlay and ignore.
    tmpfile_ignored_regions = tempfile.NamedTemporaryFile("wt")
    if args.overlay_regions:
        additional_ignore_regions: set[Region] = set()
        overlay_regions: defaultdict[str, OrderedDict[str, set[Region]]] = defaultdict(
            OrderedDict
        )
        # Pass reference of overlay regions to update.
        overlay_regions = read_overlay_regions(
            args.overlay_regions, ignored_regions=additional_ignore_regions
        )
        sys.stderr.write(f"Overlapping {len(args.overlay_regions)} bedfile(s).\n")

        # Write regions to tempfile.
        for rgn in additional_ignore_regions:
            print(rgn.as_tsv(), file=tmpfile_ignored_regions)
    else:
        overlay_regions = defaultdict(OrderedDict)

    sys.stderr.write(f"Running nucflag with {args.threads} threads.\n")
    results = run_nucflag(
        args.infile,
        args.input_regions.name if args.input_regions else None,
        args.threads,
        args.config,
    )
    tmpfile_ignored_regions.close()

    dfs_regions: list[pl.DataFrame] = []
    sys.stderr.write(f"Generating outputs with {args.processes} processes.\n")
    with ProcessPoolExecutor(max_workers=args.processes, max_tasks_per_child=1) as pool:
        plot_args = [
            (
                Interval(res.st, res.end, res.ctg),
                res.cov,
                res.regions,
                overlay_regions.get(res.ctg, OrderedDict()),
                args.output_plot_dir,
                args.output_cov_dir,
                args.add_mapq,
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
