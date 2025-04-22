#!/usr/bin/env python
import os
import sys
import logging
import time
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
from .region import Region, add_bin_overlay_region, add_mapq_overlay_region

from py_nucflag import run_nucflag

# Configure logging format to match rs-nucflag
# Set UTC
logging.Formatter.converter = time.gmtime
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s.%(msecs)03dZ \033[32m%(levelname)s\033[0m  [py_nucflag::%(name)s] %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)

# Create the logger
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use per-base read coverage to classify/plot misassemblies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--infile",
        required=True,
        help="Indexed BAM or CRAM file aligned to a genome assembly.",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        default=None,
        help="Reference fasta. Used to bin pileup using average nucleotide identity. Required if --infile is a CRAM file.",
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
        "-x",
        "--preset",
        default=None,
        choices=["ont_r9", "hifi"],
        help="Sequencing data specific preset.",
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
        "--add_builtin_tracks",
        nargs="*",
        choices=["mapq", "bin"],
        help="Add built-in tracks used in nucflag as overlay tracks.",
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
    add_builtin_tracks: set[str],
) -> pl.DataFrame:
    # Plot contig.
    if plot_dir and isinstance(df_cov, pl.DataFrame):
        if "mapq" in add_builtin_tracks:
            logger.info(f"Adding mapq track for {itv.data}.")
            overlay_regions["mapq"] = set(
                add_mapq_overlay_region(itv.data, df_cov.select("pos", "mapq"))
            )

        if "bin" in add_builtin_tracks:
            logger.info(f"Adding bin track for {itv.data}.")
            overlay_regions["bin"] = set(
                add_bin_overlay_region(itv.data, df_cov.select("pos", "bin"))
            )

        logger.info(f"Plotting {itv.data}.")
        _ = plot_coverage(
            itv=itv, df_cov=df_cov, df_misasm=df_misasm, overlay_regions=overlay_regions
        )

        output_plot = os.path.join(plot_dir, f"{itv.data}.png")
        logger.info(f"Saving plot to {output_plot}.")
        plt.savefig(output_plot, dpi=600, bbox_inches="tight")

    # Output coverage.
    if cov_dir and isinstance(df_cov, pl.DataFrame):
        logger.info(f"Saving coverage data for {itv.data}.")
        output_cov = os.path.join(cov_dir, f"{itv.data}.tsv.gz")
        df_cov.write_csv(output_cov, include_header=True)

    return df_misasm


def main() -> int:
    args = parse_args()
    if args.output_plot_dir:
        os.makedirs(args.output_plot_dir, exist_ok=True)
    if args.output_cov_dir:
        os.makedirs(args.output_cov_dir, exist_ok=True)

    # Load ignored regions.
    ignored_regions_file = None
    tmpfile_ignored_regions = tempfile.NamedTemporaryFile("wt")
    if args.ignore_regions:
        logger.info(f"Ignoring region(s) from {args.ignore_regions}.")
        with open(args.ignore_regions, "rt") as fh:
            for line in fh:
                tmpfile_ignored_regions.write(line)
        ignored_regions_file = tmpfile_ignored_regions

    # Load additional regions to overlay and ignore.
    if args.overlay_regions:
        additional_ignore_regions: set[Region] = set()
        overlay_regions: defaultdict[str, OrderedDict[str, set[Region]]] = defaultdict(
            OrderedDict
        )
        # Pass reference of overlay regions to update.
        overlay_regions = read_overlay_regions(
            args.overlay_regions, ignored_regions=additional_ignore_regions
        )
        logger.info(f"Overlapping {len(args.overlay_regions)} bedfile(s).")

        # Write regions to tempfile.
        for rgn in additional_ignore_regions:
            print(rgn.as_tsv(), file=tmpfile_ignored_regions)
        # Toggle ignored regions if additional ones found in overlay bed.
        if not ignored_regions_file and additional_ignore_regions:
            ignored_regions_file = tmpfile_ignored_regions
    else:
        overlay_regions = defaultdict(OrderedDict)

    logger.info(f"Running nucflag with {args.threads} threads.")
    results = run_nucflag(
        args.infile,
        fasta=args.fasta,
        bed=args.input_regions.name if args.input_regions else None,
        ignore_bed=ignored_regions_file,
        threads=args.threads,
        cfg=args.config,
        preset=args.preset,
    )
    tmpfile_ignored_regions.close()

    # Tracks to add that are builtin to nucflag.
    # mapq/bins/...
    added_builtin_tracks = (
        set(args.add_builtin_tracks) if args.add_builtin_tracks else set()
    )

    dfs_regions: list[pl.DataFrame] = []
    logger.info(f"Generating outputs with {args.processes} processes.")
    with ProcessPoolExecutor(max_workers=args.processes, max_tasks_per_child=1) as pool:
        plot_args = [
            (
                Interval(res.st, res.end, res.ctg),
                res.cov,
                res.regions,
                overlay_regions.get(res.ctg, OrderedDict()),
                args.output_plot_dir,
                args.output_cov_dir,
                added_builtin_tracks,
            )
            for res in results
        ]
        futures = [
            (args[0].data, pool.submit(plot_misassemblies, *args)) for args in plot_args
        ]
        for ctg, future in futures:
            if future.exception():
                logger.error(f"Failed to plot {ctg} ({future.exception()})")
                continue
            dfs_regions.append(future.result())

    logger.info(f"Writing region BED file to {args.output_regions.name}.")
    write_output(
        dfs_regions,
        args.output_regions,
        args.output_status,
    )
    logger.info("Done!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
