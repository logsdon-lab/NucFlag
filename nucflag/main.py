#!/usr/bin/env python
import os
import sys
import time
import random
import logging
import tomllib
import tempfile
import argparse
from intervaltree import Interval
from collections import OrderedDict, defaultdict
from concurrent.futures import ProcessPoolExecutor

import polars as pl
import matplotlib.pyplot as plt

from .plot import plot_coverage
from .io import (
    read_overlay_regions,
    write_output,
    write_bigwig,
    BED9_COLS,
)
from .region import Region, add_bin_overlay_region, add_mapq_overlay_region

from py_nucflag import run_nucflag_itv, get_regions, print_config_from_preset

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

DEFAULT_WG_WINDOW = 10_000_000


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
        "--output_pileup_dir",
        default=None,
        help="Output pileup dir. Generates bigWig files per region.",
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
        default=2,
        type=int,
        help="Threads for nucflag classification.",
    )
    parser.add_argument(
        "-p",
        "--processes",
        default=8,
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
        "--add_pileup_data",
        nargs="*",
        choices=["cov", "mismatch", "mapq", "indel", "softclip"],
        default=["cov", "mismatch"],
        help="Add these pileup data types as bigWigs to --output_pileup_dir.",
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
    aln: str,
    fasta: str | None,
    itv: tuple[int, int, str],
    ignore_bed: str,
    threads: int,
    config: str | None,
    preset: str | None,
    overlay_regions: OrderedDict[str, set[Region]],
    plot_dir: str | None,
    pileup_dir: str | None,
    add_pileup_data: set[str],
    add_builtin_tracks: set[str],
) -> pl.DataFrame:
    # Safer logging. Each itv should be unique so no hash collisions?
    random.seed(hash(itv))
    wait_time = random.random()
    time.sleep(wait_time)

    # Detect misassemblies and return object with cov pileup and regions.
    res = run_nucflag_itv(
        aln,
        itv=itv,
        fasta=fasta,
        ignore_bed=ignore_bed,
        threads=threads,
        cfg=config,
        preset=preset,
    )
    st, end, ctg = itv
    ctg_coords = f"{ctg}:{st}-{end}"
    ctg_coords_filesafe = f"{ctg}_{st}-{end}"

    # Plot contig.
    if plot_dir and isinstance(res.pileup, pl.DataFrame):
        if "mapq" in add_builtin_tracks:
            logger.info(f"Adding mapq track for {ctg_coords}.")
            overlay_regions["mapq"] = set(
                add_mapq_overlay_region(ctg, res.pileup.select("pos", "mapq"))
            )

        if "bin" in add_builtin_tracks:
            logger.info(f"Adding bin track for {ctg_coords}.")
            overlay_regions["bin"] = set(
                add_bin_overlay_region(ctg, res.pileup.select("pos", "bin"))
            )

        logger.info(f"Plotting {ctg_coords}.")
        _ = plot_coverage(
            itv=Interval(st, end, ctg),
            df_pileup=res.pileup,
            df_misasm=res.regions,
            overlay_regions=overlay_regions,
        )

        output_plot = os.path.join(plot_dir, f"{ctg_coords_filesafe}.png")
        logger.info(f"Saving plot to {output_plot}.")
        plt.savefig(output_plot, dpi=600, bbox_inches="tight")

    # Output coverage.
    if pileup_dir and add_pileup_data:
        logger.info(f"Saving pileup data for {ctg_coords} to {pileup_dir}.")
        write_bigwig(
            res.pileup,
            chrom_lengths=f"{fasta}.fai",
            columns=list(add_pileup_data),
            output_dir=os.path.join(pileup_dir),
        )

    return res.regions


def main() -> int:
    args = parse_args()
    if args.output_plot_dir:
        os.makedirs(args.output_plot_dir, exist_ok=True)
    if args.output_pileup_dir:
        os.makedirs(args.output_pileup_dir, exist_ok=True)

    # Load ignored regions.
    ignore_bed = None
    tmpfile_ignore_bed = tempfile.NamedTemporaryFile("wt")
    if args.ignore_regions:
        logger.info(f"Ignoring region(s) from {args.ignore_regions}.")
        with open(args.ignore_regions, "rt") as fh:
            for line in fh:
                tmpfile_ignore_bed.write(line)
        ignore_bed = tmpfile_ignore_bed

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
            print(rgn.as_tsv(), file=tmpfile_ignore_bed)
        # Toggle ignored regions if additional ones found in overlay bed.
        if not ignore_bed and additional_ignore_regions:
            ignore_bed = tmpfile_ignore_bed
    else:
        overlay_regions = defaultdict(OrderedDict)

    if args.config:
        with open(args.config, "rb") as fh:
            cfg = tomllib.load(fh)
        cfg_general: dict = cfg.get("general", {})
        window = cfg_general.get("bp_wg_window", DEFAULT_WG_WINDOW)
    else:
        window = DEFAULT_WG_WINDOW

    # Print config to stderr.
    print_config_from_preset(args.preset, args.config)

    regions: list[tuple[int, int, str]] = get_regions(
        aln=args.infile,
        fasta=args.fasta,
        bed=args.input_regions.name if args.input_regions else None,
        window=window,
    )

    # Tracks to add that are builtin to nucflag.
    # mapq/bins/...
    added_builtin_tracks = (
        set(args.add_builtin_tracks) if args.add_builtin_tracks else set()
    )
    added_pileup_data = set(args.add_pileup_data) if args.add_pileup_data else set()
    dfs_regions: list[pl.DataFrame] = []
    logger.info(
        f"Generating outputs with {args.processes} processes with nucflag running with {args.threads} threads."
    )
    all_args = [
        [
            args.infile,
            args.fasta,
            rgn,
            ignore_bed,
            args.threads,
            args.config,
            args.preset,
            overlay_regions.get(rgn[2], OrderedDict()),
            args.output_plot_dir,
            args.output_pileup_dir,
            added_pileup_data,
            added_builtin_tracks,
        ]
        for rgn in regions
    ]
    if args.processes == 1:
        for a in all_args:
            res = plot_misassemblies(*a)
            dfs_regions.append(res)
    else:
        with ProcessPoolExecutor(
            max_workers=args.processes, max_tasks_per_child=1
        ) as pool:
            futures = [(a[2][2], pool.submit(plot_misassemblies, *a)) for a in all_args]
            for ctg, future in futures:
                if future.exception():
                    logger.error(
                        f"Failed to write output for {ctg} ({future.exception()})"
                    )
                    continue
                dfs_regions.append(future.result())

    # Remove tempfile of ignored regions.
    tmpfile_ignore_bed.close()

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
