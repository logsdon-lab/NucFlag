#!/usr/bin/env python
import os
import ast
import sys
import time
import random
import logging
import tomllib
import tempfile
import argparse
from typing import TextIO
from intervaltree import Interval  # type: ignore[import-untyped]
from collections import OrderedDict, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import polars as pl
import matplotlib.pyplot as plt

from .plot import plot_coverage
from .io import (
    STATUSES,
    read_overlay_regions,
    write_output,
    write_bigwig,
    BED9_COLS,
)
from .region import (
    Region,
    add_bin_overlay_region,
    add_mapq_overlay_region,
    add_misassemblies_overlay_region,
)

from py_nucflag import run_nucflag_itv, get_regions, print_config_from_preset  # type: ignore[import-untyped]

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
    input_args = parser.add_argument_group(title="Input", description="Input files.")
    input_args.add_argument(
        "-i",
        "--infile",
        required=True,
        help="Indexed BAM or CRAM file.",
    )
    input_args.add_argument(
        "-f",
        "--fasta",
        default=None,
        help="Reference fasta. Used to bin pileup using average nucleotide identity and detect repeats.",
    )
    input_args.add_argument(
        "-b",
        "--input_regions",
        default=None,
        type=argparse.FileType("rt"),
        help="BED file with regions to check.",
    )
    input_args.add_argument(
        "--ignore_regions",
        default=None,
        type=argparse.FileType("rt"),
        help="Bed file with regions to ignore. With format: [contig, start, end]",
    )
    output_args = parser.add_argument_group(
        title="Outputs", description="Output files."
    )
    output_args.add_argument(
        "-o",
        "--output_regions",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help=f"Output bed file with checked regions. With format: {[c[0] for c in BED9_COLS]}",
    )
    output_args.add_argument(
        "-s",
        "--output_status",
        default=None,
        type=argparse.FileType("wt"),
        help="Bed file with status of contigs and percentage breakdown of each misassembly type.",
    )
    output_args.add_argument(
        "-d",
        "--output_plot_dir",
        default=None,
        help="Output plot dir.",
    )
    output_args.add_argument(
        "--output_pileup_dir",
        default=None,
        help="Output pileup dir. Generates bigWig files per region.",
    )
    output_args.add_argument(
        "--add_pileup_data",
        nargs="*",
        choices=["cov", "mismatch", "mapq", "indel", "softclip"],
        default=["cov", "mismatch"],
        help="Add these pileup data types as bigWigs to --output_pileup_dir.",
    )
    config_args = parser.add_argument_group(
        title="Config", description="Configuration."
    )
    config_args.add_argument(
        "-t",
        "--threads",
        default=2,
        type=int,
        help="Threads for nucflag classification.",
    )
    config_args.add_argument(
        "-p",
        "--processes",
        default=8,
        type=int,
        help="Processes for plotting.",
    )
    config_args.add_argument(
        "-x",
        "--preset",
        default=None,
        choices=["ont_r9", "ont_r10", "hifi"],
        help="Sequencing data specific preset.",
    )
    config_args.add_argument(
        "-c",
        "--config",
        default=None,
        type=str,
        help="Threshold/params as toml file.",
    )
    config_args.add_argument(
        "--ignore_mtypes",
        nargs="*",
        choices=[status for status in STATUSES if status != "correct"],
        help="Ignore call types from plot and output bedfile.",
    )
    plot_args = parser.add_argument_group(title="Plot", description="Plot arguments.")
    plot_args.add_argument(
        "--overlay_regions",
        nargs="*",
        type=argparse.FileType("rt"),
        help="Overlay additional regions as BED4 or BED9 alongside coverage plot.",
    )
    plot_args.add_argument(
        "--add_builtin_tracks",
        nargs="*",
        choices=["mapq", "bin"],
        help="Add built-in tracks used in nucflag as overlay tracks.",
    )
    plot_args.add_argument(
        "--ylim",
        default=3.0,
        type=ast.literal_eval,
        help="Plot y-axis limit. If float, used as a scaling factor from mean. (ex. 3.0 is mean times 3)",
    )
    # parser.add_argument("-v", "--version", action="version", version=version("nucflag"))
    return parser.parse_args()


def plot_misassemblies(
    aln: str,
    fasta: str | None,
    itv: tuple[int, int, str],
    ignore_bed: str | None,
    threads: int,
    config: str | None,
    preset: str | None,
    overlay_regions: OrderedDict[str, set[Region]],
    plot_dir: str | None,
    pileup_dir: str | None,
    add_pileup_data: set[str],
    add_builtin_tracks: set[str],
    ignore_mtypes: list[str],
    ylim: int | float,
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

    # Ignore mtypes from results.
    if ignore_mtypes:
        # Must be sorted and from one chrom.
        # Not writable.
        regions = (
            res.regions.with_columns(
                name=pl.when(pl.col("name").is_in(ignore_mtypes))
                .then(pl.lit("correct"))
                .otherwise(pl.col("name"))
            )
            .with_columns(grp=pl.col("name").rle_id())
            .group_by(["grp"])
            .agg(
                pl.col("#chrom").first(),
                pl.col("chromStart").min(),
                pl.col("chromEnd").max(),
                pl.col("name").first(),
                # May change as result of filtering.
                pl.col("score").median(),
                pl.col("strand").first(),
                pl.col("thickStart").min(),
                pl.col("thickEnd").max(),
                pl.col("itemRgb").first(),
            )
            .drop("grp")
        )
    else:
        regions = res.regions

    # Plot contig.
    if plot_dir and isinstance(res.pileup, pl.DataFrame):
        if "mapq" in add_builtin_tracks:
            logger.info(f"Adding mapq track for {ctg_coords}.")
            overlay_regions["MAPQ"] = set(
                add_mapq_overlay_region(ctg, res.pileup.select("pos", "mapq"))
            )

        if "bin" in add_builtin_tracks:
            logger.info(f"Adding bin track for {ctg_coords}.")
            overlay_regions["Structure"] = set(
                add_bin_overlay_region(
                    ctg, res.pileup.select("pos", "bin", "bin_ident")
                )
            )

        overlay_regions["Types"] = set(add_misassemblies_overlay_region(regions))

        logger.info(f"Plotting {ctg_coords}.")
        _ = plot_coverage(
            itv=Interval(st, end, ctg),
            df_pileup=res.pileup,
            overlay_regions=overlay_regions,
            plot_ylim=ylim,
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

    return regions


def main() -> int:
    args = parse_args()
    if args.output_plot_dir:
        os.makedirs(args.output_plot_dir, exist_ok=True)
    if args.output_pileup_dir:
        os.makedirs(args.output_pileup_dir, exist_ok=True)

    if not (isinstance(args.ylim, float) or isinstance(args.ylim, int)):
        raise ValueError(f"y-axis limit must be float or int. {args.ylim}")

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
        overlay_regions: defaultdict[str, OrderedDict[str, set[Region]]] = defaultdict(
            OrderedDict
        )
        # Pass reference of overlay regions to update.
        overlay_regions = read_overlay_regions(args.overlay_regions)
        logger.info(f"Overlapping {len(args.overlay_regions)} bedfile(s).")
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
        bed=args.input_regions.name if args.input_regions else None,
        window=window,
    )

    # Ignore types.
    ignore_mtypes = []
    if args.ignore_mtypes:
        ignore_mtypes = args.ignore_mtypes
        logger.info(f"Ignoring {ignore_mtypes} from output bed and plots.")

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
    save_res = args.output_regions.name == "<stdout>" and args.output_status
    all_args = [
        [
            args.infile,
            args.fasta,
            rgn,
            ignore_bed.name if ignore_bed else None,
            args.threads,
            args.config,
            args.preset,
            overlay_regions.get(rgn[2], OrderedDict()),
            args.output_plot_dir,
            args.output_pileup_dir,
            added_pileup_data,
            added_builtin_tracks,
            ignore_mtypes,
            args.ylim,
        ]
        for rgn in regions
    ]
    if args.processes == 1:
        for a in all_args:
            res = plot_misassemblies(*a)
            # Write to file as soon as done.
            res.write_csv(args.output_regions, include_header=False, separator="\t")
            dfs_regions.append(res)
    else:
        with ProcessPoolExecutor(
            max_workers=args.processes, max_tasks_per_child=1
        ) as pool:
            futures = {pool.submit(plot_misassemblies, *a): a[2] for a in all_args}
            for future in as_completed(futures):
                if future.exception():
                    st, end, chrom = futures[future]
                    raise RuntimeError(
                        f"Failed to write output for {chrom}:{st}-{end} ({future.exception()})"
                    )
                df_res = future.result()
                # Write to file as soon as done.
                df_res.write_csv(
                    args.output_regions, include_header=False, separator="\t"
                )
                # But if writing to stdout, cannot retrieve later for status so save.
                if save_res:
                    dfs_regions.append(df_res)

    if save_res:
        output_regions: list[pl.DataFrame] | TextIO = dfs_regions
    else:
        output_regions = args.output_regions

    # Remove tempfile of ignored regions.
    tmpfile_ignore_bed.close()

    logger.info(f"Writing region BED file to {args.output_regions.name}.")
    write_output(
        output_regions,
        args.output_status,
    )
    logger.info("Done!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
