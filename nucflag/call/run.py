#!/usr/bin/env python
import os
import time
import pprint
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
    read_overlay_regions,
    write_output,
    write_bigwig,
    generate_status_from_regions,
    read_identity_breakpoints,
    BED9_COLS,
)
from .region import (
    Region,
    add_ident_overlay_region,
    add_mapq_overlay_region,
    add_misassemblies_overlay_region,
)

from py_nucflag import run_nucflag_itv, get_regions, get_config_from_preset  # type: ignore[import-untyped]

# Get the logger
logger = logging.getLogger(__name__)

DEFAULT_WG_WINDOW = 10_000_000


def plot_misassemblies(
    aln: str,
    fasta: str | None,
    itv: tuple[int, int, str],
    ignore_bed: str | None,
    threads: int,
    config: str | None,
    preset: str | None,
    tracks: OrderedDict[str, set[Region]],
    ovl_tracks: OrderedDict[str, set[Region]],
    plot_dir: str | None,
    pileup_dir: str | None,
    add_pileup_data: set[str],
    add_builtin_tracks: set[str],
    ignore_mtypes: list[str],
    overlap_calls: bool,
    ylim: int | float,
    ident_breakpoints: tuple[list[float], list[str]],
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
        if "ident" in add_builtin_tracks:
            logger.info(f"Adding local self-identity track for {ctg_coords}.")
            tracks["Self-identity"] = set(
                add_ident_overlay_region(
                    ctg, res.pileup.select("pos", "bin", "bin_ident"), ident_breakpoints
                )
            )

        if "mapq" in add_builtin_tracks:
            logger.info(f"Adding mapq track for {ctg_coords}.")
            tracks["MAPQ"] = set(
                add_mapq_overlay_region(ctg, res.pileup.select("pos", "mapq"))
            )

        if overlap_calls:
            ovl_tracks["Calls"] = set(add_misassemblies_overlay_region(regions))
        else:
            tracks["Calls"] = set(add_misassemblies_overlay_region(regions))

        logger.info(f"Plotting {ctg_coords}.")
        fig, _ = plot_coverage(
            itv=Interval(st, end, ctg),
            df_pileup=res.pileup,
            tracks=tracks,
            ovl_tracks=ovl_tracks,
            plot_ylim=ylim,
        )

        output_plot = os.path.join(plot_dir, f"{ctg_coords_filesafe}.png")
        logger.info(f"Saving plot to {output_plot}.")
        plt.savefig(output_plot, dpi=600, bbox_inches="tight")
        plt.close(fig)

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


def call_misassemblies(args: argparse.Namespace) -> int:
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
    if args.tracks:
        tracks: defaultdict[str, OrderedDict[str, set[Region]]] = defaultdict(
            OrderedDict
        )
        tracks = read_overlay_regions(args.tracks)
        logger.info(f"Adding {len(args.tracks)} bedfile(s).")
    else:
        tracks = defaultdict(OrderedDict)

    if args.overlap_tracks:
        ovl_tracks: defaultdict[str, OrderedDict[str, set[Region]]] = defaultdict(
            OrderedDict
        )
        ovl_tracks = read_overlay_regions(args.overlap_tracks)
        logger.info(f"Overlapping {len(args.overlap_tracks)} bedfile(s).")
    else:
        ovl_tracks = defaultdict(OrderedDict)

    if args.config:
        with open(args.config, "rb") as fh:
            cfg = tomllib.load(fh)
        cfg_general: dict = cfg.get("general", {})
        window = cfg_general.get("bp_wg_window", DEFAULT_WG_WINDOW)
    else:
        window = DEFAULT_WG_WINDOW

    # Print config to stderr.
    config_str = get_config_from_preset(args.preset, args.config)
    cfg = tomllib.loads(config_str)
    logger.info(f"Using config:\n{pprint.pformat(cfg, underscore_numbers=True)}")

    regions: list[tuple[int, int, str]] = get_regions(
        aln=args.infile,
        bed=args.input_regions.name if args.input_regions else None,
        window=window,
    )
    # Load identity breakpoints
    ident_breakpoints = read_identity_breakpoints(args.ident_breakpoints)

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
            tracks.get(rgn[2], OrderedDict()),
            ovl_tracks.get(rgn[2], OrderedDict()),
            args.output_plot_dir,
            args.output_pileup_dir,
            added_pileup_data,
            added_builtin_tracks,
            ignore_mtypes,
            args.overlap_calls,
            args.ylim,
            ident_breakpoints,
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


def create_status(args: argparse.Namespace) -> int:
    df_regions = pl.read_csv(
        args.infile,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        schema=dict(BED9_COLS),
    )
    if df_regions.is_empty():
        raise ValueError(f"No regions to generate status for {args.infile}")

    df_status = generate_status_from_regions(df_regions)
    df_status.write_csv(file=args.outfile, include_header=True, separator="\t")
    return 0
