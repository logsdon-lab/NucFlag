#!/usr/bin/env python
import os
import io
import sys
import pysam
import pprint
import tomllib
import warnings
import argparse
import matplotlib

import polars as pl
import multiprocessing as mp

from collections import defaultdict
from typing import DefaultDict

from .classifier import classify_plot_assembly
from .io import read_ignored_regions, read_regions
from .region import Region, RegionStatus
from .config import DEF_CONFIG
from .constants import (
    PLOT_FONT_SIZE,
    RGX_REGION,
)

matplotlib.use("agg")
warnings.filterwarnings("ignore")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use per-base read coverage to classify/plot misassemblies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input_bam", required=True, help="Input bam file. Must be indexed."
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
        help="Bed file with status of contigs. With format: contig\tstart\tend\tmisassembled|good",
    )
    parser.add_argument(
        "-r",
        "--regions",
        nargs="*",
        help=f"Regions with the format: {RGX_REGION.pattern}",
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
        default=DEF_CONFIG,
        type=argparse.FileType("rb"),
        help="Additional threshold/params as toml file.",
    )
    parser.add_argument(
        "--ignore_regions",
        default=None,
        type=argparse.FileType("rt"),
        help="Bed file with regions to ignore. With format: contig|all\tstart\tend\tabsolute|relative",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    if args.output_plot_dir:
        os.makedirs(args.output_plot_dir, exist_ok=True)

    if isinstance(args.config, io.IOBase):
        # Fill missing config. Keep user config.
        config = DEF_CONFIG | tomllib.load(args.config)
    else:
        config = DEF_CONFIG | args.config

    sys.stderr.write(f"Using config:\n{pprint.pformat(config)}.\n")

    # Read regions in bam.
    sys.stderr.write(f"Reading in BAM file: {args.input_bam}\n")
    with pysam.AlignmentFile(args.input_bam, threads=args.threads) as bam:
        regions = list(read_regions(bam, args))

    sys.stderr.write(f"Loaded {len(regions)} region(s).\n")

    # Load ignored regions.
    ignored_regions: DefaultDict[str, list[Region]] = defaultdict(list)
    if args.ignore_regions:
        for region in read_ignored_regions(args.ignore_regions):
            ignored_regions[region.name].append(region)

        total_ignored_positions = sum(len(v) for _, v in ignored_regions.items())
        total_ignored_regions = len(
            regions if "all" in ignored_regions else ignored_regions
        )

        sys.stderr.write(
            f"Ignoring {total_ignored_positions} position(s) from {total_ignored_regions} region(s).\n"
        )

    # Set text size
    matplotlib.rcParams.update({"font.size": PLOT_FONT_SIZE})

    # for region in regions:
    #     classify_plot_assembly(
    #         args.input_bam,
    #         args.output_plot_dir,
    #         args.threads,
    #         *region,
    #         config,
    #         (
    #             ignored_regions["all"]
    #             if "all" in ignored_regions
    #             else ignored_regions.get(region[0])
    #         ),
    #     )

    with mp.Pool(processes=args.processes) as pool:
        results = pool.starmap(
            classify_plot_assembly,
            [
                (
                    args.input_bam,
                    args.output_plot_dir,
                    args.threads,
                    *region,
                    config,
                    # "all" takes precedence.
                    (
                        ignored_regions["all"]
                        if "all" in ignored_regions
                        else ignored_regions.get(region[0])
                    ),
                )
                for region in regions
            ],
        )

    # Save misassemblies to output bed
    dfs_misasm = []
    region_status = []
    for region_info, df_misasm in zip(regions, results):
        if not df_misasm.is_empty():
            region_status.append((*region_info, RegionStatus.MISASSEMBLED))
            dfs_misasm.append(df_misasm)
        else:
            region_status.append((*region_info, RegionStatus.GOOD))

    try:
        df_all_misasm: pl.DataFrame = pl.concat(dfs_misasm)
    except pl.exceptions.NoDataError:
        df_all_misasm = pl.DataFrame(schema=["contig", "start", "end"])

    df_all_misasm.sort(by=["contig", "start"]).write_csv(
        file=args.output_misasm, include_header=False, separator="\t"
    )

    if args.output_status:
        df_asm_status = pl.DataFrame(
            region_status, schema=["contig", "start", "end", "status"]
        )
        df_asm_status.sort(by=["contig", "start"]).write_csv(
            args.output_status, include_header=False, separator="\t"
        )


if __name__ == "__main__":
    raise SystemExit(main())
