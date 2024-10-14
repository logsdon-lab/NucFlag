#!/usr/bin/env python
import io
import os
import sys
import pprint
import argparse
from collections import defaultdict
from typing import DefaultDict
from concurrent.futures.process import ProcessPoolExecutor
from importlib.metadata import version

import tomllib

from .classifier.classifier import classify_plot_assembly
from .config import DEF_CONFIG
from .io import (
    read_asm_regions,
    read_ignored_regions,
    read_overlay_regions,
    write_misassemblies_and_status,
)
from .region import Region


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
        default=DEF_CONFIG,
        type=argparse.FileType("rb"),
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
    parser.add_argument("-v", "--version", action="version", version=version("nucflag"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.output_plot_dir:
        os.makedirs(args.output_plot_dir, exist_ok=True)
    if args.output_cov_dir:
        os.makedirs(args.output_cov_dir, exist_ok=True)

    if isinstance(args.config, io.IOBase):
        config = tomllib.load(args.config)
    else:
        config = args.config

    sys.stderr.write(f"Using config:\n{pprint.pformat(config)}.\n")

    # Read regions in bam or from input regions.
    general_cfg: dict = config.get("general", {})
    regions = list(
        read_asm_regions(
            args.infile,
            args.input_regions,
            threads=args.threads,
            window_size=general_cfg.get(
                "window_size", DEF_CONFIG["general"]["window_size"]
            ),
        )
    )
    sys.stderr.write(f"Loaded {len(regions)} region(s).\n")

    # Load ignored regions.
    if args.ignore_regions:
        ignored_regions: DefaultDict[str, set[Region]] = read_ignored_regions(
            args.ignore_regions
        )
        total_ignored_positions = sum(len(v) for _, v in ignored_regions.items())
        total_ignored_regions = len(
            regions if "all" in ignored_regions else ignored_regions
        )
        sys.stderr.write(
            f"Ignoring {total_ignored_positions} position(s) from {total_ignored_regions} region(s).\n"
        )
    else:
        total_ignored_positions = 0
        ignored_regions = defaultdict(set)

    # Load additional regions to overlay.
    if args.overlay_regions:
        overlay_regions: DefaultDict[str, DefaultDict[int, set[Region]]] = defaultdict(
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

    # res = []
    # for region in regions:
    #     res.append(classify_plot_assembly(
    #         args.infile,
    #         args.output_plot_dir,
    #         args.output_cov_dir,
    #         args.threads,
    #         *region,
    #         config,
    #         overlay_regions.get(region[0]),
    #         ignored_regions.get("all", set()).union(ignored_regions.get(region[0], set())),
    #     ))

    # Use new process pool, which doesn't cause a memory leak.
    # Also create new processes aggressively to clear resources.
    # https://stackoverflow.com/a/61492363
    with ProcessPoolExecutor(max_workers=args.processes, max_tasks_per_child=1) as pool:
        res = pool.map(
            classify_plot_assembly,
            *zip(
                *[
                    (
                        args.infile,
                        args.output_plot_dir,
                        args.output_cov_dir,
                        args.threads,
                        *region,
                        config,
                        overlay_regions.get(region[0]),
                        ignored_regions.get("all", set()).union(
                            ignored_regions.get(region[0], set())
                        ),
                    )
                    for region in regions
                ]
            ),
        )

    write_misassemblies_and_status(
        res,
        regions,
        args.output_misasm,
        args.output_status,
        bed_provided=bool(args.input_regions),
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
