#!/usr/bin/env python
import re
import os
import io
import sys
import pysam
import pprint
import tomllib
import warnings
import argparse
import matplotlib
import scipy.signal

import portion as pt
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import multiprocessing as mp

from enum import StrEnum
from typing import Generator, Any

matplotlib.use("agg")
warnings.filterwarnings("ignore")

PLOT_FONT_SIZE = 16
PLOT_HEIGHT = 6
PLOT_WIDTH = 16
PLOT_DPI = 600
PLOT_YLIM = 100
RGX_REGION = re.compile(r"(.+):(\d+)-(\d+)")

DEF_CONFIG = {
    "first": dict(
        added_region_bounds=0,
        thr_min_peak_horizontal_distance=100_000,
        thr_min_peak_width=20,
        thr_min_valley_horizontal_distance=100_000,
        thr_min_valley_width=10,
        thr_peak_height_std_above=4,
        thr_valley_height_std_below=3,
    ),
    "second": dict(
        thr_min_perc_first=0.1,
        thr_peak_height_std_above=3,
        group_distance=30_000,
        thr_min_group_size=5,
        thr_min_group_len=500,
        thr_collapse_het_ratio=0.1,
    ),
}


class Misassembly(StrEnum):
    COLLAPSE_VAR = "COLLAPSE_VAR"
    COLLAPSE = "COLLAPSE"
    MISJOIN = "MISJOIN"
    GAP = "GAP"
    FALSE_DUP = "FALSE_DUP"

    def as_color(self) -> str:
        match self:
            case self.COLLAPSE_VAR:
                return "blue"
            case self.COLLAPSE:
                return "green"
            case self.MISJOIN:
                return "yellow"
            case self.GAP:
                return "gray"
            case self.FALSE_DUP:
                return "purple"
            case _:
                raise ValueError(f"Invalid color {self}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use per-base read coverage to classify/plot misassemblies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input_bam", help="Input bam file. Must be indexed.")
    parser.add_argument(
        "-b", "--input_bed", default=None, help="Bed file with regions to plot."
    )
    parser.add_argument(
        "-d",
        "--output_plot_dir",
        default=None,
        help="Output plot dir.",
    )
    parser.add_argument(
        "-o",
        "--output_bed",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output bed file with misassembled regions.",
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
        help="Addtional threshold/params as toml file.",
    )

    return parser.parse_args()


def get_coverage_by_base(
    bam: pysam.AlignmentFile, contig: str, start: int, end: int
) -> np.ndarray:
    coverage = bam.count_coverage(
        contig, start=start, stop=end, read_callback="nofilter", quality_threshold=None
    )
    assert len(coverage) == 4
    return np.array(tuple(zip(*coverage)))


def read_regions(
    bam: pysam.AlignmentFile, args: argparse.Namespace
) -> Generator[tuple[str, int, int], None, None]:
    refs = {}

    if args.regions is not None or args.input_bed is not None:
        if args.regions is not None:
            sys.stderr.write(f"Reading in {len(args.regions)} region(s).\n")

            for region in args.regions:
                match = re.match(RGX_REGION, region)
                assert match, region + " not valid!"
                chrm, start, end = match.groups()
                refs[chrm] = [int(start), int(end)]
                yield (chrm, int(start), int(end))

        if args.input_bed is not None:
            sys.stderr.write(f"Reading in regions from {args.input_bed}.\n")

            for line in open(args.input_bed):
                if line[0] == "#":
                    continue
                chrm, start, end, *_ = line.strip().split()
                refs[chrm] = [int(start), int(end)]
                yield (chrm, int(start), int(end))

    else:
        sys.stderr.write(
            "Reading the whole bam because no region or bed argument was made.\n"
        )
        for read in bam.fetch(until_eof=True):
            ref = read.reference_name
            if ref not in refs:
                refs[ref] = [2147483648, 0]

            start = read.reference_start
            end = read.reference_end
            if refs[ref][0] > start:
                refs[ref][0] = start
            if refs[ref][1] < end:
                refs[ref][1] = end

        for contig in refs:
            yield (contig, refs[contig][0], refs[contig][1])


def plot_coverage(
    df: pl.DataFrame,
    misassemblies: dict[Misassembly, set[pt.Interval]],
    contig_name: str,
) -> tuple[plt.Figure, Any]:
    fig, ax = plt.subplots(figsize=(PLOT_WIDTH, PLOT_HEIGHT))

    (_,) = ax.plot(
        df["position"],
        df["first"],
        "o",
        color="black",
        markeredgewidth=0.0,
        markersize=2,
        label="Most Frequent Base",
    )
    (_,) = ax.plot(
        df["position"],
        df["second"],
        "o",
        color="red",
        markeredgewidth=0.0,
        markersize=2,
        label="Second Most Frequent Base",
    )

    # Add misassembly rect patches to highlight region.
    for misasm, regions in misassemblies.items():
        color = misasm.as_color()
        for region in regions:
            plt.axvspan(
                region.lower, region.upper, color=color, alpha=0.4, label=misasm
            )

    # Add legend. Deduplicate multiple labels.
    # https://stackoverflow.com/a/36189073
    handles, labels = plt.gca().get_legend_handles_labels()
    labels, ids = np.unique(labels, return_index=True)
    handles = [handles[i] for i in ids]
    plt.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.2),
        ncols=len(labels),
        borderaxespad=0,
        fancybox=True,
    )

    maxval = df["position"].max()
    minval = df["position"].min()
    subval = 0

    title = "{}:{}-{}\n".format(contig_name, minval, maxval)
    ax.set_title(title, fontweight="bold")

    if maxval < 1_000_000:
        xlabels = [format((label - subval), ",.0f") for label in ax.get_xticks()]
        lab = "bp"
    elif maxval < 10_000_000:
        xlabels = [format((label - subval) / 1000, ",.1f") for label in ax.get_xticks()]
        lab = "kbp"
    else:
        xlabels = [format((label - subval) / 1000, ",.1f") for label in ax.get_xticks()]
        lab = "kbp"

    ax.set_ylim(0, PLOT_YLIM)
    ax.set_xlabel("Assembly position ({})".format(lab), fontweight="bold")
    ax.set_ylabel("Sequence read depth", fontweight="bold")
    ax.set_xticklabels(xlabels)

    # Hide the right and top spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    fig.tight_layout()

    return fig, ax


def peak_finder(
    data: np.ndarray,
    positions: np.ndarray,
    *,
    height: int,
    added_region_bounds: int,
    distance: int,
    width: int,
) -> list[pt.Interval]:
    _, peak_info = scipy.signal.find_peaks(
        data, height=height, distance=distance, width=width
    )
    return [
        pt.open(
            positions[int(left_pos)] - added_region_bounds,
            positions[int(right_pos)] + added_region_bounds,
        )
        for left_pos, right_pos in zip(peak_info["left_ips"], peak_info["right_ips"])
    ]


# https://stackoverflow.com/a/7353335
def consecutive(data, stepsize: int = 1):
    return np.split(data, np.where((np.diff(data) <= stepsize) == False)[0] + 1)  # noqa: E712


def classify_misassemblies(
    cov_first_second: np.ndarray,
    positions: np.ndarray,
    *,
    config: dict[str, Any],
) -> tuple[pl.DataFrame, dict[Misassembly, set[pt.Interval]]]:
    df = pl.DataFrame(
        {
            "position": positions,
            "first": cov_first_second[0],
            "second": cov_first_second[1],
        }
    )
    del cov_first_second

    # Calculate std and mean for both most and second most freq read.
    # Remove gaps which would artificially lower mean.
    df_gapless = df.filter(pl.col("first") != 0)
    mean_first, stdev_first = df_gapless["first"].mean(), df_gapless["first"].std()
    mean_second, stdev_second = df_gapless["second"].mean(), df_gapless["second"].std()
    del df_gapless

    first_peak_coords = peak_finder(
        df["first"],
        positions,
        height=mean_first
        + (config["first"]["thr_peak_height_std_above"] * stdev_first),
        added_region_bounds=config["first"]["added_region_bounds"],
        distance=config["first"]["thr_min_peak_horizontal_distance"],
        width=config["first"]["thr_min_peak_width"],
    )
    first_valley_coords = peak_finder(
        -df["first"],
        positions,
        height=-(
            mean_first - (config["first"]["thr_valley_height_std_below"] * stdev_first)
        ),
        added_region_bounds=config["first"]["added_region_bounds"],
        distance=config["first"]["thr_min_valley_horizontal_distance"],
        width=config["first"]["thr_min_valley_width"],
    )

    # Remove secondary rows that don't meet minimal secondary coverage.
    second_thr = max(
        round(mean_first * config["second"]["thr_min_perc_first"]),
        round(
            mean_second + (config["second"]["thr_peak_height_std_above"] * stdev_second)
        ),
    )

    classified_second_outliers = set()
    df_second_outliers = df.filter(pl.col("second") > second_thr)
    # Group consecutive positions allowing a maximum gap of stepsize.
    # Larger stepsize groups more positions.
    second_outliers_coords = []
    for grp in consecutive(
        df_second_outliers["position"], stepsize=config["second"]["group_distance"]
    ):
        if len(grp) < config["second"]["thr_min_group_size"]:
            continue
        coords = pt.open(grp[0], grp[-1])
        coords_len = coords.upper - coords.lower
        if coords_len > config["second"]["thr_min_group_len"]:
            second_outliers_coords.append(coords)

    misassemblies: dict[Misassembly, set[pt.Interval]] = {m: set() for m in Misassembly}

    # Intersect intervals and classify collapses.
    for peak in first_peak_coords:
        for second_outlier in second_outliers_coords:
            if second_outlier in peak:
                misassemblies[Misassembly.COLLAPSE_VAR].add(peak)
                classified_second_outliers.add(second_outlier)

        if peak not in misassemblies[Misassembly.COLLAPSE_VAR]:
            misassemblies[Misassembly.COLLAPSE].add(peak)

    # Classify gaps.
    df_gaps = df.filter(pl.col("first") == 0)
    misassemblies[Misassembly.GAP] = set(
        pt.open(grp[0], grp[-1])
        for grp in consecutive(df_gaps["position"], stepsize=1)
        if len(grp) > 1
    )

    # Classify misjoins.
    for valley in first_valley_coords:
        for second_outlier in second_outliers_coords:
            if second_outlier in valley:
                misassemblies[Misassembly.MISJOIN].add(valley)
                classified_second_outliers.add(second_outlier)

    # Check remaining secondary regions not categorized.
    for second_outlier in second_outliers_coords:
        if second_outlier in classified_second_outliers:
            continue

        df_second_outlier = df.filter(
            (pl.col("position") >= second_outlier.lower)
            & (pl.col("position") <= second_outlier.upper)
            & (pl.col("second") != 0)
        )
        df_second_outlier_het_ratio = df_second_outlier.mean().with_columns(
            het_ratio=pl.col("second") / (pl.col("first") + pl.col("second"))
        )
        # Use het ratio to classfiy.
        if (
            df_second_outlier_het_ratio["het_ratio"][0]
            > config["second"]["thr_collapse_het_ratio"]
        ):
            misassemblies[Misassembly.COLLAPSE_VAR].add(second_outlier)
        else:
            misassemblies[Misassembly.MISJOIN].add(second_outlier)

    # Annotate df with misassembliy.
    df = df.with_columns(status=pl.lit("Good"))
    for mtype, regions in misassemblies.items():
        for region in regions:
            df = df.with_columns(
                status=pl.when(
                    (pl.col("position") >= region.lower)
                    & (pl.col("position") <= region.upper)
                )
                .then(pl.lit(mtype))
                .otherwise(pl.col("status"))
            )

    # TODO: false dupes
    return df, misassemblies


def classify_plot_assembly(
    infile: str,
    output_dir: str | None,
    threads: int,
    contig: str,
    start: int,
    end: int,
    config: dict[str, Any],
) -> pl.DataFrame:
    bam = pysam.AlignmentFile(infile, threads=threads)
    contig_name = f"{contig}:{start}-{end}"

    sys.stderr.write(f"Reading in NucFreq from region: {contig_name}\n")

    df_group_labeled, miassemblies = classify_misassemblies(
        np.flip(
            np.sort(get_coverage_by_base(bam, contig, start, end), axis=1)
        ).transpose(),
        np.arange(start, end),
        config=config,
    )

    if output_dir:
        _ = plot_coverage(df_group_labeled, miassemblies, contig)

        sys.stderr.write(f"Plotted {contig_name}.\n")

        output_plot = os.path.join(output_dir, f"{contig_name}.png")
        plt.tight_layout()
        plt.savefig(output_plot, dpi=PLOT_DPI)

    df_misassemblies = pl.DataFrame(
        [
            (contig_name, interval.lower, interval.upper, misasm)
            for misasm, intervals in miassemblies.items()
            for interval in intervals
        ],
        schema=["contig", "start", "stop", "misassembly"],
    )
    return df_misassemblies


def main():
    args = parse_args()
    if args.output_plot_dir:
        os.makedirs(args.output_plot_dir, exist_ok=True)

    if isinstance(args.config, io.IOBase):
        # Fill missing config. Keep user config.
        config = DEF_CONFIG | tomllib.load(args.config)
    else:
        config = args.config

    sys.stderr.write(f"Using config:\n{pprint.pformat(config)}.\n")

    # Read regions and close file handle to bam.
    bam = pysam.AlignmentFile(args.input_bam, threads=args.threads)
    sys.stderr.write(f"Reading in BAM file: {args.input_bam}\n")

    regions = list(read_regions(bam, args))
    sys.stderr.write(f"Loaded {len(regions)} region(s).\n")

    bam.close()

    # Set text size
    matplotlib.rcParams.update({"font.size": PLOT_FONT_SIZE})

    with mp.Pool(processes=args.processes) as pool:
        results = pool.starmap(
            classify_plot_assembly,
            [
                (args.input_bam, args.output_plot_dir, args.threads, *region, config)
                for region in regions
            ],
        )

    # Save misassemblies to output bed
    all_misasm = pl.concat(result for result in results if not result.is_empty()).sort(
        by=["contig", "start"]
    )
    all_misasm.write_csv(file=args.output_bed, include_header=False, separator="\t")


if __name__ == "__main__":
    raise SystemExit(main())