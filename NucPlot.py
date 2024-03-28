#!/usr/bin/env python
import re
import os
import sys
import pysam
import warnings
import argparse
import matplotlib
import scipy.signal

matplotlib.use("agg")
import portion as pt
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import multiprocessing as mp

from enum import StrEnum, auto
from typing import Generator, Any

warnings.filterwarnings("ignore")

PLOT_COLORS = sns.color_palette()
PLOT_FONT_SIZE = 16
PLOT_REGION_HEIGHT = 5
PLOT_WIDTH = 16
PLOT_DPI = 600


class Misassembly(StrEnum):
    COLLAPSE_VAR = auto()
    COLLAPSE = auto()
    MISJOIN = auto()
    GAP = auto()
    FALSE_DUP = auto()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use per-base read coverage to classify/plot misassemblies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input_bam", help="Input bam file.")
    parser.add_argument(
        "-b", "--input_bed", default=None, help="Bed file with regions to plot."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        default=None,
        help="Output plot dir.",
    )
    parser.add_argument(
        "-m",
        "--output_bed",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output bed file with misassembled regions.",
    )
    parser.add_argument(
        "-r", "--regions", nargs="*", help="Regions with the format: (.*):(\\d+)-(\\d+)"
    )
    parser.add_argument(
        "--input_repeatmasker",
        help="Input RepeatMasker output file to add to plot.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "-t", "--threads", default=4, help="Threads for reading bam file."
    )
    parser.add_argument(
        "-p", "--processes", default=4, help="Processes for classifying/plotting."
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
        sys.stderr.write("Reading in the region or bed argument(s).\n")
        if args.regions is not None:
            for region in args.regions:
                match = re.match(r"(.+):(\d+)-(\d+)", region)
                assert match, region + " not valid!"
                chrm, start, end = match.groups()
                refs[chrm] = [int(start), int(end)]
                yield (chrm, int(start), int(end))

        if args.input_bed is not None:
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


def read_repeatmasker(input_rm: str) -> pl.DataFrame | None:
    rm_output = None
    if input_rm is not None:
        names = [
            "score",
            "perdiv",
            "perdel",
            "perins",
            "qname",
            "start",
            "end",
            "left",
            "strand",
            "repeat",
            "family",
            "rstart",
            "rend",
            "rleft",
            "ID",
        ]

        rm_output = (
            pl.scan_csv(
                input_rm,
                skip_rows=3,
                has_header=False,
                new_columns=names,
                dtypes={"start": pl.Int64, "end": pl.Int64},
            )
            .with_columns(label=pl.col("family").str.replace("/.*", ""))
            .collect()
        )

        cmap = {
            lab: PLOT_COLORS[idx % len(PLOT_COLORS)]
            for idx, lab in enumerate(sorted(rm_output["label"].unique()))
        }

        rm_output = rm_output.with_columns(color=pl.col("label").replace(cmap))

    return rm_output


def plot_coverage(
    df: pl.DataFrame,
    contig_name: str,
    rm_output: pl.DataFrame | None,
) -> tuple[plt.Figure, Any]:
    ylim = 100

    # get the correct axis
    fig, ax = plt.subplots(figsize=(PLOT_WIDTH, PLOT_REGION_HEIGHT))

    if rm_output:
        rmax = ax
        sys.stderr.write("Subsetting the repeatmasker file.\n")
        rm = rm_output.filter(
            (pl.col("qname") == contig_name)
            & (pl.col("start") >= df["position"].min())
            & (pl.col("end") <= df["position"].max())
        )
        assert len(rm.shape[0]) != 0, "No matching RM contig"

        rmlength = len(rm.shape[0]) * 1.0
        height_offset = ylim / 20
        for rmcount, row in enumerate(rm.iter_rows(named=True)):
            sys.stderr.write(
                "\rDrawing the {} repeatmasker rectangles:\t{:.2%}".format(
                    rmlength, rmcount / rmlength
                )
            )
            width = row["end"] - row["start"]
            rect = patches.Rectangle(
                (row["start"], ylim - height_offset),
                width,
                height_offset,
                linewidth=1,
                edgecolor="none",
                facecolor=row["color"],
                alpha=0.75,
            )
            rmax.add_patch(rect)

        sys.stderr.write("\nPlotting the repeatmasker rectangles.\n")
        plt.show()
        sys.stderr.write("Done plotting the repeatmasker rectangles.\n")

    (prime,) = ax.plot(
        df["position"],
        df["first"],
        "o",
        color="black",
        markeredgewidth=0.0,
        markersize=2,
        label="most frequent base pair",
    )
    (sec,) = ax.plot(
        df["position"],
        df["second"],
        "o",
        color="red",
        markeredgewidth=0.0,
        markersize=2,
        label="second most frequent base pair",
    )

    maxval = df["position"].max()
    minval = df["position"].min()
    subval = 0

    title = "{}:{}-{}\n".format(contig_name, minval, maxval)
    ax.set_title(title, fontweight="bold")

    if maxval < 1000000:
        xlabels = [format((label - subval), ",.0f") for label in ax.get_xticks()]
        lab = "bp"
    elif maxval < 10000000:
        xlabels = [format((label - subval) / 1000, ",.1f") for label in ax.get_xticks()]
        lab = "kbp"
    else:
        xlabels = [format((label - subval) / 1000, ",.1f") for label in ax.get_xticks()]
        lab = "kbp"

    ax.set_ylim(0, ylim)

    ax.set_xlabel("Assembly position ({})".format(lab), fontweight="bold")
    ax.set_ylabel("Sequence read depth", fontweight="bold")

    # Including this causes some internal bug in matplotlib when the font-size changes
    # ylabels = [format(label, ",.0f") for label in ax.get_yticks()]
    # ax.set_yticklabels(ylabels)
    ax.set_xticklabels(xlabels)

    # Hide the right and top spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")

    return fig, ax


def peak_finder(
    data: np.ndarray,
    positions: np.ndarray,
    *,
    height: int,
    added_region_bounds: int,
    distance: int = 100_000,
    width: int = 20,
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
    added_region_bounds: int = 0,
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
    mean_first, stdev_first = df["first"].mean(), df["first"].std()
    mean_second, stdev_second = df["second"].mean(), df["second"].std()

    first_peak_coords = peak_finder(
        df["first"],
        positions,
        height=mean_first * 2,
        added_region_bounds=added_region_bounds,
    )
    first_valley_coords = peak_finder(
        -df["first"],
        positions,
        height=-(mean_first - (3 * stdev_first)),
        added_region_bounds=added_region_bounds,
    )

    # Remove secondary rows that don't meet minimal secondary coverage.
    second_thr = max(round(mean_first * 0.1), round(mean_second + (3 * stdev_second)))

    classified_second_outliers = set()
    df_second_outliers = df.filter(pl.col("second") > second_thr)
    # Group consecutive positions allowing a maximum gap of stepsize.
    # Larger stepsize groups more positions.
    second_outliers_coords = [
        pt.open(grp[0], grp[-1])
        for grp in consecutive(df_second_outliers["position"], stepsize=25_000)
        if len(grp) > 5
    ]

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
        if df_second_outlier_het_ratio["het_ratio"][0] > 0.1:
            misassemblies[Misassembly.COLLAPSE_VAR].add(second_outlier)
        elif df_second_outlier_het_ratio["het_ratio"][0] < 0.4:
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
    rm_output: pl.DataFrame,
) -> pl.DataFrame:
    bam = pysam.AlignmentFile(infile, threads=threads)
    contig_name = f"{contig}:{start}-{end}"

    sys.stderr.write(f"Reading in NucFreq from region: {contig_name}\n")

    df_group_labeled, miassemblies = classify_misassemblies(
        np.flip(
            np.sort(get_coverage_by_base(bam, contig, start, end), axis=1)
        ).transpose(),
        np.arange(start, end),
    )

    if output_dir:
        _ = plot_coverage(df_group_labeled, contig, rm_output)

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
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)

    # Read regions and close file handle to bam.
    bam = pysam.AlignmentFile(args.input_bam, threads=args.threads)
    regions = list(read_regions(bam, args))
    bam.close()

    # set text size
    matplotlib.rcParams.update({"font.size": PLOT_FONT_SIZE})

    # Read repeatmasker. UNTESTED.
    rm_output = read_repeatmasker(args.input_repeatmasker)

    with mp.Pool(processes=args.processes) as pool:
        results = pool.starmap(
            classify_plot_assembly,
            [
                (args.input_bam, args.output_dir, args.threads, *region, rm_output)
                for region in regions
            ],
        )

    # Save misassemblies to output bed
    all_misasm = pl.concat(results).sort(by=["contig", "start"])
    all_misasm.write_csv(file=args.output_bed, include_header=False, separator="\t")


if __name__ == "__main__":
    raise SystemExit(main())
