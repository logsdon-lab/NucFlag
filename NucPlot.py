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

from typing import Generator, Any

warnings.filterwarnings("ignore")

PLOT_COLORS = sns.color_palette()
PLOT_FONT_SIZE = 16
PLOT_REGION_HEIGHT = 4
PLOT_WIDTH = 16
PLOT_DPI = 600


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use per-base read coverage to classify/plot misassemblies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--infile", help="Input bam file.")
    parser.add_argument(
        "-o",
        "--output",
        help="Output plot file or plot dir. If dir, each contig is generated as an individual plot.",
        required=False,
        default="nucplot.png",
    )
    parser.add_argument(
        "-r",
        "--repeatmasker",
        help="RepeatMasker output to add to plot.",
        type=argparse.FileType("r"),
        default=None,
    )
    parser.add_argument(
        "--regions", nargs="*", help="Regions with the format: (.*):(\\d+)-(\\d+)"
    )
    parser.add_argument("--bed", default=None, help="Bed file with regions to plot.")
    parser.add_argument(
        "--obed", default=None, help="Bed file with misassembled regions."
    )
    parser.add_argument("--threads", default=4, help="Threads for reading bam file.")

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

    if args.regions is not None or args.bed is not None:
        sys.stderr.write("Reading in the region or bed argument(s).\n")
        if args.regions is not None:
            for region in args.regions:
                match = re.match(r"(.+):(\d+)-(\d+)", region)
                assert match, region + " not valid!"
                chrm, start, end = match.groups()
                refs[chrm] = [int(start), int(end)]
                yield (chrm, int(start), int(end))

        if args.bed is not None:
            for line in open(args.bed):
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


def read_repeatmasker(args: argparse.Namespace) -> pl.DataFrame | None:
    rm_output = None
    if args.repeatmasker is not None:
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
                args.repeatmasker,
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

        args.repeatmasker.close()

    return rm_output


def plot_coverage(
    axs: Any | None,
    output: str | None,
    group_n: int,
    contig_name: str,
    first: np.ndarray,
    second: np.ndarray,
    truepos: np.ndarray,
    rm_output: pl.DataFrame | None,
):
    ylim = int(first.max() * 1.05)

    # get the correct axis
    if axs:
        ax = axs[group_n]
    else:
        fig, ax = plt.subplots(figsize=(PLOT_WIDTH, PLOT_REGION_HEIGHT))

    if rm_output:
        rmax = ax
        sys.stderr.write("Subsetting the repeatmasker file.\n")
        rm = rm_output.filter(
            (pl.col("qname") == contig_name)
            & (pl.col("start") >= min(truepos))
            & (pl.col("end") <= max(truepos))
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
        truepos,
        first,
        "o",
        color="black",
        markeredgewidth=0.0,
        markersize=2,
        label="most frequent base pair",
    )
    (sec,) = ax.plot(
        truepos,
        second,
        "o",
        color="red",
        markeredgewidth=0.0,
        markersize=2,
        label="second most frequent base pair",
    )

    maxval = max(truepos)
    minval = min(truepos)
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

    sys.stderr.write(f"Added contig, {title}, to plot.\n")


def peak_finder(
    data: np.ndarray,
    positions: np.ndarray,
    *,
    height: int,
    added_region_bounds: int,
    distance: int = 100_000,
    width: int = 500,
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
    added_region_bounds: int = 100_000,
) -> pl.DataFrame:
    first: np.ndarray = cov_first_second[0]
    second: np.ndarray = cov_first_second[1]
    del cov_first_second

    # Calculate std and mean for both most and second most freq read.
    mean_first, stdev_first = first.mean(), first.std()
    mean_second, stdev_second = second.mean(), second.std()

    first_peak_coords = peak_finder(
        first,
        positions,
        height=mean_first + (3 * stdev_first),
        added_region_bounds=added_region_bounds,
    )
    first_valley_coords = peak_finder(
        -first,
        positions,
        height=-(mean_first - (3 * stdev_first)),
        added_region_bounds=added_region_bounds,
    )

    # Remove secondary rows that don't meet minimal secondary coverage.
    second_thr = max(round(mean_first * 0.1), round(mean_second + (3 * stdev_second)))

    second_outliers, *_ = np.where(second > second_thr)
    # Group consecutive positions allowing a maximum gap of stepsize.
    second_outliers_coords = [
        pt.open(positions[grp[0]], positions[grp[-1]])
        for grp in consecutive(second_outliers, stepsize=10_000)
        if len(grp) > 5
    ]

    # Intersect intervals and classify.
    collapse_w_no_variant = []
    collapse_w_variant: set[pt.Interval] = set()
    for peak in first_peak_coords:
        for second_outlier in second_outliers_coords:
            if second_outlier in peak:
                collapse_w_variant.add(peak)

        if peak not in collapse_w_variant:
            collapse_w_no_variant.append(peak)

    # TODO: Gaps here?
    raise NotImplementedError

    misjoins: list[pt.Interval] = [
        valley
        for valley in first_valley_coords
        for second_outlier in second_outliers_coords
        if second_outlier in valley
    ]
    print(misjoins)
    # TODO: false dupes

    return pl.DataFrame()


def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.infile, threads=args.threads)
    regions = list(read_regions(bam, args))

    # SET up the plot based on the number of regions
    num_regions = len(regions)
    height = num_regions * PLOT_REGION_HEIGHT

    # set text size
    matplotlib.rcParams.update({"font.size": PLOT_FONT_SIZE})

    # TODO: Change so plot fn returns figure and generate subplot at end.
    # make axes
    fig, axs = None, None
    output_dir = args.output if os.path.isdir(args.output) else None
    if not output_dir:
        fig, axs = plt.subplots(
            nrows=num_regions, ncols=1, figsize=(PLOT_WIDTH, height)
        )
        if num_regions == 1:
            axs = [axs]
        figsize = fig.get_size_inches() * PLOT_DPI

        if (figsize > (2**16)).any():
            raise ValueError("Too many regions provided. Reduce the number of regions.")

    # rm_output = read_repeatmasker(args)

    for group_n, (contig, start, end) in enumerate(regions):
        sys.stderr.write(
            "Reading in NucFreq from region: {}:{}-{}\n".format(contig, start, end)
        )

        df_group_labeled = classify_misassemblies(
            np.flip(
                np.sort(get_coverage_by_base(bam, contig, start, end), axis=1)
            ).transpose(),
            np.arange(start, end),
        )
        print(df_group_labeled)
        # plot_coverage(axs, output_dir, group_n, contig, first, second, position, rm_output)

    # plt.tight_layout()
    # plt.savefig(args.outfile, dpi=PLOT_DPI)


if __name__ == "__main__":
    raise SystemExit(main())
