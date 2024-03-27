#!/usr/bin/env python
import re
import sys
import pysam
import warnings
import argparse
import matplotlib

matplotlib.use("agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

from typing import Generator

warnings.filterwarnings("ignore")

PLOT_COLORS = sns.color_palette()
PLOT_FONT_SIZE = 16
PLOT_REGION_HEIGHT = 4
PLOT_WIDTH = 16
PLOT_DPI = 600


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use nucleotide coverage to classify misassemblies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--infile", help="Input bam file.")
    parser.add_argument(
        "-o",
        "--outfile",
        help="Output plot file.",
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


def read_repeatmasker(args: argparse.Namespace) -> pd.DataFrame | None:
    RM = None
    cmap = {}
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
        lines = []
        for idx, line in enumerate(args.repeatmasker):
            if idx > 2:
                lines.append(line.strip().split()[0:15])

        RM = pd.DataFrame(lines, columns=names)
        RM.start = RM.start.astype(int)
        RM.end = RM.end.astype(int)
        RM["label"] = RM.family.str.replace("/.*", "")
        for idx, lab in enumerate(sorted(RM.label.unique())):
            cmap[lab] = PLOT_COLORS[idx % len(PLOT_COLORS)]
        RM["color"] = RM.label.map(cmap)

        args.repeatmasker.close()

    return RM


def plot_coverage(
    axs,
    group_n: int,
    contig_name: str,
    first: np.ndarray,
    second: np.ndarray,
    truepos: np.ndarray,
    rm_output: pd.DataFrame | None,
):
    ylim = int(first.max() * 1.05)

    # if args.obed:
    #     tmp = group.loc[
    #         group.second >= args.minobed,
    #         ["contig", "position", "position", "first", "second"],
    #     ]
    #     if counter == 0:
    #         tmp.to_csv(
    #             args.obed,
    #             header=["#contig", "start", "end", "first", "second"],
    #             sep="\t",
    #             index=False,
    #         )
    #     else:
    #         tmp.to_csv(args.obed, mode="a", header=None, sep="\t", index=False)

    # get the correct axis
    ax = axs[group_n]

    if rm_output:
        rmax = ax
        sys.stderr.write("Subsetting the repeatmakser file.\n")
        rm = rm_output[
            (rm_output.qname == contig_name)
            & (rm_output.start >= min(truepos))
            & (rm_output.end <= max(truepos))
        ]
        assert len(rm.index) != 0, "No matching RM contig"

        rmlength = len(rm.index) * 1.0
        height_offset = ylim / 20
        for rmcount, row in enumerate(rm.iterrows()):
            sys.stderr.write(
                "\rDrawing the {} repeatmasker rectangles:\t{:.2%}".format(
                    rmlength, rmcount / rmlength
                )
            )
            width = row.end - row.start
            rect = patches.Rectangle(
                (row.start, ylim - height_offset),
                width,
                height_offset,
                linewidth=1,
                edgecolor="none",
                facecolor=row.color,
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


def main():
    args = parse_args()
    bam = pysam.AlignmentFile(args.infile, threads=args.threads)
    regions = list(read_regions(bam, args))

    # SET up the plot based on the number of regions
    num_regions = len(regions)
    height = num_regions * PLOT_REGION_HEIGHT

    # set text size
    matplotlib.rcParams.update({"font.size": PLOT_FONT_SIZE})

    # make axes
    fig, axs = plt.subplots(nrows=num_regions, ncols=1, figsize=(PLOT_WIDTH, height))
    if num_regions == 1:
        axs = [axs]
    figsize = fig.get_size_inches() * PLOT_DPI

    if (figsize > (2**16)).any():
        raise ValueError("Too many regions provided. Reduce the number of regions.")

    # # make space for the bottom label of the plot
    # # fig.subplots_adjust(bottom=0.2)
    # # set figure YLIM
    # YLIM = int(max(df["first"]) * 1.05)

    rm_output = read_repeatmasker(args)

    for group_n, (contig, start, end) in enumerate(regions):
        sys.stderr.write(
            "Reading in NucFreq from region: {}:{}-{}\n".format(contig, start, end)
        )

        sorted_cov = np.flip(
            np.sort(get_coverage_by_base(bam, contig, start, end), axis=1)
        ).transpose()

        first: np.ndarray = sorted_cov[0]
        second: np.ndarray = sorted_cov[1]
        position = np.arange(start, end)

        plot_coverage(axs, group_n, contig, first, second, position, rm_output)

    plt.tight_layout()
    plt.savefig(args.outfile, dpi=PLOT_DPI)


if __name__ == "__main__":
    raise SystemExit(main())
