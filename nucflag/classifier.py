import os
import sys
import pysam
import scipy.signal

import numpy as np
import polars as pl
import portion as pt
import matplotlib.pyplot as plt

from .io import get_coverage_by_base
from .plot import plot_coverage
from .constants import PLOT_DPI
from .misassembly import Misassembly
from .region import Region

from typing import Any
from collections import defaultdict


def peak_finder(
    data: np.ndarray,
    positions: np.ndarray,
    *,
    height: int,
    distance: int,
    width: int,
    group_distance: int = 5_000,
) -> list[pt.Interval]:
    _, peak_info = scipy.signal.find_peaks(
        data, height=height, distance=distance, width=width
    )

    # Peaks are passed in sorted order.
    intervals: list[pt.Interval] = []
    for left_pos, right_pos in zip(peak_info["left_ips"], peak_info["right_ips"]):
        new_peak = pt.closed(positions[int(left_pos)], positions[int(right_pos)])
        try:
            prev_peak = intervals.pop()
            if prev_peak == new_peak:
                intervals.append(new_peak)
                continue

            intersection = prev_peak.intersection(new_peak)
            dst_between = new_peak.lower - prev_peak.upper

            # Merge peaks if there are any intersections.
            if not intersection.empty:
                new_peak = prev_peak.union(new_peak)
                intervals.append(new_peak)
            # Merge peaks if within a set distance.
            elif dst_between < group_distance:
                new_peak = pt.closed(prev_peak.lower, new_peak.upper)
                intervals.append(new_peak)
            else:
                intervals.append(prev_peak)
                intervals.append(new_peak)
        # First peak.
        except IndexError:
            intervals.append(new_peak)

    return intervals


# https://stackoverflow.com/a/7353335
def consecutive(data, stepsize: int = 1):
    return np.split(data, np.where((np.diff(data) <= stepsize) == False)[0] + 1)  # noqa: E712


def filter_interval_expr(interval: pt.Interval, *, col: str = "position") -> pl.Expr:
    return (pl.col(col) >= interval.lower) & (pl.col(col) <= interval.upper)


def classify_misassemblies(
    cov_first_second: np.ndarray,
    positions: np.ndarray,
    *,
    config: dict[str, Any],
    ignored_regions: list[Region] | None,
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

    first_peak_height_thr = mean_first + (
        config["first"]["thr_peak_height_std_above"] * stdev_first
    )
    first_valley_height_thr = mean_first - (
        config["first"]["thr_valley_height_std_below"] * stdev_first
    )

    first_peak_coords = peak_finder(
        df["first"],
        positions,
        height=first_peak_height_thr,
        distance=config["first"]["thr_min_peak_horizontal_distance"],
        width=config["first"]["thr_min_peak_width"],
        group_distance=config["first"]["peak_group_distance"],
    )
    first_valley_coords = peak_finder(
        -df["first"],
        positions,
        height=-first_valley_height_thr,
        distance=config["first"]["thr_min_valley_horizontal_distance"],
        width=config["first"]["thr_min_valley_width"],
        group_distance=config["first"]["valley_group_distance"],
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
        second_outliers_coords.append(pt.open(grp[0], grp[-1]))

    misassemblies: dict[Misassembly, set[pt.Interval]] = {m: set() for m in Misassembly}

    # Intersect intervals and classify collapses.
    for peak in first_peak_coords:
        for second_outlier in second_outliers_coords:
            if second_outlier in peak:
                misassemblies[Misassembly.COLLAPSE_VAR].add(peak)
                classified_second_outliers.add(second_outlier)

        if peak not in misassemblies[Misassembly.COLLAPSE_VAR]:
            local_mean_collapse_first = (
                df.filter(filter_interval_expr(peak)).mean().get_column("first")[0]
            )
            # If local mean of suspected collapsed region is greater than thr, is a collapse.
            if local_mean_collapse_first > first_peak_height_thr:
                misassemblies[Misassembly.COLLAPSE].add(peak)

    # Classify gaps.
    df_gaps = df.filter(pl.col("first") == 0)
    gaps = set()
    for grp in consecutive(df_gaps["position"], stepsize=1):
        if len(grp) < 2:
            continue

        gap_len = grp[-1] - grp[0]

        if gap_len < config["gaps"]["thr_max_allowed_gap_size"]:
            continue

        gap = pt.open(grp[0], grp[-1])
        gaps.add(gap)

    misassemblies[Misassembly.GAP] = gaps

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
            filter_interval_expr(second_outlier) & (pl.col("second") != 0)
        )
        df_second_outlier_het_ratio = df_second_outlier.mean().with_columns(
            het_ratio=pl.col("second") / (pl.col("first") + pl.col("second"))
        )
        # Use het ratio to classify.
        # Low ratio consider collapse with var.
        if (
            df_second_outlier_het_ratio["het_ratio"][0]
            < config["second"]["thr_collapse_het_ratio"]
        ):
            misassemblies[Misassembly.COLLAPSE_VAR].add(second_outlier)
        else:
            misassemblies[Misassembly.MISJOIN].add(second_outlier)

    # Annotate df with misassembly.
    lf = df.lazy().with_columns(status=pl.lit("Good"))

    contig_region = pt.open(positions[0], positions[-1])
    if not ignored_regions:
        ignored_regions = []

    filtered_misassemblies = defaultdict(set)
    for mtype, regions in misassemblies.items():
        remove_regions = set()
        for region in regions:
            # Remove ignored regions.
            if any(
                ignored_region.contains(region, full=contig_region)
                for ignored_region in ignored_regions
            ):
                remove_regions.add(region)
                continue

            lf = lf.with_columns(
                status=pl.when(filter_interval_expr(region))
                .then(pl.lit(mtype))
                .otherwise(pl.col("status"))
            )
        filtered_misassemblies[mtype] = regions - remove_regions

    # TODO: false dupes
    return lf.collect(), filtered_misassemblies


def classify_plot_assembly(
    input_bam: str,
    output_dir: str | None,
    threads: int,
    contig: str,
    start: int,
    end: int,
    config: dict[str, Any],
    ignored_regions: list[Region] | None,
) -> pl.DataFrame:
    bam = pysam.AlignmentFile(input_bam, threads=threads)
    contig_name = f"{contig}:{start}-{end}"

    sys.stderr.write(f"Reading in NucFreq from region: {contig_name}\n")

    df_group_labeled, misassemblies = classify_misassemblies(
        np.flip(
            np.sort(get_coverage_by_base(bam, contig, start, end), axis=1)
        ).transpose(),
        np.arange(start, end),
        config=config,
        ignored_regions=ignored_regions,
    )

    if output_dir:
        _ = plot_coverage(df_group_labeled, misassemblies, contig)

        sys.stderr.write(f"Plotted {contig_name}.\n")

        output_plot = os.path.join(output_dir, f"{contig_name}.png")
        plt.tight_layout()
        plt.savefig(output_plot, dpi=PLOT_DPI)

    df_misassemblies = pl.DataFrame(
        [
            (contig_name, interval.lower, interval.upper, misasm)
            for misasm, intervals in misassemblies.items()
            for interval in intervals
        ],
        schema=["contig", "start", "stop", "misassembly"],
    )
    return df_misassemblies
