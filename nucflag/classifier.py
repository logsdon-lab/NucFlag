import gzip
import os
import shutil
import sys
from collections import defaultdict
from typing import Any, DefaultDict

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import portion as pt
import pysam
import scipy.signal

from .constants import PLOT_DPI
from .io import get_coverage_by_base
from .misassembly import Misassembly
from .plot import plot_coverage
from .region import Region, update_relative_ignored_regions


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
    df_cov: pl.DataFrame,
    *,
    config: dict[str, Any],
    ignored_regions: list[Region],
) -> tuple[pl.DataFrame, dict[Misassembly, set[pt.Interval]]]:
    # Filter ignored regions.
    exprs_ignored_region_positions: list[pl.Expr] = [
        ~(
            (pl.col("position") >= r.region.lower)
            & (pl.col("position") <= r.region.upper)
        )
        for r in ignored_regions
        if r.region
    ]

    if df_cov.filter(pl.col("first") != 0).is_empty():
        return (
            df_cov.with_columns(status=pl.lit(Misassembly.GAP)),
            {Misassembly.GAP: {pt.open(df_cov["position"][0], df_cov["position"][-1])}},
        )

    # Calculate std and mean for both most and second most freq read.
    # Remove ignored regions and gaps in coverage which would affect mean.
    df_summary = (
        df_cov.lazy()
        .filter(pl.col("first") != 0)
        .filter(
            exprs_ignored_region_positions if exprs_ignored_region_positions else True
        )
        .select("first", "second")
        .describe()
    )

    mean_first, mean_second = (
        df_summary.filter(pl.col("statistic") == "mean")
        .select("first", "second")
        .row(0)
    )
    stdev_first, stdev_second = (
        df_summary.filter(pl.col("statistic") == "std").select("first", "second").row(0)
    )

    # Calculate misjoin height threshold. Filters for some percent of the mean or static value.
    misjoin_height_thr = (
        mean_first * config["first"]["thr_misjoin_valley"]
        if isinstance(config["first"]["thr_misjoin_valley"], float)
        else config["first"]["thr_misjoin_valley"]
    )
    collapse_height_thr = (
        mean_first * config["first"]["thr_collapse_peak"]
        if isinstance(config["first"]["thr_collapse_peak"], float)
        else config["first"]["thr_collapse_peak"]
    )

    first_peak_height_thr = mean_first + (
        config["first"]["thr_peak_height_std_above"] * stdev_first
    )
    first_valley_height_thr = mean_first - (
        config["first"]["thr_valley_height_std_below"] * stdev_first
    )

    first_peak_coords = peak_finder(
        df_cov["first"],
        df_cov["position"],
        height=first_peak_height_thr,
        distance=config["first"]["thr_min_peak_horizontal_distance"],
        width=config["first"]["thr_min_peak_width"],
        group_distance=config["first"]["peak_group_distance"],
    )
    first_valley_coords = peak_finder(
        -df_cov["first"],
        df_cov["position"],
        # Account for when thr goes negative.
        height=-(
            misjoin_height_thr
            if first_valley_height_thr < 0
            else first_valley_height_thr
        ),
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
    df_second_outliers = df_cov.filter(pl.col("second") > second_thr)
    # Group consecutive positions allowing a maximum gap of stepsize.
    # Larger stepsize groups more positions.
    second_outliers_coords: list[pt.Interval] = []
    for grp in consecutive(
        df_second_outliers["position"], stepsize=config["second"]["group_distance"]
    ):
        if len(grp) < config["second"]["thr_min_group_size"]:
            continue
        second_outliers_coords.append(pt.open(grp[0], grp[-1]))

    misassemblies: dict[Misassembly, set[pt.Interval]] = {m: set() for m in Misassembly}

    # Intersect intervals and classify collapses.
    for peak in first_peak_coords:
        local_max_collapse_first = (
            df_cov.filter(filter_interval_expr(peak)).max().get_column("first")[0]
        )
        includes_het = False

        for second_outlier in second_outliers_coords:
            # If local max of suspected collapsed region is greater than thr, is a collapse.
            if local_max_collapse_first < first_peak_height_thr:
                continue

            if peak.overlaps(second_outlier):
                misassemblies[Misassembly.COLLAPSE_VAR].add(peak.union(second_outlier))
                classified_second_outliers.add(second_outlier)
                includes_het = True

        if local_max_collapse_first >= collapse_height_thr and not includes_het:
            misassemblies[Misassembly.COLLAPSE].add(peak)

    # Classify gaps.
    df_gaps = df_cov.filter(pl.col("first") == 0)
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
        includes_het = False
        for second_outlier in second_outliers_coords:
            if valley.overlaps(second_outlier):
                # Merge intervals.
                misassemblies[Misassembly.MISJOIN].add(valley.union(second_outlier))
                includes_het = True
                classified_second_outliers.add(second_outlier)

        # Otherwise, check if valley's median falls below threshold.
        # This means that while no overlapping secondary reads, low coverage means likely elsewhere.
        # Treat as a misjoin.
        if not includes_het:
            # Filter first to get general region.
            df_valley = df_cov.filter(filter_interval_expr(valley)).filter(
                pl.col("first") <= misjoin_height_thr
            )
            # Skip if fewer than 2 points found.
            if df_valley.shape[0] < config["first"]["thr_min_valley_width"]:
                continue

            # Get bounds of region and calculate median.
            # Avoid flagging if intersects gap region.
            df_valley = (
                df_valley
                if df_valley.shape[0] == 1
                else df_cov.filter(
                    filter_interval_expr(
                        pt.open(
                            df_valley["position"].min(), df_valley["position"].max()
                        )
                    )
                )
            )
            if df_valley["first"].min() <= misjoin_height_thr and not any(
                g.overlaps(valley) for g in misassemblies[Misassembly.GAP]
            ):
                misassemblies[Misassembly.MISJOIN].add(valley)

    # Assign heterozygous site for remaining secondary regions not categorized.
    for second_outlier in second_outliers_coords:
        if second_outlier in classified_second_outliers:
            continue

        misassemblies[Misassembly.HET].add(second_outlier)

    # Annotate df with misassembly.
    lf = df_cov.lazy().with_columns(status=pl.lit("Good"))

    filtered_misassemblies = defaultdict(set)
    for mtype, regions in misassemblies.items():
        remove_regions = set()
        for region in regions:
            # Remove ignored regions.
            # TODO: This still could be better. O(n). Would need to rewrite and use an interval tree.
            if any(
                ignored_region.region.overlaps(region)
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
    infile: str,
    output_dir: str | None,
    output_cov_dir: str | None,
    threads: int,
    contig: str,
    start: int,
    end: int,
    config: dict[str, Any],
    overlay_regions: DefaultDict[int, set[Region]] | None,
    ignored_regions: set[Region],
) -> pl.DataFrame:
    contig_name = f"{contig}:{start}-{end}"
    sys.stderr.write(f"Reading in NucFreq from region: {contig_name}\n")

    try:
        bam = pysam.AlignmentFile(infile, threads=threads)
        cov_first_second = np.flip(
            np.sort(get_coverage_by_base(bam, contig, start, end), axis=1).transpose(),
            axis=0,
        )
        df = pl.DataFrame(
            {
                "position": np.arange(start, end),
                "first": cov_first_second[0],
                "second": cov_first_second[1],
            },
            schema={"position": pl.Int64, "first": pl.Int16, "second": pl.Int16},
        )
        del cov_first_second
    except ValueError:
        df = pl.read_csv(
            infile,
            separator="\t",
            has_header=True,
            dtypes={"position": pl.Int64, "first": pl.Int16, "second": pl.Int16},
        )

    # Update ignored regions if relative.
    updated_ignored_regions = list(
        update_relative_ignored_regions(ignored_regions, ctg_start=start, ctg_end=end)
    )

    df_group_labeled, misassemblies = classify_misassemblies(
        df,
        config=config,
        ignored_regions=updated_ignored_regions,
    )

    if output_dir:
        _ = plot_coverage(df_group_labeled, misassemblies, contig, overlay_regions)

        sys.stderr.write(f"Plotting {contig_name}.\n")

        output_plot = os.path.join(output_dir, f"{contig_name}.png")
        plt.savefig(output_plot, dpi=PLOT_DPI, bbox_inches="tight")

    if output_cov_dir:
        sys.stderr.write(f"Writing coverage bed file for {contig_name}.\n")

        output_bed = os.path.join(output_cov_dir, f"{contig_name}.bed")
        df_group_labeled.write_csv(output_bed, separator="\t")

        sys.stderr.write(f"Compressing coverage bed file for {contig_name}.\n")
        with open(output_bed, "rb") as f_in:
            with gzip.open(f"{output_bed}.gz", "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(output_bed)

    df_misassemblies = pl.DataFrame(
        [
            (contig_name, interval.lower, interval.upper, misasm)
            for misasm, intervals in misassemblies.items()
            for interval in intervals
        ],
        schema=["contig", "start", "stop", "misassembly"],
    )
    return df_misassemblies
