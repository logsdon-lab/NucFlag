import gzip
import os
import shutil
import sys
from collections import defaultdict
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pysam
import scipy.signal
from intervaltree import Interval, IntervalTree

from .constants import PLOT_DPI, PROP_INCLUDE
from .io import get_coverage_by_base
from .misassembly import Misassembly
from .plot import plot_coverage
from .region import Region, update_relative_ignored_regions


def peak_finder(
    df: pl.DataFrame,
    *,
    height: int,
    width: int,
    invert: bool = False,
    group_distance: int = 5_000,
) -> IntervalTree:
    # Including zeroes will alter the calculation of valley heights as will always be the max.
    df_subset = df.filter(pl.col("first") != 0)
    data = df_subset["first"] if not invert else -df_subset["first"]
    positions = df_subset["position"]
    # Use height as first threshold.
    peaks, peak_info = scipy.signal.find_peaks(data, height=height, width=width)
    # Calculate median avoiding peaks.
    median_no_peaks = data.filter(~positions.is_in(peaks)).median()

    # Create interval tree adjusted to group distance.
    intervals = IntervalTree()
    for left_idx, right_idx, peak_ht in zip(
        peak_info["left_ips"], peak_info["right_ips"], peak_info["peak_heights"]
    ):
        left_idx, right_idx = int(left_idx), int(right_idx)
        left_pos, right_pos = positions[left_idx], positions[right_idx]
        # Calculate relative height
        peak_rel_ht = peak_ht - median_no_peaks
        intervals.add(
            Interval(
                left_pos - group_distance,
                right_pos + group_distance,
                peak_rel_ht,
            )
        )
    # Merge taking largest height.
    intervals.merge_overlaps(strict=False, data_reducer=lambda x, y: max(x, y))

    return IntervalTree(
        Interval(
            interval.begin + group_distance,
            interval.end - group_distance,
            interval.data,
        )
        for interval in intervals.iter()
    )


# https://stackoverflow.com/a/7353335
def consecutive(data, stepsize: int = 1):
    return np.split(data, np.where((np.diff(data) <= stepsize) == False)[0] + 1)  # noqa: E712


def filter_interval_expr(interval: Interval, *, col: str = "position") -> pl.Expr:
    return pl.col(col).is_between(interval.begin, interval.end)


def calculate_het_ratio(
    df: pl.DataFrame, interval: Interval, second_thr: float
) -> float:
    df_het = df.filter(filter_interval_expr(interval)).filter(
        pl.col("second") > second_thr
    )
    # Apply median filter with small window size to remove outliers while not oversmoothing.
    first_signal: np.ndarray = scipy.signal.medfilt(df_het["first"], 3)
    second_signal: np.ndarray = scipy.signal.medfilt(df_het["second"], 3)
    # Use max to get general amplitude of data.
    het_first_median = first_signal.max()
    het_second_median = second_signal.max()
    return het_second_median / (het_first_median + het_second_median)


def classify_misassemblies(
    df_cov: pl.DataFrame,
    *,
    config: dict[str, Any],
    ignored_regions: list[Region],
) -> tuple[pl.DataFrame, dict[Misassembly, IntervalTree]]:
    ignored_regions_intervals = IntervalTree(
        r.region for r in ignored_regions if r.region
    )
    # Filter ignored regions.
    exprs_ignored_region_positions: list[pl.Expr] = [
        ~pl.col("position").is_between(r.begin, r.end)
        for r in ignored_regions_intervals.iter()
    ]

    df_cov = df_cov.with_columns(
        include=pl.when(
            exprs_ignored_region_positions if exprs_ignored_region_positions else True
        )
        .then(True)
        .otherwise(False)
    )

    # Calculate std and mean for both most and second most freq read.
    # Remove ignored regions and gaps in coverage which would affect mean.
    df_summary = (
        df_cov.lazy()
        .filter((pl.col("first") != 0) & pl.col("include"))
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

    if any(
        metric is None
        for metric in (mean_first, mean_second, stdev_first, stdev_second)
    ):
        return (
            df_cov.with_columns(status=pl.lit(Misassembly.GAP)),
            {
                Misassembly.GAP: {
                    Interval(df_cov["position"][0], df_cov["position"][-1])
                }
            },
        )

    collapse_height_thr = round(mean_first * config["first"]["thr_collapse_peak"])
    misjoin_height_thr = round(mean_first * config["first"]["thr_misjoin_valley"])
    first_peak_coords = peak_finder(
        df_cov,
        height=collapse_height_thr,
        width=config["first"]["thr_min_peak_width"],
        group_distance=config["first"]["peak_group_distance"],
    )
    first_valley_coords = peak_finder(
        df_cov,
        invert=True,
        height=-mean_first,
        width=config["first"]["thr_min_valley_width"],
        group_distance=config["first"]["valley_group_distance"],
    )

    # Remove secondary rows that don't meet minimal secondary coverage.
    second_thr = round(mean_first * config["second"]["thr_min_perc_first"])

    classified_second_outliers = set()
    df_second_outliers = df_cov.filter(pl.col("second") > second_thr)
    # Group consecutive positions allowing a maximum gap of stepsize.
    # Larger stepsize groups more positions.
    second_outliers_coords: IntervalTree = IntervalTree()
    for grp in consecutive(
        df_second_outliers["position"], stepsize=config["second"]["group_distance"]
    ):
        if len(grp) < config["second"]["thr_min_group_size"]:
            continue
        second_outliers_coords.add(Interval(grp[0], grp[-1]))

    misassemblies: defaultdict[Misassembly, IntervalTree] = defaultdict(IntervalTree)

    # Intersect intervals and classify collapses.
    for peak in first_peak_coords.iter():
        # If height of suspected collapsed region is greater than thr, is a collapse.
        if peak.data < collapse_height_thr:
            continue

        overlaps: set[Interval] = second_outliers_coords.overlap(peak)
        for overlap in overlaps:
            new_overlap_interval = Interval(
                min(peak.begin, overlap.begin), max(peak.end, overlap.end)
            )
            misassemblies[Misassembly.COLLAPSE_VAR].add(new_overlap_interval)
            classified_second_outliers.add(overlap)

        if not overlaps:
            misassemblies[Misassembly.COLLAPSE].add(peak)

    # Classify gaps.
    df_gaps = df_cov.filter(pl.col("first") == 0)
    for grp in consecutive(df_gaps["position"], stepsize=1):
        if grp.size == 0:
            continue
        # In case of 1 len interval
        if grp[0] == grp[-1]:
            gap_len = 1
            gap = Interval(grp[0], grp[0] + 1)
        else:
            gap_len = grp[-1] - grp[0]
            gap = Interval(grp[0], grp[-1])

        if gap_len < config["gaps"]["thr_max_allowed_gap_size"]:
            continue

        misassemblies[Misassembly.GAP].add(gap)

    # Classify misjoins.
    for valley in first_valley_coords:
        second_overlaps = second_outliers_coords.overlap(valley)
        overlaps_gap = misassemblies[Misassembly.GAP].overlaps(valley)
        overlaps_collapse = misassemblies[Misassembly.COLLAPSE].overlaps(
            valley
        ) or misassemblies[Misassembly.COLLAPSE_VAR].overlaps(valley)

        # Ignore if overlaps existing gap, collapse, or collapse with variant.
        if overlaps_gap or overlaps_collapse:
            continue

        # Calculate relative height of valley.
        # And only take those that meet criteria.
        valley_below_thr = valley.data > misjoin_height_thr

        for overlap in second_overlaps:
            # Merge intervals.
            new_overlap_interval = Interval(
                min(valley.begin, overlap.begin), max(valley.end, overlap.end)
            )
            if valley_below_thr:
                misassemblies[Misassembly.MISJOIN].add(new_overlap_interval)
                classified_second_outliers.add(overlap)
            else:
                het_ratio = calculate_het_ratio(df_cov, overlap, second_thr)
                if het_ratio >= config["second"]["thr_het_ratio"]:
                    misassemblies[Misassembly.ERROR].add(new_overlap_interval)
                    classified_second_outliers.add(overlap)

        if not second_overlaps and valley_below_thr:
            misassemblies[Misassembly.MISJOIN].add(valley)

    # Assign heterozygous site for remaining secondary regions not categorized.
    for second_outlier in second_outliers_coords:
        if second_outlier in classified_second_outliers:
            continue

        het_ratio = calculate_het_ratio(df_cov, second_outlier, second_thr)

        # If under het ratio, treat as het. Otherwise, some error.
        if het_ratio >= config["second"]["thr_het_ratio"]:
            misassemblies[Misassembly.ERROR].add(second_outlier)
        else:
            misassemblies[Misassembly.HET].add(second_outlier)

    # Annotate df with misassembly.
    lf = df_cov.lazy().with_columns(status=pl.lit("Good"))
    for mtype, regions in misassemblies.items():
        # Remove any overlapping regions.
        regions.merge_overlaps(strict=False)

        removed_regions = set()
        for region in regions.iter():
            region_len = region.begin - region.end
            expr_filter_region = filter_interval_expr(region)

            # Remove ignored regions.
            if ignored_regions_intervals.overlaps(region):
                df_include = (
                    df_cov.filter(expr_filter_region)
                    .get_column("include")
                    .value_counts()
                )
                try:
                    include_rows = df_include.filter(pl.col("include")).get_column(
                        "count"
                    )[0]
                except IndexError:
                    include_rows = 0

                # Include misassemblies that overlap with greater than 50% non-ignored region.
                prop_include = include_rows / region_len
                if prop_include > PROP_INCLUDE:
                    # And update its status.
                    lf = lf.with_columns(
                        status=pl.when(expr_filter_region)
                        .then(pl.lit(mtype))
                        .otherwise(pl.col("status"))
                    )
                else:
                    # Otherwise, remove it.
                    removed_regions.add(region)

        for rm_region in removed_regions:
            regions.remove_overlap(rm_region.begin, rm_region.end)

    # TODO: false dupes
    return lf.collect(), misassemblies


def classify_plot_assembly(
    infile: str,
    output_dir: str | None,
    output_cov_dir: str | None,
    threads: int,
    contig: str,
    start: int,
    end: int,
    config: dict[str, Any],
    overlay_regions: defaultdict[int, set[Region]] | None,
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
    del df

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

    del df_group_labeled

    return pl.DataFrame(
        [
            (contig_name, interval.begin, interval.end, misasm)
            for misasm, intervals in misassemblies.items()
            for interval in intervals
        ],
        schema=["contig", "start", "stop", "misassembly"],
        orient="row",
    )
