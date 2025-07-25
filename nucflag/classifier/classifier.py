import os
import sys
from collections import defaultdict
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pysam
from intervaltree import Interval, IntervalTree

from .het import identify_hets
from .collapse import get_secondary_allele_coords, identify_collapses
from .misjoin import identify_misjoins, identify_zero_cov_regions
from .common import peak_finder, filter_interval_expr
from ..utils import check_indexed
from ..io import get_coverage_by_base, write_bigwig
from ..misassembly import Misassembly
from ..plot import plot_coverage
from ..region import Region, update_relative_ignored_regions


PLOT_DPI = 600


def filter_misassembled_ignored_regions(
    df_cov: pl.DataFrame,
    misassemblies: defaultdict[Misassembly, IntervalTree],
    ignored_regions: IntervalTree,
    ignored_mtypes: set[Misassembly],
) -> pl.DataFrame:
    """
    Annotates input coverage dataframe with status column. Either good or misassembly type.
    Also filters misassemblies within ignored regions.
    """
    # Annotate df with misassembly.
    lf = df_cov.lazy().with_columns(status=pl.lit("Good"))
    # TODO: Might overflow stack.
    for mtype, regions in misassemblies.items():
        # Remove any overlapping regions.
        regions.merge_overlaps(strict=False)

        # Copy the misassembled intervals as we need to modify the intervaltree.
        for region in regions.items():
            # Trim misassembly removing ignored region.
            for ignored_region in ignored_regions.overlap(region):
                regions.chop(ignored_region.begin, ignored_region.end)

        for region in regions.iter():
            # And update its status.
            lf = lf.with_columns(
                status=pl.when(filter_interval_expr(region))
                .then(pl.lit(mtype))
                .otherwise(pl.col("status"))
            )

    # Remove misassembly types.
    for mtype in ignored_mtypes:
        lf = lf.with_columns(
            status=pl.when(pl.col("status") == str(mtype))
            .then(pl.lit("good"))
            .otherwise(pl.col("status"))
        )
        try:
            misassemblies.pop(mtype)
        except KeyError:
            pass

    return lf.collect()


def classify_misassemblies(
    df_cov: pl.DataFrame,
    *,
    config: dict[str, Any],
    ignored_regions: list[Region],
    ignored_mtypes: set[Misassembly],
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

    # If can't calculate metrics, means entire region has 0 coverage or is ignored.
    # Just call misjoin.
    if any(
        metric is None
        for metric in (mean_first, mean_second, stdev_first, stdev_second)
    ):
        return (
            df_cov.with_columns(status=pl.lit(Misassembly.MISJOIN)),
            {
                Misassembly.MISJOIN: {
                    Interval(df_cov["position"][0], df_cov["position"][-1])
                }
            },
        )

    collapse_height_thr = round(mean_first * config["first"]["thr_collapse_peak"])
    misjoin_height_thr = round(mean_first * config["first"]["thr_misjoin_valley"])

    # Including zeroes will alter the calculation of valley heights as will always be the max.
    # Set ignored regions to median. This avoids calling misassemblies outside that can extend to valid region.
    df_subset = df_cov.with_columns(
        first=pl.when(~pl.col("include"))
        .then(pl.col("first").median())
        .otherwise(pl.col("first")),
        second=pl.when(~pl.col("include"))
        .then(pl.col("second").median())
        .otherwise(pl.col("second")),
    )
    first_data = df_subset["first"]
    positions = df_subset["position"]

    mean_no_peaks, first_peak_coords = peak_finder(
        data=first_data,
        positions=positions,
        height_thr=collapse_height_thr - mean_first,
        abs_height_thr=collapse_height_thr,
        width=config["first"]["thr_min_peak_width"],
        group_distance=config["first"]["peak_group_distance"],
    )
    mean_no_valleys, first_valley_coords = peak_finder(
        data=-first_data,
        positions=positions,
        height_thr=misjoin_height_thr,
        abs_height_thr=-(mean_first - misjoin_height_thr),
        width=config["first"]["thr_min_valley_width"],
        group_distance=config["first"]["valley_group_distance"],
    )
    # Recalculate thresholds based on means without peaks/valleys.
    collapse_height_thr = round(mean_no_peaks * config["first"]["thr_collapse_peak"])
    misjoin_height_thr = round(-mean_no_valleys * config["first"]["thr_misjoin_valley"])

    # Remove secondary rows that don't meet minimal secondary coverage.
    second_thr = round(mean_first * config["second"]["thr_min_perc_first"])
    second_outliers_coords: IntervalTree = get_secondary_allele_coords(
        df_cov,
        second_thr=second_thr,
        group_distance=config["second"]["group_distance"],
        min_group_size=config["second"]["thr_min_group_size"],
        thr_het_ratio=config["second"]["thr_het_ratio"],
    )

    classified_second_outliers: set[Interval] = set()
    misassemblies: defaultdict[Misassembly, IntervalTree] = defaultdict(IntervalTree)

    identify_zero_cov_regions(
        df_cov,
        misassemblies,
    )
    identify_collapses(
        first_peak_coords,
        second_outliers_coords,
        classified_second_outliers,
        misassemblies,
        collapse_height_thr=collapse_height_thr - mean_no_peaks,
    )
    identify_misjoins(
        first_valley_coords,
        second_outliers_coords,
        classified_second_outliers,
        misassemblies,
        misjoin_height_thr=misjoin_height_thr,
    )
    identify_hets(
        second_outliers_coords,
        classified_second_outliers,
        misassemblies,
    )

    # Annotate df with misassembly and remove misassemblies in ignored regions.
    df_cov = filter_misassembled_ignored_regions(
        df_cov, misassemblies, ignored_regions_intervals, ignored_mtypes
    )

    # TODO: false dupes
    return df_cov, misassemblies


def classify_plot_assembly(
    infile: str,
    chrom_sizes: str | None,
    output_dir: str | None,
    output_cov_dir: str | None,
    threads: int,
    contig: str,
    start: int,
    end: int,
    config: dict[str, Any],
    overlay_regions: defaultdict[int, set[Region]] | None,
    ignored_regions: set[Region],
    ignored_mtypes: set[Misassembly],
    ylim: float | int,
) -> pl.DataFrame:
    contig_name = f"{contig}_{start}-{end}"
    sys.stderr.write(f"Reading in NucFreq from region: {contig_name}\n")

    # Check file is indexed to prevent silent failure to read alignment file.
    check_indexed(infile)
    try:
        aln = pysam.AlignmentFile(infile, threads=threads)
        cov_first_second = np.flip(
            np.sort(get_coverage_by_base(aln, contig, start, end), axis=1).transpose(),
            axis=0,
        )
        df = pl.DataFrame(
            {
                "position": np.arange(start, end),
                "first": cov_first_second[0],
                "second": cov_first_second[1],
            },
            schema={"position": pl.UInt64, "first": pl.UInt32, "second": pl.UInt32},
        )
        del cov_first_second
    except ValueError:
        df = pl.read_csv(
            infile,
            separator="\t",
            has_header=True,
            dtypes={"position": pl.UInt64, "first": pl.UInt32, "second": pl.UInt32},
        )

    # Update ignored regions if relative.
    updated_ignored_regions = list(
        update_relative_ignored_regions(ignored_regions, ctg_start=start, ctg_end=end)
    )

    df_group_labeled, misassemblies = classify_misassemblies(
        df,
        config=config,
        ignored_regions=updated_ignored_regions,
        ignored_mtypes=ignored_mtypes,
    )
    del df

    if output_dir:
        _ = plot_coverage(
            df_group_labeled, misassemblies, contig, overlay_regions, ylim
        )

        sys.stderr.write(f"Plotting {contig_name}.\n")

        output_plot = os.path.join(output_dir, f"{contig_name}.png")
        plt.savefig(output_plot, dpi=PLOT_DPI, bbox_inches="tight")

    if output_cov_dir:
        sys.stderr.write(f"Writing coverage files for {contig_name}.\n")
        write_bigwig(
            contig,
            df_group_labeled,
            chrom_sizes,
            columns=["first", "second"],
            output_prefix=os.path.join(output_cov_dir, contig_name),
        )

    del df_group_labeled

    return pl.DataFrame(
        [
            (contig, int(interval.begin), int(interval.end), misasm)
            for misasm, intervals in misassemblies.items()
            for interval in intervals
        ],
        schema=["contig", "start", "stop", "misassembly"],
        orient="row",
    )
