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
from intervaltree import Interval, IntervalTree

from .gap import identify_gaps
from .het import identify_hets
from .collapse import identify_collapses
from .misjoin import identify_misjoins
from .common import peak_finder, consecutive, filter_interval_expr
from ..utils import check_bam_indexed
from ..constants import PLOT_DPI, PROP_INCLUDE
from ..io import get_coverage_by_base
from ..misassembly import Misassembly
from ..plot import plot_coverage
from ..region import Region, update_relative_ignored_regions


def get_secondary_allele_coords(
    df_cov: pl.DataFrame, *, second_thr: int, group_distance: int, min_group_size: int
) -> IntervalTree:
    df_second_outliers = df_cov.filter(pl.col("second") > second_thr)
    second_outliers_coords: IntervalTree = IntervalTree()
    for grp in consecutive(df_second_outliers["position"], stepsize=group_distance):
        if len(grp) < min_group_size:
            continue
        second_outliers_coords.add(Interval(grp[0], grp[-1]))

    return second_outliers_coords


def filter_misassembled_ignored_regions(
    df_cov: pl.DataFrame,
    misassemblies: defaultdict[Misassembly, IntervalTree],
    ignored_regions: IntervalTree,
) -> pl.DataFrame:
    """
    Annotates input coverage dataframe with status column. Either good or misassembly type.
    Also filters misassemblies within ignored regions.
    """
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
            if ignored_regions.overlaps(region):
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

    return lf.collect()


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

    # If can't calculate metrics, means entire region is a gap.
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
    second_outliers_coords: IntervalTree = get_secondary_allele_coords(
        df_cov,
        second_thr=second_thr,
        group_distance=config["second"]["group_distance"],
        min_group_size=config["second"]["thr_min_group_size"],
    )

    classified_second_outliers: set[Interval] = set()
    misassemblies: defaultdict[Misassembly, IntervalTree] = defaultdict(IntervalTree)

    identify_collapses(
        first_peak_coords,
        second_outliers_coords,
        classified_second_outliers,
        misassemblies,
        collapse_height_thr=collapse_height_thr,
    )
    identify_gaps(
        df_cov,
        misassemblies,
        max_allowed_gap_size=config["gaps"]["thr_max_allowed_gap_size"],
    )
    identify_misjoins(
        df_cov,
        first_valley_coords,
        second_outliers_coords,
        classified_second_outliers,
        misassemblies,
        misjoin_height_thr=misjoin_height_thr,
        second_thr=second_thr,
        het_ratio_thr=config["second"]["thr_het_ratio"],
    )
    identify_hets(
        df_cov,
        second_outliers_coords,
        classified_second_outliers,
        misassemblies,
        second_thr=second_thr,
        het_ratio_thr=config["second"]["thr_het_ratio"],
    )

    # Annotate df with misassembly and remove misassemblies in ignored regions.
    df_cov = filter_misassembled_ignored_regions(
        df_cov, misassemblies, ignored_regions_intervals
    )

    # TODO: false dupes
    return df_cov, misassemblies


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

    # Check bamfile is indexed to prevent silent failure to read alignment file.
    check_bam_indexed(infile)
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
