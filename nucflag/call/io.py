import os
import gzip
import logging
from collections import OrderedDict, defaultdict
from typing import DefaultDict, Generator, Iterable, TextIO

import pyBigWig  # type: ignore
import numpy as np
import polars as pl
from matplotlib.colors import rgb2hex
from intervaltree import Interval  # type: ignore[import-untyped]

from .region import Action, ActionOpt, Region
from .status import generate_status_from_regions

logger = logging.getLogger(__name__)


IDENT_BREAKPOINTS = [
    85.0,
    90.0,
    95.0,
    97.5,
    98.0,
    98.5,
    98.75,
    99.0,
    99.25,
    99.5,
    99.75,
    100.0,
]
IDENT_COLORS = [
    "#4b3991",
    "#2974af",
    "#4a9da8",
    "#57b894",
    "#9dd893",
    "#e1f686",
    "#ffffb2",
    "#fdda79",
    "#fb9e4f",
    "#ee5634",
    "#c9273e",
    "#8a0033",
]


def read_bed_file(
    bed_file: TextIO,
) -> Generator[tuple[str, int, int, list[str]], None, None]:
    for line in bed_file.readlines():
        if line[0] == "#":
            continue
        chrm, start, end, *other = line.strip().split("\t")
        yield (chrm, int(start), int(end), other)


def read_regions(bed_file: TextIO, action: ActionOpt) -> Generator[Region, None, None]:
    for i, line in enumerate(read_bed_file(bed_file)):
        ctg, start, end, other = line
        try:
            desc = other[0]
        except IndexError:
            desc = None
        try:
            item_rgb = other[5]
        except IndexError:
            item_rgb = None

        # Allow item_rgb to be hexcode or not.
        if action == ActionOpt.PLOT:
            if item_rgb and not item_rgb.startswith("#"):
                rgb = tuple(int(v) / 255 for v in item_rgb.split(","))
                assert len(rgb) == 3, "Invalid RGB value."
                action_desc = rgb2hex(rgb)
            else:
                action_desc = item_rgb
        else:
            action_desc = None

        yield Region(
            name=ctg,
            region=Interval(start, end),
            desc=desc,
            action=Action(action, action_desc),
        )


def read_ignored_regions(infile: TextIO) -> DefaultDict[str, set[Region]]:
    ignored_regions: DefaultDict[str, set[Region]] = defaultdict(set)
    for region in read_regions(infile, action=ActionOpt.IGNORE):
        ignored_regions[region.name].add(region)

    return ignored_regions


def read_overlay_regions(
    infiles: Iterable[TextIO],
) -> DefaultDict[str, OrderedDict[str, set[Region]]]:
    """
    Read input overlay BED files.
    """
    overlay_regions: defaultdict[str, OrderedDict[str, set[Region]]] = defaultdict(
        OrderedDict
    )
    # 1 as zero idx reserved for mapq
    for i, bed in enumerate(infiles, 1):
        track_label = f"Track {i}"
        for region in read_regions(bed, action=ActionOpt.PLOT):
            if track_label not in overlay_regions[region.name]:
                overlay_regions[region.name][track_label] = set()

            overlay_regions[region.name][track_label].add(region)

    return overlay_regions


def write_bigwig(
    df_pileup: pl.DataFrame, chrom_lengths: str, columns: list[str], output_dir: str
):
    chrom = df_pileup["chrom"][0]
    start = max(0, df_pileup["pos"][0] - 1)
    end = df_pileup["pos"][-1]
    chrom_coords = f"{chrom}_{start}-{end}"

    if not os.path.exists(chrom_lengths):
        logging.warning(
            f"Chromosome lengths (-f) are required to generate bigWig files for {chrom}. Generating wig files."
        )
        for col in columns:
            outfile = os.path.join(output_dir, f"{chrom_coords}_{col}.wig.gz")
            col = "bin_ident" if col == "ident" else col
            with gzip.open(outfile, "wb") as fh:
                df_values = (
                    df_pileup.select(col)
                    .cast({col: pl.Float64})
                    .rename({col: f"fixedStep chrom={chrom} start={start} step=1"})
                )
                df_values.write_csv(fh, include_header=True)  # type: ignore[call-overload]
    else:
        df_chrom_lengths = pl.read_csv(
            chrom_lengths,
            separator="\t",
            has_header=False,
            new_columns=["chrom", "length"],
            columns=[0, 1],
        ).filter(pl.col("chrom") == chrom)

        header = list(df_chrom_lengths.iter_rows())
        for col in columns:
            outfile = os.path.join(output_dir, f"{chrom_coords}_{col}.bw")
            col = "bin_ident" if col == "ident" else col
            with pyBigWig.open(outfile, "w") as bw:
                # https://github.com/deeptools/pyBigWig/issues/126
                bw.addHeader(header)
                values = df_pileup[col].to_numpy().astype(np.float64)
                bw.addEntries(chrom, start, values=values, span=1, step=1)


def write_output(
    dfs_regions: list[pl.DataFrame],
    output_regions: TextIO | None,
    output_status: TextIO | None,
    *,
    status_by_region: bool,
) -> None:
    df_region = pl.concat(dfs_regions)

    if output_regions:
        # Erase file and then rewrite in sorted order.
        output_regions.seek(0)
        output_regions.truncate(0)
        df_region.unique().sort(by=["#chrom", "chromStart"]).write_csv(
            file=output_regions, include_header=True, separator="\t"
        )

    if not output_status:
        return

    if status_by_region:
        df_status = pl.concat(
            [
                generate_status_from_regions(df_region, groupby="region")
                for df_region in dfs_regions
            ]
        ).sort(by=["#chrom", "chromStart"])
    else:
        df_status = generate_status_from_regions(df_region, groupby="region")

    df_status.write_csv(file=output_status, include_header=True, separator="\t")


def read_identity_breakpoints(infile: TextIO | None) -> tuple[list[float], list[str]]:
    if not infile:
        return IDENT_BREAKPOINTS, IDENT_COLORS

    idents, colors = [], []
    for line in infile:
        ident, hexcode_color = line.strip().split("\t")
        idents.append(float(ident))
        colors.append(hexcode_color)

    return idents, colors
