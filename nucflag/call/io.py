import io
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

from ..common import BED9_COLS, STATUSES, add_group_columns
from .region import Action, ActionOpt, Region

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
            f"Chromosome lengths are required to generate bigWig files for {chrom}. Generating wig files."
        )
        for col in columns:
            with gzip.open(
                os.path.join(output_dir, f"{chrom_coords}_{col}.wig.gz"), "wb"
            ) as fh:
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
            with pyBigWig.open(outfile, "w") as bw:
                # https://github.com/deeptools/pyBigWig/issues/126
                bw.addHeader(header)
                values = df_pileup[col].to_numpy().astype(np.float64)
                bw.addEntries(chrom, start, values=values, span=1, step=1)


def generate_status_from_regions(df_region: pl.DataFrame) -> pl.DataFrame:
    # Regions don't have coordinates so need to group by break in contiguity of adjacent intervals.
    df_region_grp = add_group_columns(df_region)
    return (
        df_region_grp.group_by(["#chrom", "name", "group"])
        .agg(
            chromStart=pl.col("minStart").first(),
            chromEnd=pl.col("maxEnd").first(),
            perc=(
                pl.col("length").sum()
                / (pl.col("maxEnd").first() - pl.col("minStart").first())
            )
            * 100.0,
        )
        .pivot(
            on="name",
            index=["#chrom", "chromStart", "chromEnd"],
            values="perc",
            maintain_order=True,
        )
        # Ensure column exists.
        # https://github.com/pola-rs/polars/issues/18372#issuecomment-2390371173
        .with_columns(
            **{status: pl.coalesce(pl.col(f"^{status}$"), 0.0) for status in STATUSES},
        )
        .with_columns(
            status=pl.when(pl.col("correct") == 100.0)
            .then(pl.lit("correct"))
            .otherwise(pl.lit("misassembled")),
        )
        .select(
            pl.col("#chrom"),
            pl.col("chromStart"),
            pl.col("chromEnd"),
            pl.col("status"),
            *[pl.col(status) for status in STATUSES],
        )
        .fill_null(0.0)
        # chromStart and chromEnd get cast to float for some reason.
        .cast({"chromStart": pl.Int64, "chromEnd": pl.Int64})
        .sort(by=["#chrom", "chromStart"])
    )


def write_output(
    output_regions: TextIO | list[pl.DataFrame],
    output_status: TextIO | None,
) -> None:
    # If written to stdout, read saved inputs.
    if isinstance(output_regions, list):
        try:
            df_region = pl.concat(output_regions)
        except Exception:
            # No assembly errors.
            df_region = pl.DataFrame(schema=BED9_COLS)
    elif isinstance(output_regions, io.TextIOBase):
        # Load out-of-order file.
        try:
            df_region = pl.read_csv(
                output_regions.name,
                has_header=False,
                separator="\t",
                schema=dict(BED9_COLS),
                raise_if_empty=False,
            )
            # Erase file and then rewrite in sorted order.
            output_regions.seek(0)
            output_regions.truncate(0)
            df_region.unique().sort(by=["#chrom", "chromStart"]).write_csv(
                file=output_regions, include_header=True, separator="\t"
            )
        except pl.exceptions.ShapeError:
            df_region = pl.DataFrame(schema=BED9_COLS)
        except FileNotFoundError:
            return
    else:
        raise ValueError(f"Invalid misasm output. {output_regions}")
    if not output_status:
        return

    df_status = generate_status_from_regions(df_region)
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
