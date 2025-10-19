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

from .region import Action, ActionOpt, Region

logger = logging.getLogger(__name__)

BED9_COLS = [
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
]
STATUSES = [
    "good",
    "indel",
    "softclip",
    "het_mismap",
    "collapse",
    "misjoin",
    "mismatch",
    "false_dup",
    "homopolymer",
    "dinucleotide",
    "simple",
    "other",
    "scaffold",
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
        for region in read_regions(bed, action=ActionOpt.PLOT):
            if str(i) not in overlay_regions[region.name]:
                overlay_regions[region.name][str(i)] = set()

            overlay_regions[region.name][str(i)].add(region)

    return overlay_regions


def write_bigwig(
    df_pileup: pl.DataFrame, chrom_lengths: str, columns: list[str], output_dir: str
):
    chrom = df_pileup["chrom"][0]
    start = df_pileup["pos"][0]
    if not os.path.exists(chrom_lengths):
        logging.warning(
            f"Chromosome lengths are required to generate bigWig files for {chrom}. Generating wig files."
        )
        for col in columns:
            with gzip.open(
                os.path.join(output_dir, f"{chrom}_{col}.wig.gz"), "wb"
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
            outfile = os.path.join(output_dir, f"{chrom}_{col}.bw")
            with pyBigWig.open(outfile, "w") as bw:
                # https://github.com/deeptools/pyBigWig/issues/126
                bw.addHeader(header)
                values = df_pileup[col].to_numpy().astype(np.float64)
                bw.addEntries(chrom, start, values=values, span=1, step=1)


def write_output(
    output_regions: TextIO | list[pl.DataFrame],
    output_status: TextIO | None,
) -> None:
    output_cols = ["chrom", "chromStart", "chromEnd", "name"]

    # If written to stdout, read saved inputs.
    if isinstance(output_regions, list):
        try:
            df_region = pl.concat(output_regions)
        except Exception:
            # No assembly errors.
            df_region = pl.DataFrame(schema=output_cols)
    elif isinstance(output_regions, io.TextIOBase):
        # Load out-of-order file.
        try:
            df_region = pl.read_csv(
                output_regions.name,
                has_header=False,
                separator="\t",
                new_columns=output_cols,
                raise_if_empty=False,
            )
            # Erase file and then rewrite in sorted order.
            output_regions.truncate(0)
            df_region.unique().sort(by=["chrom", "chromStart"]).write_csv(
                file=output_regions, include_header=False, separator="\t"
            )
        except pl.exceptions.ShapeError:
            df_region = pl.DataFrame(schema=output_cols)
        except FileNotFoundError:
            return
    else:
        raise ValueError(f"Invalid misasm output. {output_regions}")
    if not output_status:
        return

    df_status = (
        df_region.with_columns(
            length=pl.col("chromEnd") - pl.col("chromStart"),
            ctg_length=pl.col("chromEnd").max().over("chrom")
            - pl.col("chromStart").min().over("chrom"),
            chromStart=pl.col("chromStart").min().over("chrom"),
            chromEnd=pl.col("chromEnd").max().over("chrom"),
        )
        .group_by(["chrom", "name"])
        .agg(
            chromStart=pl.col("chromStart").first(),
            chromEnd=pl.col("chromEnd").first(),
            perc=(pl.col("length").sum() / pl.col("ctg_length").first()) * 100.0,
        )
        .pivot(
            on="name",
            index=["chrom", "chromStart", "chromEnd"],
            values="perc",
            maintain_order=True,
        )
        .with_columns(
            status=pl.when(pl.col("good") == 100.0)
            .then(pl.lit("good"))
            .otherwise(pl.lit("misassembled")),
            # Ensure column exists.
            # https://github.com/pola-rs/polars/issues/18372#issuecomment-2390371173
            **{status: pl.coalesce(pl.col(f"^{status}$"), 0.0) for status in STATUSES},
        )
        .select(
            pl.col("chrom"),
            pl.col("chromStart"),
            pl.col("chromEnd"),
            pl.col("status"),
            *[pl.col(status) for status in STATUSES],
        )
        .fill_null(0.0)
        # chromStart and chromEnd get cast to float for some reason.
        .cast({"chromStart": pl.Int64, "chromEnd": pl.Int64})
    )
    df_status.sort(by="chrom").write_csv(
        file=output_status, include_header=True, separator="\t"
    )
