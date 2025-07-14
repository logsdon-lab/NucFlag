import os
import gzip
import logging
from collections import OrderedDict, defaultdict
from typing import DefaultDict, Generator, Iterable, TextIO

import pyBigWig
import numpy as np
import polars as pl
from intervaltree import Interval

from .region import Action, ActionOpt, IgnoreOpt, Region

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
    "low_quality",
    "collapse",
    "misjoin",
    "false_dupe",
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


def read_regions(
    bed_file: TextIO, default_actions_str: str | None = None
) -> Generator[Region, None, None]:
    for i, line in enumerate(read_bed_file(bed_file)):
        ctg, start, end, other = line
        try:
            desc = other[0]
        except IndexError:
            desc = None
        try:
            actions_str: str | None = other[1]
        except IndexError:
            actions_str = default_actions_str

        if not actions_str:
            continue

        # Split actions column.
        for action_str in actions_str.split(","):
            action_desc: IgnoreOpt | str | None
            action_opt, _, action_desc = action_str.partition(":")
            try:
                action_opt = ActionOpt(action_opt)
            except ValueError:
                logging.warning(
                    f"Unknown action option ({action_opt}) on line {i} in {bed_file.name}.\n"
                )
                action_opt = ActionOpt.NOOP

            if action_opt == ActionOpt.IGNORE:
                action_desc = IgnoreOpt(action_desc)
            elif action_opt == ActionOpt.PLOT:
                action_desc = action_desc
            else:
                action_desc = None

            yield Region(
                name=ctg,
                region=Interval(start, end),
                desc=desc,
                action=Action(action_opt, action_desc),
            )


def read_ignored_regions(infile: TextIO) -> DefaultDict[str, set[Region]]:
    ignored_regions: DefaultDict[str, set[Region]] = defaultdict(set)
    for region in read_regions(infile, default_actions_str="ignore:absolute"):
        if region.action and region.action.opt == ActionOpt.IGNORE:
            ignored_regions[region.name].add(region)

    return ignored_regions


def read_overlay_regions(
    infiles: Iterable[TextIO],
    *,
    ignored_regions: set[Region] | None = None,
) -> DefaultDict[str, OrderedDict[str, set[Region]]]:
    """
    Read input overlay BED files and optionally updated ignored regions if any are specified.
    """
    overlay_regions: defaultdict[str, OrderedDict[str, set[Region]]] = defaultdict(
        OrderedDict
    )
    # 1 as zero idx reserved for mapq
    for i, bed in enumerate(infiles, 1):
        for region in read_regions(bed):
            # Add region to ignored regions.
            if (
                isinstance(ignored_regions, set)
                and region.action
                and (region.action.opt == ActionOpt.IGNORE)
            ):
                ignored_regions.add(region)
            else:
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
                df_values.write_csv(fh, include_header=True)
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
    dfs_regions: list[pl.DataFrame],
    output_regions: TextIO,
    output_status: TextIO | None,
) -> None:
    try:
        df_region = pl.concat(df for df in dfs_regions if not df.is_empty())
    except ValueError:
        df_region = pl.DataFrame(schema=["chrom", "chromStart", "chromEnd", "name"])

    df_region.sort(by=["chrom", "chromStart"]).write_csv(
        file=output_regions, include_header=False, separator="\t"
    )

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
