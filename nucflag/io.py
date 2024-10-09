import sys
from collections import defaultdict
from typing import DefaultDict, Generator, Iterable, TextIO

import polars as pl
import numpy as np
import pysam
from intervaltree import Interval

from .config import DEF_CONFIG
from .utils import check_bam_indexed
from .misassembly import Misassembly
from .region import Action, ActionOpt, IgnoreOpt, Region, RegionStatus


def get_coverage_by_base(
    bam: pysam.AlignmentFile, contig: str, start: int, end: int
) -> np.ndarray:
    coverage = bam.count_coverage(
        contig, start=start, stop=end, read_callback="nofilter", quality_threshold=None
    )
    assert len(coverage) == 4
    return np.array(tuple(zip(*coverage)))


def read_bed_file(
    bed_file: TextIO,
) -> Generator[tuple[str, int, int, list[str]], None, None]:
    for line in bed_file.readlines():
        if line[0] == "#":
            continue
        chrm, start, end, *other = line.strip().split("\t")
        yield (chrm, int(start), int(end), other)


def read_asm_regions(
    infile: str,
    input_regions: TextIO | None,
    *,
    threads: int = 4,
    window_size: int = DEF_CONFIG["general"]["window_size"],
) -> Generator[tuple[str, int, int], None, None]:
    if input_regions:
        sys.stderr.write(f"Reading in regions from {input_regions.name}.\n")

        yield from (
            (ctg, start, stop) for ctg, start, stop, *_ in read_bed_file(input_regions)
        )
    else:
        check_bam_indexed(infile)
        if not infile.endswith(".bam"):
            raise NotImplementedError(
                "Reading regions from coverage file not supported."
            )

        with pysam.AlignmentFile(infile, threads=threads) as bam:
            sys.stderr.write(
                f"Reading entire {bam.filename} in {window_size:,} bp intervals because no bedfile was provided.\n"
            )
            for ref in bam.references:
                ref_len = bam.get_reference_length(ref)
                num, rem = divmod(ref_len, window_size)
                for i in range(1, num + 1):
                    yield (ref, (i - 1) * window_size, i * window_size)

                final_start = num * window_size
                yield (ref, final_start, final_start + rem)


def read_regions(bed_file: TextIO) -> Generator[Region, None, None]:
    for i, line in enumerate(read_bed_file(bed_file)):
        ctg, start, end, other = line
        try:
            desc = other[0]
        except IndexError:
            desc = None
        try:
            actions_str = other[1]
            # Split actions column.
            for action_str in actions_str.split(","):
                action_desc: IgnoreOpt | str | None
                action_opt, _, action_desc = action_str.partition(":")
                try:
                    action_opt = ActionOpt(action_opt)
                except ValueError:
                    sys.stderr.write(
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
        except IndexError:
            continue


def read_ignored_regions(infile: TextIO) -> DefaultDict[str, set[Region]]:
    ignored_regions: DefaultDict[str, set[Region]] = defaultdict(set)
    for region in read_regions(infile):
        if region.action and region.action.opt == ActionOpt.IGNORE:
            ignored_regions[region.name].add(region)

    return ignored_regions


def read_overlay_regions(
    infiles: Iterable[TextIO],
    *,
    ignored_regions: DefaultDict[str, set[Region]] | None = None,
) -> DefaultDict[str, DefaultDict[int, set[Region]]]:
    """
    Read input overlay BED files and optionally updated ignored regions if any are specified.
    """
    overlay_regions: DefaultDict[str, DefaultDict[int, set[Region]]] = defaultdict(
        lambda: defaultdict(set)
    )
    for i, bed in enumerate(infiles):
        for region in read_regions(bed):
            # Add region to ignored regions.
            if (
                isinstance(ignored_regions, defaultdict)
                and region.action
                and (region.action.opt == ActionOpt.IGNORE)
            ):
                ignored_regions[region.name].add(region)
            else:
                overlay_regions[region.name][i].add(region)

    return overlay_regions


def write_misassemblies_and_status(
    dfs_misasm: Iterable[pl.DataFrame],
    regions: Iterable[tuple[str, int, int]],
    output_misasm: TextIO,
    output_status: TextIO | None,
    *,
    bed_provided: bool,
):
    try:
        df_misasm = pl.concat(df for df in dfs_misasm if not df.is_empty())
    except ValueError:
        df_misasm = pl.DataFrame(schema=["contig", "start", "end", "misassembly"])

    df_misasm = df_misasm.with_columns(
        contig=pl.when(bed_provided)
        .then(pl.col("contig"))
        .otherwise(pl.col("contig").str.replace(r":\d+-\d+$", ""))
    )

    df_misasm.sort(by=["contig", "start"]).write_csv(
        file=output_misasm, include_header=False, separator="\t"
    )

    if output_status:
        # Save misassemblies to output bed
        region_status = []
        # If only HETs, consider correctly assembled.
        for region, start, end in regions:
            df_misasm_ctg = df_misasm.filter(
                pl.col("contig") == f"{region}:{start}-{end}"
            )
            if df_misasm_ctg.filter(
                pl.col("misassembly") != str(Misassembly.HET)
            ).is_empty():
                region_status.append((region, start, end, RegionStatus.GOOD))
            else:
                region_status.append((region, start, end, RegionStatus.MISASSEMBLED))

        df_asm_status = pl.DataFrame(
            region_status, schema=["contig", "start", "end", "status"]
        )
        df_asm_status.sort(by=["contig", "start"]).write_csv(
            output_status, include_header=False, separator="\t"
        )
