import io
import os
import sys
import gzip
from collections import defaultdict
from typing import DefaultDict, Generator, Iterable, TextIO

import polars as pl
import numpy as np
import pysam
import pyBigWig
from intervaltree import Interval, IntervalTree

from .config import DEF_CONFIG
from .utils import check_indexed
from .region import Action, ActionOpt, IgnoreOpt, Region, RegionStatus


def get_coverage_by_base(
    aln: pysam.AlignmentFile, contig: str, start: int, end: int
) -> np.ndarray:
    coverage = aln.count_coverage(
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
        check_indexed(infile)
        if not (infile.endswith(".bam") or infile.endswith(".cram")):
            raise NotImplementedError(
                "Reading regions from coverage file not supported."
            )
        with pysam.AlignmentFile(infile, threads=threads) as aln:
            sys.stderr.write(
                f"Reading entire {aln.filename} in {window_size:,} bp intervals because no bedfile was provided.\n"
            )
            for ref in aln.references:
                ref_len = aln.get_reference_length(ref)
                num, rem = divmod(ref_len, window_size)
                for i in range(1, num + 1):
                    yield (ref, (i - 1) * window_size, i * window_size)

                final_start = num * window_size
                yield (ref, final_start, final_start + rem)


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


def read_ignored_regions(infile: TextIO) -> DefaultDict[str, set[Region]]:
    ignored_regions: DefaultDict[str, set[Region]] = defaultdict(set)
    for region in read_regions(infile, default_actions_str="ignore:absolute"):
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


def write_aggregated_misassemblies_and_status(
    regions: Iterable[tuple[str, int, int]],
    output_misasm: TextIO | list[pl.DataFrame],
    output_status: TextIO | None,
):
    output_cols = ["contig", "start", "end", "misassembly"]

    # If written to stdout, read saved inputs.
    if isinstance(output_misasm, list):
        try:
            df_misasm = pl.concat(output_misasm)
        except Exception:
            df_misasm = pl.DataFrame(schema=output_cols)
    elif isinstance(output_misasm, io.TextIOBase):
        # Load out-of-order file.
        try:
            df_misasm = pl.read_csv(
                output_misasm.name,
                has_header=False,
                separator="\t",
                new_columns=output_cols,
                raise_if_empty=False,
            )
            # Erase file and then rewrite in sorted order.
            output_misasm.truncate(0)
            df_misasm.unique().sort(by=["contig", "start"]).write_csv(
                file=output_misasm, include_header=False, separator="\t"
            )
        except pl.exceptions.ShapeError:
            df_misasm = pl.DataFrame(schema=output_cols)
        except FileNotFoundError:
            return
    else:
        raise ValueError(f"Invalid misasm output. {output_misasm}")

    if output_status:
        # Save misassemblies to output bed
        region_status = []
        # Create interval tree per contig ignoring HETs
        itrees_misasm: defaultdict[str, IntervalTree] = defaultdict(IntervalTree)
        for contig, start, end, _ in df_misasm.filter(
            pl.col("misassembly") != "HET"
        ).iter_rows():
            itrees_misasm[contig].add(Interval(start, end))

        for region, start, end in regions:
            if not itrees_misasm[region].overlaps(start, end):
                region_status.append((region, start, end, str(RegionStatus.GOOD)))
            else:
                region_status.append(
                    (region, start, end, str(RegionStatus.MISASSEMBLED))
                )

        df_asm_status = pl.DataFrame(
            region_status, orient="row", schema=["contig", "start", "end", "status"]
        )
        df_asm_status.sort(by=["contig", "start"]).write_csv(
            output_status, include_header=False, separator="\t"
        )


def write_bigwig(
    contig: str,
    df_pileup: pl.DataFrame,
    contig_lengths: str | None,
    columns: list[str],
    output_prefix: str,
):
    start = df_pileup["position"][0]
    if not contig_lengths or not os.path.exists(contig_lengths):
        sys.stderr.write(
            f"Chrom lengths needed to generate bigWig for {contig}. Generating wig files.\n"
        )
        for col in columns:
            outfile = f"{output_prefix}_{col}.wig.gz"
            with gzip.open(outfile, "wb") as fh:
                df_values = (
                    df_pileup.select(col)
                    .cast({col: pl.Float64})
                    .rename({col: f"fixedStep chrom={contig} start={start} step=1"})
                )
                df_values.write_csv(fh, include_header=True)
    else:
        df_contig_lengths = pl.read_csv(
            contig_lengths,
            separator="\t",
            has_header=False,
            new_columns=["contig", "length"],
            columns=[0, 1],
        ).filter(pl.col("contig") == contig)

        header = list(df_contig_lengths.iter_rows())
        for col in columns:
            outfile = f"{output_prefix}_{col}.bw"
            with pyBigWig.open(outfile, "w") as bw:
                # https://github.com/deeptools/pyBigWig/issues/126
                bw.addHeader(header)
                values = df_pileup[col].to_numpy().astype(np.float64)
                bw.addEntries(contig, start, values=values, span=1, step=1)
