import sys
from collections import defaultdict
from typing import DefaultDict, Generator, Iterable, TextIO

import numpy as np
import portion as pt
import pysam

from .region import Action, ActionOpt, IgnoreOpt, Region
from .constants import WINDOW_SIZE


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
        chrm, start, end, *other = line.strip().split()
        yield (chrm, int(start), int(end), other)


def read_asm_regions(
    bamfile: str,
    input_regions: TextIO | None,
    *,
    threads: int = 4,
    window_size: int = WINDOW_SIZE,
) -> Generator[tuple[str, int, int], None, None]:
    if input_regions:
        sys.stderr.write(f"Reading in regions from {input_regions.name}.\n")

        yield from (
            (ctg, start, stop) for ctg, start, stop, *_ in read_bed_file(input_regions)
        )
    else:
        with pysam.AlignmentFile(bamfile, threads=threads) as bam:
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
                    region=pt.open(start, end),
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
