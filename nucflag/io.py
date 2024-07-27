import sys
from collections import defaultdict
from typing import DefaultDict, Generator, Iterable, TextIO

import numpy as np
import portion as pt
import pysam

from .region import Action, ActionOpt, IgnoreOpt, Region


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
    bamfile: str, input_regions: TextIO | None, *, threads: int = 4
) -> Generator[tuple[str, int, int], None, None]:
    if input_regions:
        sys.stderr.write(f"Reading in regions from {input_regions.name}.\n")

        yield from (
            (ctg, start, stop) for ctg, start, stop, *_ in read_bed_file(input_regions)
        )
    else:
        refs = {}
        with pysam.AlignmentFile(bamfile, threads=threads) as bam:
            sys.stderr.write(f"Reading entire {bam} because no bedfile was provided.\n")
            for read in bam.fetch(until_eof=True):
                ref = read.reference_name
                if ref not in refs:
                    refs[ref] = [2147483648, 0]

                start = read.reference_start
                end = read.reference_end
                if refs[ref][0] > start:
                    refs[ref][0] = start
                if refs[ref][1] < end:
                    refs[ref][1] = end

        for contig in refs:
            yield (contig, refs[contig][0], refs[contig][1])


def read_regions(bed_file: TextIO) -> Generator[Region, None, None]:
    for line in read_bed_file(bed_file):
        ctg, start, end, other = line
        try:
            desc = other[0]
        except IndexError:
            desc = None
        try:
            actions_str = other[1]
            # Split actions column.
            # TODO: Or use multi-cols?
            for action_str in actions_str.split(","):
                action_opt, _, action_desc = action_str.partition(":")
                action_opt = ActionOpt(action_opt)
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
