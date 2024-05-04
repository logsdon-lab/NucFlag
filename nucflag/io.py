import re
import sys
import pysam
import argparse

import numpy as np
import portion as pt

from .region import Region, RegionMode
from .constants import RGX_REGION

from typing import TextIO, Generator


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


def read_regions(
    bam: pysam.AlignmentFile, args: argparse.Namespace
) -> Generator[tuple[str, int, int], None, None]:
    if args.regions is not None or args.input_regions is not None:
        if args.regions is not None:
            sys.stderr.write(f"Reading in {len(args.regions)} region(s).\n")

            for region in args.regions:
                match = re.match(RGX_REGION, region)
                assert match, region + " not valid!"
                chrm, start, end = match.groups()
                yield (chrm, int(start), int(end))

        if args.input_regions is not None:
            sys.stderr.write(f"Reading in regions from {args.input_regions.name}.\n")

            yield from (
                (ctg, start, stop)
                for ctg, start, stop, *_ in read_bed_file(args.input_regions)
            )
    else:
        refs = {}

        sys.stderr.write(
            "Reading the whole bam because no region or bed argument was made.\n"
        )
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


def read_ignored_regions(bed_file: TextIO) -> Generator[Region, None, None]:
    for i, line in enumerate(read_bed_file(bed_file)):
        ctg, start, end, other = line
        try:
            mode = RegionMode(other[0])
        except IndexError:
            sys.stderr.write(
                f"Line {i} ({line}) in {bed_file.name} doesn't have a mode. Skipping."
            )
            continue

        yield Region(name=ctg, region=pt.open(start, end), mode=mode)
