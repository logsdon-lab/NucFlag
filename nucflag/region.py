from enum import StrEnum, auto
from typing import Generator, Iterable, NamedTuple

import numpy as np
from intervaltree import Interval


class RegionStatus(StrEnum):
    MISASSEMBLED = auto()
    GOOD = auto()


class IgnoreOpt(StrEnum):
    ABSOLUTE = auto()
    RELATIVE = auto()


class ActionOpt(StrEnum):
    IGNORE = auto()
    PLOT = auto()
    NOOP = auto()


class Action(NamedTuple):
    opt: ActionOpt
    desc: IgnoreOpt | str | None


class Region(NamedTuple):
    name: str
    region: Interval
    desc: str | None
    action: Action | None


def update_relative_ignored_regions(
    ignored_regions: Iterable[Region], *, ctg_start: int, ctg_end: int
) -> Generator[Region, None, None]:
    for region in ignored_regions:
        if (region.action and region.action.opt != ActionOpt.IGNORE) or (
            region.action and region.action.desc != IgnoreOpt.RELATIVE
        ):
            yield region
            continue

        if region.region.begin > region.region.end:
            raise ValueError(
                f"Region lower bound cannot be larger than upper bound. ({region})"
            )
        if region.region.begin < 0:
            rel_start = ctg_end
        else:
            rel_start = ctg_start

        lower = rel_start + region.region.begin
        upper = rel_start + region.region.end

        yield Region(
            region.name,
            Interval(max(lower, 0), np.clip(upper, 0, ctg_end)),
            region.desc,
            region.action,
        )
