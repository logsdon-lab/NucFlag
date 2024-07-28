from enum import StrEnum, auto
from typing import Any, Generator, Iterable, NamedTuple

import numpy as np
import portion as pt


class RegionStatus(StrEnum):
    MISASSEMBLED = auto()
    GOOD = auto()


class IgnoreOpt(StrEnum):
    ABSOLUTE = auto()
    RELATIVE = auto()


class ActionOpt(StrEnum):
    IGNORE = auto()
    PLOT = auto()


class Action(NamedTuple):
    opt: ActionOpt
    desc: Any | None


class Region(NamedTuple):
    name: str
    region: pt.Interval
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

        if region.region.lower > region.region.upper:
            raise ValueError(
                f"Region lower bound cannot be larger than upper bound. ({region})"
            )
        if region.region.lower < 0:
            rel_start = ctg_end
        else:
            rel_start = ctg_start

        lower = rel_start + region.region.lower
        upper = rel_start + region.region.upper

        yield Region(
            region.name,
            pt.open(max(lower, 0), np.clip(upper, 0, ctg_end)),
            region.desc,
            region.action,
        )
