import portion as pt
import numpy as np

from enum import StrEnum, auto
from typing import NamedTuple, Any


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

    def contains(self, other: pt.Interval, *, full: pt.Interval | None = None) -> bool:
        if not self.action or (self.action and self.action.opt != ActionOpt.IGNORE):
            return other in self.region

        ignore_opt = self.action.desc

        if ignore_opt == IgnoreOpt.ABSOLUTE or ignore_opt is None:
            region = self.region
        else:
            if not full:
                raise ValueError("Missing full interval.")

            if self.region.lower > self.region.upper:
                raise ValueError(
                    "Region lower bound cannot be larger than upper bound."
                )
            if self.region.lower < 0:
                start = full.upper
            else:
                start = full.lower

            lower = start + self.region.lower
            upper = start + self.region.upper

            region = pt.open(max(lower, 0), np.clip(upper, 0, full.upper))

        return other in region
