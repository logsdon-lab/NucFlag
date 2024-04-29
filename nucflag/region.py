import portion as pt
import numpy as np

from enum import StrEnum, auto
from typing import NamedTuple


class RegionStatus(StrEnum):
    MISASSEMBLED = auto()
    GOOD = auto()


class RegionMode(StrEnum):
    ABSOLUTE = auto()
    RELATIVE = auto()


class Region(NamedTuple):
    name: str
    region: pt.Interval
    mode: RegionMode

    def contains(self, other: pt.Interval, *, full: pt.Interval | None = None) -> bool:
        if self.mode == RegionMode.ABSOLUTE:
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
