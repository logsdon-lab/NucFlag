from intervaltree import Interval
import pytest

from nucflag.region import (
    Action,
    ActionOpt,
    IgnoreOpt,
    Region,
    update_relative_ignored_regions,
)


def test_update_relative_region():
    updated_regions = sorted(
        update_relative_ignored_regions(
            [
                Region(
                    name="",
                    region=Interval(0, 50),
                    desc=None,
                    action=Action(ActionOpt.IGNORE, desc=IgnoreOpt.RELATIVE),
                ),
                Region(
                    name="",
                    region=Interval(-50, 0),
                    desc=None,
                    action=Action(ActionOpt.IGNORE, desc=IgnoreOpt.RELATIVE),
                ),
            ],
            ctg_start=10,
            ctg_end=100,
        )
    )
    expected_regions = sorted(
        [
            Region(
                name="",
                region=Interval(10, 60),
                desc=None,
                action=Action(opt=ActionOpt.IGNORE, desc=IgnoreOpt.RELATIVE),
            ),
            Region(
                name="",
                region=Interval(50, 100),
                desc=None,
                action=Action(opt=ActionOpt.IGNORE, desc=IgnoreOpt.RELATIVE),
            ),
        ]
    )

    assert updated_regions == expected_regions


def test_invalid_update_relative_region():
    # lower bound must be smaller than upper bound
    with pytest.raises(ValueError):
        _ = list(
            update_relative_ignored_regions(
                [
                    Region(
                        name="",
                        region=Interval(0, -250),
                        desc=None,
                        action=Action(ActionOpt.IGNORE, desc=IgnoreOpt.RELATIVE),
                    )
                ],
                ctg_start=0,
                ctg_end=100,
            )
        )

    with pytest.raises(ValueError):
        _ = list(
            update_relative_ignored_regions(
                [
                    Region(
                        name="",
                        region=Interval(250, 0),
                        desc=None,
                        action=Action(ActionOpt.IGNORE, desc=IgnoreOpt.RELATIVE),
                    )
                ],
                ctg_start=0,
                ctg_end=100,
            )
        )
