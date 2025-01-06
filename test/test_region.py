from intervaltree import Interval
import pytest

from nucflag.classifier.false_duplication import merge_itvs
from nucflag.region import (
    Action,
    ActionOpt,
    IgnoreOpt,
    Region,
    update_relative_ignored_regions,
)


@pytest.mark.parametrize(
    ["itvs", "exp_itvs", "kwargs"],
    [
        # Check book-ended merging.
        (
            [
                Interval(1, 3),
                Interval(3, 6),
                Interval(6, 9),
            ],
            [Interval(1, 9)],
            {},
        ),
        # Check dst.
        (
            [
                Interval(1, 5),
                Interval(10, 15),
                Interval(25, 30),
            ],
            [Interval(1, 15), Interval(25, 30)],
            {"dst": 5},
        ),
        # Check fn_merge_itv
        (
            [
                Interval(1, 3, 1),
                Interval(3, 6, 2),
                Interval(6, 9, 4),
            ],
            # Sum itv.data. Add prev itv.data at every merge.
            [Interval(1, 9, 8)],
            {
                "fn_merge_itv": lambda i1, i2: Interval(
                    i1.begin, i2.end, i1.data * i2.data
                )
            },
        ),
        # check fn_cmp
        (
            [
                Interval(1, 5, 1),
                Interval(5, 10, 2),
                # Cannot merge as itv.data would exceed 5.
                Interval(10, 150, 3),
            ],
            [Interval(1, 10, 3), Interval(10, 150, 3)],
            {
                "fn_merge_itv": lambda i1, i2: Interval(
                    i1.begin, i2.end, i1.data + i2.data
                ),
                "fn_cmp": lambda i1, i2: (i1.data + i2.data) <= 5,
            },
        ),
    ],
)
def test_merge_itvs(itvs, exp_itvs, kwargs):
    res = merge_itvs(itvs, **kwargs)
    assert res == exp_itvs


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
