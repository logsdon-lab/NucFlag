import pytest
import portion as pt
from nucflag.region import Region, RegionMode


# full : 350 750
# rel : 0 250 -> 350 600
# rel : 0 -250 -> x
# rel : -250 0 -> 500, 750
# rel : -250 50 -> 500, 750
# rel : -250 -50 -> 500, 700
# rel : 250 0 -> x
# rel : 50 100 -> 400
def test_check_relative_region():
    left_region = Region(name="left", region=pt.open(0, 5), mode=RegionMode.RELATIVE)
    right_region = Region(name="right", region=pt.open(-5, 0), mode=RegionMode.RELATIVE)
    assert left_region.contains(pt.open(0, 3), full=pt.open(0, 100))
    assert not left_region.contains(pt.open(0, 7), full=pt.open(0, 100))
    assert right_region.contains(pt.open(97, 99), full=pt.open(0, 100))
    assert not right_region.contains(pt.open(0, 5), full=pt.open(0, 100))


def test_invalid_check_relative_region():
    # lower bound must be smaller than upper bound
    with pytest.raises(ValueError):
        Region(name="", region=pt.open(0, -250), mode=RegionMode.RELATIVE).contains(
            pt.open(0, 0), full=pt.open(0, 100)
        )

    with pytest.raises(ValueError):
        Region(name="", region=pt.open(250, 0), mode=RegionMode.RELATIVE).contains(
            pt.open(0, 0), full=pt.open(0, 100)
        )


def test_check_absolute_region():
    region_1 = Region(name="start", region=pt.open(1, 5), mode=RegionMode.ABSOLUTE)

    assert region_1.contains(pt.open(1, 3))
    assert not region_1.contains(pt.open(0, 5))
