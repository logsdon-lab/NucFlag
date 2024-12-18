from math import floor
from typing import Generator

import scipy.signal
import scipy.stats
import numpy as np
import polars as pl

from intervaltree import Interval, IntervalTree


def peak_finder(
    data: pl.Series,
    positions: pl.Series,
    *,
    abs_height_thr: float,
    height_thr: float,
    width: int,
    group_distance: int = 5_000,
) -> tuple[float, IntervalTree]:
    """
    Finds peaks using `scipy.signal.find_peaks` and fits a 2nd order polynomial to estimate peak height.

    # Arguments
    * `data`
    * `positions`
        * Positions of data.
    * `abs_height_thr`
        * Threshold absolute height of peak within data.
    * `height_thr`
        * Threshold estimated height of peak.
    * `width`
        * Threshold minimal width of peak.
    * `group_distance`
        * Group peaks by some distance.

    # Returns
    * Tuple of mean without peaks and interval tree of peaks.
    """
    # Use height as first threshold.
    # Ensure by flooring value that valley regions with cov 1 are found.
    peaks, peak_info = scipy.signal.find_peaks(
        data, height=floor(abs_height_thr), width=width
    )
    # Calculate mean avoiding peaks.
    mean_no_peaks = data.filter(~positions.is_in(peaks)).mean()
    detect_valleys = data.mean() < 0
    # Need to use adj mean for valleys as true height calculated by subtracting min coverage.
    mean_adj_no_peaks = mean_no_peaks if not detect_valleys else -data.min()
    # Create interval tree adjusted to group distance.
    intervals = IntervalTree()
    for left_idx, right_idx, peak_ht in zip(
        peak_info["left_ips"], peak_info["right_ips"], peak_info["peak_heights"]
    ):
        left_idx, right_idx = int(np.ceil(left_idx)), int(np.floor(right_idx))
        left_pos, right_pos = positions[left_idx], positions[right_idx]

        # Calculate relative height
        if detect_valleys:
            peak_rel_ht = -(mean_no_peaks + peak_ht)
        else:
            peak_rel_ht = peak_ht - mean_no_peaks

        if peak_rel_ht < height_thr:
            continue

        peak = Interval(
            left_pos - group_distance,
            right_pos + group_distance,
            (left_idx, right_idx),
        )
        intervals.add(peak)

    # Merge taking largest height and number of positives above threshold.
    intervals.merge_overlaps(
        strict=False, data_reducer=lambda x, y: (min(x[0], y[0]), max(x[1], y[1]))
    )

    new_intervals = IntervalTree()
    for interval in intervals.iter():
        left_idx, right_idx = interval.data
        # Fit second-order polynomial and predict new set of values. This should remove outliers.
        pos = np.arange(left_idx, right_idx)
        poly = np.polynomial.polynomial.Polynomial.fit(
            pos, data[left_idx:right_idx], deg=2
        )
        coef = poly.convert().coef
        vals: np.ndarray = np.polynomial.polynomial.polyval(pos, coef)
        # Calculate relative height with predicted values.
        # Take max or min depending on polynomial ort.
        # Negative a: /\
        # Positive a: \/
        # If linear: Just take max value.
        try:
            approx_ht = vals.max() if coef[2] < 0 else vals.min()
        except IndexError:
            approx_ht = vals.max()
        peak_rel_ht = abs(approx_ht - mean_adj_no_peaks)
        if peak_rel_ht < height_thr:
            continue

        new_intervals.add(
            Interval(
                interval.begin + group_distance,
                interval.end - group_distance,
                peak_rel_ht,
            )
        )

    return mean_no_peaks, new_intervals


# https://stackoverflow.com/a/7353335
def consecutive(data, stepsize: int = 1) -> list[np.ndarray]:
    """
    Group consecutive positions allowing a maximum gap of stepsize.

    # Arguments
    * `data`
        * Array of positions
    * `stepsize`
        * Split by distance.
        * Larger stepsize groups more positions.

    # Returns
    * Split arrays
    """
    return np.split(data, np.where((np.diff(data) <= stepsize) == False)[0] + 1)  # noqa: E712


def filter_interval_expr(interval: Interval, *, col: str = "position") -> pl.Expr:
    """
    Create `pl.Expr` to filter `col` by interval.

    # Arguments
    * `interval`
        * `intervaltree.Interval` to filter by.
    * `col`
        * Column name.

    # Returns
    * `pl.Expr` to filter column.
    """
    return pl.col(col).is_between(interval.begin, interval.end)


def subtract_interval(
    interval: Interval, by: set[Interval]
) -> Generator[Interval, None, None]:
    """
    Subtract `by` from `interval`.

    # Arguments
    * `interval`
    * `by`
        * Intervals to subtract from `interval`

    # Returns
    * Iterator of intervals within `interval` that aren't in `by`.
    """
    intervals = IntervalTree(Interval(i.begin, i.end) for i in by)
    intervals.add(Interval(interval.begin, interval.end))
    intervals.split_overlaps()

    for interval in intervals.iter():
        if interval in by:
            continue
        yield interval
