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
    abs_height_thr: float,
    height_thr: float,
    width: int,
    group_distance: int,
) -> tuple[float, IntervalTree]:
    """
    Finds peaks using `scipy.signal.find_peaks`.

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
            * Each interval will contain data with the following:
                * peak height relative to the mean
                * peak approximate height
    """
    # Use height as first threshold.
    # Ensure by flooring value that valley regions with cov 1 are found.
    peaks, peak_info = scipy.signal.find_peaks(
        data, height=floor(abs_height_thr), width=width
    )
    # Calculate mean avoiding peaks.
    # Correct for starting at zero.
    adj_positions = positions - positions[0]
    mean_no_peaks = data.filter(~adj_positions.is_in(peaks)).mean()
    if not mean_no_peaks:
        raise ValueError("Cannot calculate mean without peaks.")

    # Create interval tree adjusted to group distance.
    intervals = IntervalTree()
    for left_idx, right_idx, peak_height in zip(
        peak_info["left_ips"], peak_info["right_ips"], peak_info["peak_heights"]
    ):
        left_idx, right_idx = int(np.ceil(left_idx)), int(np.floor(right_idx))
        left_pos, right_pos = positions[left_idx], positions[right_idx]

        # Calculate relative height
        peak_rel_height = abs(abs(peak_height) - abs(mean_no_peaks))
        if peak_rel_height < height_thr:
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
        # Fit second-order polynomial and predict new set of values.
        # This should reduce the effect of outliers.
        pos = np.arange(left_idx, right_idx)
        poly = np.polynomial.polynomial.Polynomial.fit(
            pos, data[left_idx:right_idx], deg=2
        )
        coef = poly.convert().coef
        vals: np.ndarray = np.polynomial.polynomial.polyval(pos, coef)

        # Calculate relative height with predicted values.
        approx_height = abs(vals.max())
        peak_rel_height = abs(approx_height - abs(mean_no_peaks))
        if peak_rel_height < height_thr:
            continue
        new_intervals.add(
            Interval(
                interval.begin + group_distance,
                interval.end - group_distance,
                (peak_rel_height, approx_height),
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
