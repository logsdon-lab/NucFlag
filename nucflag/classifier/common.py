import scipy.signal
import numpy as np
import polars as pl
from intervaltree import Interval, IntervalTree
import scipy.stats


def peak_finder(
    data: pl.Series,
    positions: pl.Series,
    *,
    abs_height_thr: float,
    height_thr: float,
    width: int,
    group_distance: int = 5_000,
) -> tuple[float, IntervalTree]:
    # Use height as first threshold.
    peaks, peak_info = scipy.signal.find_peaks(data, height=abs_height_thr, width=width)
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
    intervals.merge_overlaps(strict=False, data_reducer=lambda x, y: max(x, y))

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
        # Calculate relative height.
        approx_ht = vals.max() if coef[-1] < 0 else vals.min()
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
def consecutive(data, stepsize: int = 1):
    """
    Group consecutive positions allowing a maximum gap of stepsize.
    Larger stepsize groups more positions.
    """
    return np.split(data, np.where((np.diff(data) <= stepsize) == False)[0] + 1)  # noqa: E712


def filter_interval_expr(interval: Interval, *, col: str = "position") -> pl.Expr:
    return pl.col(col).is_between(interval.begin, interval.end)


def calculate_het_ratio(interval: Interval) -> float:
    het_first_max, het_second_max = interval.data
    return het_second_max / (het_first_max + het_second_max)
