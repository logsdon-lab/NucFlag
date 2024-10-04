import scipy.signal
import numpy as np
import polars as pl
from intervaltree import Interval, IntervalTree


def peak_finder(
    df: pl.DataFrame,
    *,
    height: int,
    width: int,
    invert: bool = False,
    group_distance: int = 5_000,
) -> IntervalTree:
    # Including zeroes will alter the calculation of valley heights as will always be the max.
    df_subset = df.filter(pl.col("first") != 0)
    data = df_subset["first"] if not invert else -df_subset["first"]
    positions = df_subset["position"]
    # Use height as first threshold.
    peaks, peak_info = scipy.signal.find_peaks(data, height=height, width=width)
    # Calculate median avoiding peaks.
    median_no_peaks = data.filter(~positions.is_in(peaks)).median()

    # Create interval tree adjusted to group distance.
    intervals = IntervalTree()
    for left_idx, right_idx, peak_ht in zip(
        peak_info["left_ips"], peak_info["right_ips"], peak_info["peak_heights"]
    ):
        left_idx, right_idx = int(left_idx), int(right_idx)
        left_pos, right_pos = positions[left_idx], positions[right_idx]
        # Calculate relative height
        peak_rel_ht = peak_ht - median_no_peaks
        intervals.add(
            Interval(
                left_pos - group_distance,
                right_pos + group_distance,
                peak_rel_ht,
            )
        )
    # Merge taking largest height.
    intervals.merge_overlaps(strict=False, data_reducer=lambda x, y: max(x, y))

    return IntervalTree(
        Interval(
            interval.begin + group_distance,
            interval.end - group_distance,
            interval.data,
        )
        for interval in intervals.iter()
    )


# https://stackoverflow.com/a/7353335
def consecutive(data, stepsize: int = 1):
    """
    Group consecutive positions allowing a maximum gap of stepsize.
    Larger stepsize groups more positions.
    """
    return np.split(data, np.where((np.diff(data) <= stepsize) == False)[0] + 1)  # noqa: E712


def filter_interval_expr(interval: Interval, *, col: str = "position") -> pl.Expr:
    return pl.col(col).is_between(interval.begin, interval.end)


def calculate_het_ratio(
    df: pl.DataFrame, interval: Interval, second_thr: float
) -> float:
    df_het = df.filter(filter_interval_expr(interval)).filter(
        pl.col("second") > second_thr
    )
    # Apply median filter with small window size to remove outliers while not oversmoothing.
    first_signal: np.ndarray = scipy.signal.medfilt(df_het["first"], 3)
    second_signal: np.ndarray = scipy.signal.medfilt(df_het["second"], 3)
    # Use max to get general amplitude of data.
    het_first_median = first_signal.max()
    het_second_median = second_signal.max()
    return het_second_median / (het_first_median + het_second_median)
