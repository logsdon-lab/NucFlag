import scipy.signal
import numpy as np
import polars as pl
from intervaltree import Interval, IntervalTree
import scipy.stats


def peak_finder(
    data: pl.Series,
    positions: pl.Series,
    *,
    abs_height_thr: int,
    height_thr: int,
    width: int,
    group_distance: int = 5_000,
) -> tuple[float, IntervalTree]:
    # Use height as first threshold.
    peaks, peak_info = scipy.signal.find_peaks(data, height=abs_height_thr, width=width)
    # Calculate median avoiding peaks.
    mean_no_peaks = data.filter(~positions.is_in(peaks)).mean()
    # Create interval tree adjusted to group distance.
    intervals = IntervalTree()
    for left_idx, right_idx, peak_ht in zip(
        peak_info["left_ips"], peak_info["right_ips"], peak_info["peak_heights"]
    ):
        left_idx, right_idx = int(left_idx), int(right_idx)
        left_pos, right_pos = positions[left_idx], positions[right_idx]

        # Calculate change is y along peak.
        dy = np.diff(data[left_idx:right_idx])

        # Calculate relative height
        peak_rel_ht = peak_ht - mean_no_peaks
        if peak_rel_ht < height_thr:
            continue

        peak = Interval(
            left_pos - group_distance,
            right_pos + group_distance,
            (peak_rel_ht, dy.max()),
        )
        intervals.add(peak)

    # Merge taking largest height and number of positives above threshold.
    intervals.merge_overlaps(
        strict=False, data_reducer=lambda x, y: (max(x[0], y[0]), max(x[1], y[1]))
    )

    return mean_no_peaks, IntervalTree(
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
    het_first_max = first_signal.max()
    het_second_max = second_signal.max()
    return het_second_max / (het_first_max + het_second_max)
