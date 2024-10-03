DEF_CONFIG = {
    "first": dict(
        thr_min_peak_width=20,
        thr_min_valley_width=5,
        thr_misjoin_valley=0.8,
        thr_collapse_peak=3.0,
        valley_group_distance=1,
        peak_group_distance=1,
    ),
    "second": dict(
        thr_min_perc_first=0.05,
        group_distance=10_000,
        thr_min_group_size=5,
        thr_het_ratio=0.25,
    ),
    "gaps": dict(thr_max_allowed_gap_size=0),
    "general": dict(window_size=5_000_000),
}
