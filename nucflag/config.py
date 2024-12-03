DEF_CONFIG: dict[str, dict] = {
    "first": dict(
        thr_min_peak_width=20,
        thr_min_valley_width=20,
        thr_misjoin_valley=0.95,
        thr_collapse_peak=2.25,
        valley_group_distance=5_000,
        peak_group_distance=5_000,
    ),
    "second": dict(
        thr_min_perc_first=0.05,
        group_distance=10_000,
        thr_min_group_size=5,
        thr_het_ratio=0.2,
    ),
    "general": dict(window_size=5_000_000),
}
