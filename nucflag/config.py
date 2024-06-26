DEF_CONFIG = {
    "first": dict(
        thr_min_peak_horizontal_distance=1,
        thr_min_peak_width=20,
        thr_min_valley_horizontal_distance=1,
        thr_min_valley_width=5,
        thr_peak_height_std_above=3.2,
        thr_valley_height_std_below=3,
        thr_misjoin_valley_height_perc_below=0.001,
        valley_group_distance=500,
        peak_group_distance=500,
    ),
    "second": dict(
        thr_min_perc_first=0.05,
        thr_peak_height_std_above=3,
        group_distance=30_000,
        thr_min_group_size=20,
        thr_collapse_het_ratio=0.2,
    ),
    "gaps": dict(thr_max_allowed_gap_size=0),
}
