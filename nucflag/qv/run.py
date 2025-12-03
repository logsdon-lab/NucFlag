import math
import argparse

import polars as pl

from ..common import BED9_COLS, add_group_columns


def calculate_qv(args: argparse.Namespace) -> int:
    df_calls = pl.read_csv(
        args.infile,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        schema=dict(BED9_COLS),
    )

    if args.ignore_calls:
        df_calls = df_calls.with_columns(
            name=pl.when(pl.col("name").is_in(args.ignore_calls))
            .then(pl.lit("correct"))
            .otherwise(pl.col("name"))
        )

    # Group regions if not bookended.
    df_calls_grouped = add_group_columns(df_region=df_calls)

    print("#chrom", "chromStart", "chromEnd", "QV", sep="\t")
    for grp, df_grp in df_calls_grouped.group_by(
        ["#chrom", "group"], maintain_order=True
    ):
        chrom = grp[0]
        min_start = df_grp["minStart"][0]
        max_end = df_grp["maxEnd"][0]

        bp_region = max_end - min_start
        bp_err = df_grp.filter(pl.col("name") != "correct")["length"].sum()

        assert (
            bp_err <= bp_region
        ), f"Length of error region ({bp_err}) exceeds length of region ({bp_region}) for {chrom}"

        try:
            # Similar to Inspector's QV metric
            # https://github.com/ChongLab/Inspector/blob/0e08f882181cc0e0e0fa749cd87fb74a278ea0f0/inspector.py#L184
            qv = abs(-10 * math.log10(bp_err / bp_region))
        except ValueError:
            # No errors in region.
            qv = math.inf

        print(chrom, min_start, max_end, qv, sep="\t")

    return 0
