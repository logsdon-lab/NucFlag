import math
import argparse

import polars as pl

from typing import Collection

from ..common import (
    STATUSES,
    BED9P_COLS,
    NON_ERROR_STATUSES,
    group_dataframe_by_contiguous_itvs,
)


def expr_qv() -> pl.Expr:
    """
    QV expression. Requires bpError and bpTotal column.
    """
    return (
        pl.when(pl.col("bpError").eq(pl.lit(0)))
        .then(pl.lit(math.inf))
        # Null interval
        .when(pl.col("bpTotal").eq(pl.lit(0)))
        .then(pl.lit(0.0))
        # Same as Inspector's QV metric.
        # https://github.com/ChongLab/Inspector/blob/0e08f882181cc0e0e0fa749cd87fb74a278ea0f0/inspector.py#L184
        # $ python -c "import math,sys; e,t = sys.argv[1:3]; print(-10 * math.log10(int(e) / int(t)))" 1 1670
        #
        # For merqury.
        # $ echo "1 1670" | awk -v k=31 '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'
        #
        # Both increase as numerator drops or denominator increases.
        .otherwise(-10 * (pl.col("bpError") / pl.col("bpTotal")).log10())
    )


def add_qv_to_df(
    df: pl.DataFrame, ignore_calls: Collection[str] = NON_ERROR_STATUSES
) -> pl.DataFrame:
    error_statuses = set(STATUSES).difference((*ignore_calls, "correct"))
    return (
        df.group_by(["#chrom", "group"], maintain_order=True)
        .agg(
            chromStart=pl.col("minStart").first(),
            chromEnd=pl.col("maxEnd").first(),
            bpTotal=pl.when(pl.col("name").is_in([*ignore_calls, "correct"]))
            .then(pl.col("length"))
            .otherwise(pl.lit(0))
            .sum(),
            bpError=pl.when(pl.col("name").is_in(error_statuses))
            .then(pl.col("length"))
            .otherwise(pl.lit(0))
            .sum(),
        )
        .with_columns(QV=expr_qv())
    )


def calculate_qv(args: argparse.Namespace) -> int:
    df_calls = pl.read_csv(
        args.infile,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        columns=list(range(9)),
        schema=dict(BED9P_COLS[0:9]),
        truncate_ragged_lines=True,
    )

    # Group regions if not bookended.
    df_calls_grouped = group_dataframe_by_contiguous_itvs(df_region=df_calls)
    # Also do all.
    df_calls_all = df_calls_grouped.with_columns(
        group=pl.lit(0),
        minStart=pl.lit(0).cast(pl.UInt64),
        maxEnd=pl.lit(0).cast(pl.UInt64),
        group_length=pl.col("group_length").unique().sum(),
        group_rows=pl.col("#chrom").count(),
        **{"#chrom": pl.lit("all")},
    )
    df_qv_grouped = add_qv_to_df(
        df_calls_grouped, ignore_calls=args.ignore_calls
    ).select("#chrom", "chromStart", "chromEnd", "QV", "bpError", "bpTotal")
    df_qv_all = add_qv_to_df(df_calls_all, ignore_calls=args.ignore_calls).select(
        "#chrom", "chromStart", "chromEnd", "QV", "bpError", "bpTotal"
    )
    df_qv = pl.concat([df_qv_all, df_qv_grouped])
    df_qv.write_csv(
        args.outfile,
        separator="\t",
        include_header=True,
    )

    return 0
