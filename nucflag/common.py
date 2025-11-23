import polars as pl


BED9_COLS = [
    ("#chrom", pl.String),
    ("chromStart", pl.UInt64),
    ("chromEnd", pl.UInt64),
    ("name", pl.String),
    ("score", pl.UInt64),
    ("strand", pl.String),
    ("thickStart", pl.UInt64),
    ("thickEnd", pl.UInt64),
    ("itemRgb", pl.String),
]
STATUSES = (
    "correct",
    "indel",
    "softclip",
    "het_mismap",
    "collapse",
    "misjoin",
    "mismatch",
    "false_dup",
    "homopolymer",
    "dinucleotide",
    "simple_repeat",
    "other_repeat",
    "scaffold",
)
