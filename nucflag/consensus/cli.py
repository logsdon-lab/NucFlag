import sys
import argparse

from typing import TYPE_CHECKING, Any

from .constants import GOOD_REGIONS

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def add_consensus_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "consensus",
        description="Generate consensus from putative misassembly bedfiles.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-i",
        "--infiles",
        required=True,
        nargs="+",
        type=argparse.FileType("rb"),
        help=f"Input misassembly calls as BED4+. Any value in the fourth column is expected to be an error with the exceptions of: {GOOD_REGIONS}",
    )
    ap.add_argument(
        "-n",
        "--number_ovl",
        type=int,
        default=2,
        help="Number of overlaps across all --infiles to be considered true.",
    )
    ap.add_argument(
        "-p",
        "--perc_ovl",
        type=float,
        default=None,
        help="Required overlap as percentage (0.8 => 80 percent overlap) to be considered true. Provide value for stricter callsets. If None, any overlap is considered true.",
    )
    ap.add_argument(
        "-t",
        "--type_ovl",
        help="Which type of overlap. Similar to bedtools: -f (Fraction of A), -F (Fraction of B), r (Reciprocal fraction A_B).",
        default="A",
        type=str,
        choices=["A", "B", "A_B"],
    )
    ap.add_argument(
        "-f",
        "--filter_calls",
        type=str,
        nargs="*",
        help=f"Calls to filter out. Typically good calls: {GOOD_REGIONS}.",
    )
    ap.add_argument(
        "-r",
        "--relative_to",
        default=None,
        type=int,
        help="Index of infile to use as reference. Defaults to file with the largest number of calls.",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("wt"),
        default=sys.stdout,
        help="Output file. Produces concensus BED4 file with overlaps delimited by ',' in 4th column.",
    )

    return None
