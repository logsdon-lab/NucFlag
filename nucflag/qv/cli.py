import sys
import argparse

from typing import TYPE_CHECKING, Any


from ..common import STATUSES

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def add_qv_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "qv",
        description="Calculate QV from NucFlag misassemblies based on formula: -10 * log10(bp_err / bp_region)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-i",
        "--infile",
        required=True,
        type=argparse.FileType("rb"),
        help="Input NucFlag misassembly calls as BED9.",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Output BED4 file with QV as 4th column.",
    )
    ap.add_argument(
        "-c",
        "--ignore-calls",
        help="Ignore specific calls in QV calculation.",
        nargs="*",
        default=["scaffold"],
        type=str,
        choices=[s for s in STATUSES if s != "correct"],
    )

    return None
