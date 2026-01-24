import sys
import argparse

from typing import TYPE_CHECKING, Any

from ..common import PRESETS

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def add_config_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "config",
        description="Generate template configfile.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-x",
        "--preset",
        type=str,
        default=None,
        choices=PRESETS,
        help="Configuration preset to use.",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("wt"),
        default=sys.stdout,
        help="Output config file.",
    )

    return None
