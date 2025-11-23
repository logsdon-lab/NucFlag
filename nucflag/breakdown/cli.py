import argparse

from typing import TYPE_CHECKING, Any


if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def add_breakdown_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "breakdown",
        description="Generate bar plot breakdown of whole-genome NucFlag misassemblies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-i",
        "--infile",
        required=True,
        type=argparse.FileType("rb"),
        help="Input whole-genome NucFlag misassembly calls as BED9.",
    )
    ap.add_argument(
        "-l",
        "--filter_length",
        type=int,
        default=10_000_000,
        help="Only plot contigs at least this long.",
    )
    ap.add_argument(
        "-t",
        "--type",
        default="percent",
        choices=["percent", "length"],
        help="Type of plot. Either percent or length breakdown.",
    )
    ap.add_argument(
        "-o", "--output_prefix", default="breakdown", help="Output file prefix."
    )

    return None
