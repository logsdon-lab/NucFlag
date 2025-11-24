import argparse

from typing import TYPE_CHECKING, Any


if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def add_ideogram_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "ideogram",
        description="Generate ideogram from NucFlag misassemblies",
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
        "-c",
        "--cytobands",
        default=None,
        help="Cytoband file. Format: [chrom, start, end, name, btype]. See https://github.com/Balthasar-eu/pyideogram.",
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
        "--track_height",
        default=2.0,
        type=float,
        help="Height of individual track.",
    )
    ap.add_argument(
        "-o", "--output_prefix", default="ideogram", help="Output file prefix."
    )

    return None
