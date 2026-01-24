import sys
import ast
import argparse

from typing import TYPE_CHECKING, Any

from ..common import BED9_COLS, STATUSES, PRESETS

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


def add_call_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "call",
        description="Use per-base read coverage to classify/plot misassemblies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    input_args = ap.add_argument_group(title="Input", description="Input files.")
    input_args.add_argument(
        "-i",
        "--infile",
        required=True,
        help="Indexed BAM or CRAM file.",
    )
    input_args.add_argument(
        "-f",
        "--fasta",
        default=None,
        help="Reference fasta. Used to bin pileup using average nucleotide identity and detect repeats.",
    )
    input_args.add_argument(
        "-b",
        "--input_regions",
        default=None,
        type=argparse.FileType("rt"),
        help="BED file with regions to check.",
    )
    input_args.add_argument(
        "--ignore_regions",
        default=None,
        type=argparse.FileType("rt"),
        help="Bed file with regions to ignore. With format: [contig, start, end]",
    )
    output_args = ap.add_argument_group(title="Outputs", description="Output files.")
    output_args.add_argument(
        "-o",
        "--output_regions",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help=f"Output bed file with checked regions. With format: {[c[0] for c in BED9_COLS]}",
    )
    output_args.add_argument(
        "-s",
        "--output_status",
        default=None,
        type=argparse.FileType("wt"),
        help="Bed file with status of contigs and percentage breakdown of each misassembly type.",
    )
    output_args.add_argument(
        "-d",
        "--output_plot_dir",
        default=None,
        help="Output plot dir.",
    )
    output_args.add_argument(
        "--output_pileup_dir",
        default=None,
        help="Output pileup dir. Generates bigWig files per region.",
    )
    output_args.add_argument(
        "--add_pileup_data",
        nargs="*",
        choices=["cov", "mismatch", "mapq", "insertion", "deletion", "softclip"],
        default=["cov", "mismatch"],
        help="Add these pileup data types as bigWigs to --output_pileup_dir.",
    )
    config_args = ap.add_argument_group(title="Config", description="Configuration.")
    config_args.add_argument(
        "-t",
        "--threads",
        default=2,
        type=int,
        help="Threads for nucflag classification.",
    )
    config_args.add_argument(
        "-p",
        "--processes",
        default=8,
        type=int,
        help="Processes for plotting.",
    )
    config_args.add_argument(
        "-x",
        "--preset",
        default=None,
        choices=PRESETS,
        help="Sequencing data specific preset.",
    )
    config_args.add_argument(
        "-c",
        "--config",
        default=None,
        type=str,
        help="Threshold/params as toml file.",
    )
    config_args.add_argument(
        "--ignore_mtypes",
        nargs="*",
        choices=[status for status in STATUSES if status != "correct"],
        help="Ignore call types from plot and output bedfile.",
    )
    plot_args = ap.add_argument_group(title="Plot", description="Plot arguments.")
    plot_args.add_argument(
        "--tracks",
        nargs="*",
        type=argparse.FileType("rt"),
        help="Add additional regions as BED4 or BED9 on top of coverage plot.",
    )
    plot_args.add_argument(
        "--overlap_tracks",
        nargs="*",
        type=argparse.FileType("rt"),
        help="Add additional regions as BED4 or BED9 overlapping coverage plot.",
    )
    plot_args.add_argument(
        "--add_builtin_tracks",
        nargs="*",
        type=str,
        choices=["mapq", "ident"],
        help="Add built-in tracks used in nucflag as overlay tracks.",
    )
    plot_args.add_argument(
        "--ident_breakpoints",
        type=argparse.FileType("rt"),
        default=None,
        help="Use TSV of breakpoint and hexcode color for identity track.",
    )
    plot_args.add_argument(
        "--overlap_calls",
        dest="overlap_calls",
        action="store_true",
        help="Overlap NucFlag calls with coverage plot instead of adding as separate track. Similar to original NucFlag < v1.0. Correct regions are not displayed.",
    )
    plot_args.add_argument(
        "--ylim",
        default=3.0,
        type=ast.literal_eval,
        help="Plot y-axis limit. If float, used as a scaling factor from mean. (ex. 3.0 is mean times 3)",
    )
    return None


def add_status_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "status",
        description="Generate status output from misassembly calls.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "-i",
        "--infile",
        type=argparse.FileType("rb"),
        required=True,
        help="Input NucFlag misassembly calls as BED9.",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Bed file with status of contigs and percentage breakdown of each misassembly type.",
    )
    return None
