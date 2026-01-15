#!/usr/bin/env python
import sys
import time
import logging
import argparse

from importlib.metadata import version

from .call import add_call_cli, add_status_cli, call_misassemblies, create_status
from .ideogram import add_ideogram_cli, create_ideogram
from .breakdown import add_breakdown_cli, create_breakdown_plot
from .qv import add_qv_cli, calculate_qv
from .consensus import add_consensus_cli, get_consensus_calls
from .config import add_config_cli, get_config

# Configure logging format to match rs-nucflag
# Set UTC
logging.Formatter.converter = time.gmtime
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s.%(msecs)03dZ \033[32m%(levelname)s\033[0m  [py_nucflag::%(name)s] %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)

# Create the logger
logger = logging.getLogger(__name__)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Use per-base read coverage to classify/plot misassemblies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-v", "--version", action="version", version=version("nucflag"))

    sub_ap = ap.add_subparsers(dest="cmd")
    add_call_cli(sub_ap)
    add_status_cli(sub_ap)
    add_ideogram_cli(sub_ap)
    add_breakdown_cli(sub_ap)
    add_qv_cli(sub_ap)
    add_consensus_cli(sub_ap)
    add_config_cli(sub_ap)

    args = ap.parse_args()

    if args.cmd == "call":
        return call_misassemblies(args)
    elif args.cmd == "status":
        return create_status(args)
    elif args.cmd == "ideogram":
        return create_ideogram(args)
    elif args.cmd == "breakdown":
        return create_breakdown_plot(args)
    elif args.cmd == "qv":
        return calculate_qv(args)
    elif args.cmd == "consensus":
        return get_consensus_calls(args)
    elif args.cmd == "config":
        return get_config(args)
    else:
        ap.print_help(sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
