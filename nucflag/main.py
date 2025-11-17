#!/usr/bin/env python
import time
import logging
import argparse

from .call import add_call_cli, call_assemblies

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
    # ap.add_argument("-v", "--version", action="version", version=version("nucflag"))

    sub_ap = ap.add_subparsers(dest="cmd")
    add_call_cli(sub_ap)

    args = ap.parse_args()

    if args.cmd == "call":
        return call_assemblies(args)
    else:
        raise ValueError(f"Invalid command: {args.cmd}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
