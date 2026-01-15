import argparse
from py_nucflag import get_config_from_preset  # type: ignore[import-untyped]


def get_config(args: argparse.Namespace) -> int:
    cfg = get_config_from_preset(args.preset)
    print(cfg, file=args.outfile)
    return 0
