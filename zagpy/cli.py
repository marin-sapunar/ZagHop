""" Command-line interface for zagpy. """
import argparse
from zagpy.sample import wigner


def main():
    parser = argparse.ArgumentParser(
        prog="zagpy",
        description="ZagHop utilities for sampling and interface tasks.")
    subparsers = parser.add_subparsers(dest="command")

    # sample subcommand with its own subparsers
    sample_parser = subparsers.add_parser(
        "sample", help="Sample initial conditions.")
    sample_sub = sample_parser.add_subparsers(dest="method")
    wigner.add_subparser(sample_sub)

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    elif args.command == "sample" and args.method is None:
        sample_parser.print_help()
    else:
        parser.print_help()
