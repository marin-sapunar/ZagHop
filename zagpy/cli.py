""" Command-line interface for zagpy. """
import argparse
from zagpy.sample import state, wigner
from zagpy.plot import trajectory, population
from zagpy.ensemble import population as ensemble_population


def main():
    parser = argparse.ArgumentParser(
        prog="zagpy",
        description="ZagHop utilities for sampling and interface tasks.")
    subparsers = parser.add_subparsers(dest="command")

    # sample subcommand with its own subparsers
    sample_parser = subparsers.add_parser(
        "sample", help="Sample initial conditions.")
    sample_sub = sample_parser.add_subparsers(dest="method")
    state.add_subparser(sample_sub)
    wigner.add_subparser(sample_sub)

    # ensemble subcommand with its own subparsers
    ensemble_parser = subparsers.add_parser(
        "ensemble", help="Calculate ensemble quantities from trajectory data.")
    ensemble_sub = ensemble_parser.add_subparsers(dest="ensemble_type")
    ensemble_population.add_subparser(ensemble_sub)

    # plot subcommand with its own subparsers
    plot_parser = subparsers.add_parser(
        "plot", help="Plot trajectory data.")
    plot_sub = plot_parser.add_subparsers(dest="plot_type")
    trajectory.add_subparser(plot_sub)
    population.add_subparser(plot_sub)

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    elif args.command == "sample" and args.method is None:
        sample_parser.print_help()
    elif args.command == "ensemble" and args.ensemble_type is None:
        ensemble_parser.print_help()
    elif args.command == "plot" and args.plot_type is None:
        plot_parser.print_help()
    else:
        parser.print_help()
