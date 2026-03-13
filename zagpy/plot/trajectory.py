
""" Plot trajectory energy data from a ZagHop simulation. """
import os
import numpy as np
import matplotlib.pyplot as plt
from zagpy.units import ENERGY


def add_subparser(subparsers):
    """ Register the 'trajectory' subcommand under the 'plot' subparser. """
    parser = subparsers.add_parser(
        "trajectory",
        formatter_class=__import__('argparse').ArgumentDefaultsHelpFormatter,
        description="Plot trajectory energy data from a ZagHop simulation.")
    parser.add_argument(
        "directory",
        type=str,
        metavar="DIR",
        help="Path to the directory containing Results/energy.dat.")
    parser.add_argument(
        "--energy-unit",
        type=str,
        default="eV",
        choices=ENERGY.keys(),
        #metavar="UNIT",
        help="Energy unit for the plot.")
    parser.add_argument(
        "--energy-zero",
        type=float,
        default=0.0,
        metavar="E0",
        help="Energy value (in Hartree) subtracted before unit conversion.")
    parser.add_argument(
        "--stride",
        type=int,
        default=10,
        metavar="N",
        help="Plot every N-th point for the kinetic/potential energy scatter.")
    parser.add_argument(
        "--save",
        type=str,
        default=None,
        metavar="FILE",
        help="Save the figure to FILE instead of displaying it.")
    parser.set_defaults(func=run)


def run(args):
    energy_file = os.path.join(args.directory, "Results", "energy.dat")
    traj = np.loadtxt(energy_file, comments="#").T
    jumps = np.where(traj[1, 1:] != traj[1, :-1])
    traj[2:] = (traj[2:] - args.energy_zero) * ENERGY[args.energy_unit]

    plt.plot(traj[0], traj[2], 'k')
    plt.scatter(traj[0, jumps], traj[2, jumps])

    for state in traj[4:]:
        plt.plot(traj[0], state)
    plt.scatter(traj[0][::args.stride], traj[3][::args.stride])

    plt.xlabel("Time (fs)")
    plt.ylabel(f"Energy ({args.energy_unit})")

    if args.save:
        plt.savefig(args.save)
    else:
        plt.show()
