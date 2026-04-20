"""Calculate adiabatic or diabatic populations from ZagHop trajectories and write to file."""
import argparse
import os
import numpy as np
import scipy.stats


def add_subparser(subparsers):
    parser = subparsers.add_parser(
        "population",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Calculate populations from an ensemble of ZagHop trajectories.")
    parser.add_argument(
        "dirs",
        type=str,
        nargs="+",
        metavar="DIR",
        help="Paths to completed trajectory directories containing Results/energy.dat.")
    parser.add_argument(
        "--bootstrap",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Compute bootstrap confidence intervals.")
    parser.add_argument(
        "--representation",
        type=str,
        default="adiabatic",
        choices=["adiabatic", "diabatic"],
        help="Population representation to calculate.")
    parser.add_argument(
        "--nstate",
        type=int,
        nargs=3,
        default=None,
        metavar=("NSINGLET", "NDOUBLET", "NTRIPLET"),
        help="Number of singlet, doublet, and triplet states. Needed for --sum-mult option.")
    parser.add_argument(
        "--sum-mult",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Sum populations of states with the same multiplicity.")
    parser.add_argument(
        "--output",
        type=str,
        default="pop",
        metavar="NAME",
        help="Base name for the output files.")
    parser.set_defaults(func=run)

def run(args):
    efile = os.path.join(args.dirs[0], "Results", "energy.dat")
    time = np.loadtxt(efile, comments="#", usecols=0)
    ntraj = len(args.dirs)
    ntime = time.shape[0]

    if args.representation == "adiabatic":
        states = []
        for cdir in args.dirs:
            efile = os.path.join(cdir, "Results", "energy.dat")
            states.append(np.loadtxt(efile, comments="#", usecols=1, dtype=int))
        cstate = np.array(states)  # shape: (ntraj, ntime)
        if args.nstate is not None:
            nstate = sum(args.nstate * np.array([1, 2, 3]))
        else:
            nstate = np.max(cstate)
        cpop = np.zeros((ntraj, ntime, nstate))
        for i in range(nstate):
            cpop[..., i] = cstate == i + 1
    elif args.representation == "diabatic":
        cpop = []
        for cfile in args.dirs:
            cstate_file = os.path.join(cfile, "Results", "cstate_diab")
            cpop.append(np.loadtxt(cstate_file, comments="#")**2)
        cpop = np.array(cpop)

    if args.sum_mult:
        i0 = args.nstate[0]
        for mult, mult_ns in zip([2, 3], args.nstate[1:]):
            for i in range(1, mult):
                cpop[..., i0:i0+mult_ns] += cpop[..., i0+mult_ns*i:i0+mult_ns*(i+1)]
            i0 += mult_ns * mult
        cpop = cpop[..., :sum(args.nstate)]

    population = cpop.mean(axis=0)
    if args.bootstrap:
        conf = scipy.stats.bootstrap([cpop], np.mean, vectorized=True, batch=100, n_resamples=50)
        lbound = conf.confidence_interval.low
        ubound = conf.confidence_interval.high
    else:
        lbound = None
        ubound = None

    cdir = os.path.dirname(args.output)
    fname = os.path.basename(args.output)
    # Write population
    cfile = os.path.join(cdir, f"{fname}_{args.representation}.dat")
    np.savetxt(cfile, np.column_stack([time, population]), fmt="%.8f")
    if args.bootstrap:
        cfile = os.path.join(cdir, f"{fname}_{args.representation}_lbound.dat")
        np.savetxt(cfile, np.column_stack([time, lbound]), fmt="%.8f")
        cfile = os.path.join(cdir, f"{fname}_{args.representation}_ubound.dat")
        np.savetxt(cfile, np.column_stack([time, ubound]), fmt="%.8f")
