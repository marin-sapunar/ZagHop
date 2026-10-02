"""Calculate adiabatic or diabatic populations from ZagHop trajectories and write to file."""
import argparse
import os
import sys
import numpy as np
import scipy.stats


def add_subparser(subparsers):
    parser = subparsers.add_parser(
        "population",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Calculate populations from an ensemble of ZagHop trajectories.
          Trajectories which ended early are continued until the end of the longest trajectory
          based on their Results/status file: in the ground state after an intersection with
          the ground state (4), in their final state after reaching the target state (2, 3) and
          in an additional 'unknown' state after an error or unfinished run (<= 0).""")
    parser.add_argument(
        "dirs",
        type=str,
        nargs="+",
        metavar="DIR",
        help="Paths to trajectory directories containing Results/energy.dat and Results/status.")
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
    unknown_group = parser.add_mutually_exclusive_group()
    unknown_group.add_argument(
        "--discard-unknown-traj",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Completely ignore trajectories which ended in an error.")
    unknown_group.add_argument(
        "--discard-unknown-step",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="""Include trajectories which ended in an error until their termination, then
          neglect them and normalize the population of the remaining trajectories.""")
    parser.add_argument(
        "--output",
        type=str,
        default="pop",
        metavar="NAME",
        help="Base name for the output files.")
    parser.set_defaults(func=run)


def read_status(cdir):
    """Read the termination status of a trajectory from its Results/status file.

    0: still running or ended abruptly, < 0: error, 1: reached max_time, 2: reached target
    state, 3: reached max_time after target state, 4: intersection with ground state.
    """
    sfile = os.path.join(cdir, "Results", "status")
    if not os.path.isfile(sfile):
        print(f"Warning: {sfile} not found, assuming status 0.", file=sys.stderr)
        return 0
    with open(sfile, "r") as f:
        return int(f.read().split()[0])


def run(args):
    dirs = np.array(args.dirs, dtype=str)
    status = np.array([read_status(cdir) for cdir in dirs])
    invalid = dirs[status > 4]
    if invalid.size > 0:
        raise ValueError(f"Unrecognized status in trajectories: {invalid}")
    counts = zip(*np.unique(status, return_counts=True))
    print("Trajectory status counts: " + ", ".join(f"{s}: {c}" for s, c in counts))
    unknown = status <= 0
    if args.discard_unknown_traj:
        dirs = dirs[~unknown]
        status = status[~unknown]
        unknown = unknown[~unknown]
        if not dirs:
            raise ValueError("All trajectories ended in an error.")

    times = []
    pops = []
    final_states = []
    for cdir in dirs:
        efile = os.path.join(cdir, "Results", "energy.dat")
        edata = np.loadtxt(efile, comments="#", usecols=(0, 1), ndmin=2)
        times.append(edata[:, 0])
        final_states.append(int(edata[-1, 1]))
        if args.representation == "adiabatic":
            pops.append(edata[:, 1].astype(int))
        elif args.representation == "diabatic":
            cstate_file = os.path.join(cdir, "Results", "cstate_diab")
            pops.append(np.loadtxt(cstate_file, comments="#", ndmin=2)**2)

    time = max(times, key=len)
    ntraj = len(dirs)
    ntime = time.shape[0]

    if args.representation == "adiabatic":
        if args.nstate is not None:
            nstate = sum(args.nstate * np.array([1, 2, 3]))
        else:
            nstate = max(np.max(cstate) for cstate in pops)
        pops = [cstate[:, None] == np.arange(1, nstate + 1) for cstate in pops]
    elif args.representation == "diabatic":
        nstate = pops[0].shape[1]

    # Populations of trajectories which ended early are continued until the end of the longest
    # trajectory. Unknown populations are either marked as unknown or excluded (NaN).
    cpop = np.zeros((ntraj, ntime, nstate))
    cpop_unknown = np.zeros((ntraj, ntime))
    trajs = zip(dirs, times, pops, final_states, status)
    for i, (cdir, ctime, pop, final_state, cstatus) in enumerate(trajs):
        n = min(len(ctime), len(pop), ntime)
        if not np.allclose(ctime[:n], time[:n], atol=1.0e-4):
            raise ValueError(f"Time steps in {cdir} do not match the other trajectories.")
        cpop[i, :n] = pop[:n]
        if cstatus in (2, 3):
            cpop[i, n:, final_state - 1] = 1
        elif cstatus == 4:
            cpop[i, n:, 0] = 1
        elif cstatus <= 0 and args.discard_unknown_step:
            cpop[i, n:] = np.nan
        elif cstatus <= 0:
            cpop_unknown[i, n:] = 1
        elif n < ntime:
            raise ValueError(f"Trajectory {cdir} reached max_time at {ctime[n-1]} fs, before "
                             f"the end of the longest trajectory at {time[-1]} fs.")

    if args.sum_mult:
        i0 = args.nstate[0]
        for mult, mult_ns in zip([2, 3], args.nstate[1:]):
            for i in range(1, mult):
                cpop[..., i0:i0+mult_ns] += cpop[..., i0+mult_ns*i:i0+mult_ns*(i+1)]
            i0 += mult_ns * mult
        cpop = cpop[..., :sum(args.nstate)]

    columns = ["time"] + [str(i + 1) for i in range(cpop.shape[-1])]
    if np.any(unknown) and not args.discard_unknown_step:
        cpop = np.concatenate([cpop, cpop_unknown[..., None]], axis=-1)
        columns.append("unknown")
    header = " ".join(columns)

    population = np.nanmean(cpop, axis=0)
    if args.bootstrap:
        conf = scipy.stats.bootstrap([cpop], np.nanmean, vectorized=True, batch=100, n_resamples=50)
        lbound = conf.confidence_interval.low
        ubound = conf.confidence_interval.high
    else:
        lbound = None
        ubound = None

    cdir = os.path.dirname(args.output)
    fname = os.path.basename(args.output)
    # Write population
    cfile = os.path.join(cdir, f"{fname}_{args.representation}.dat")
    np.savetxt(cfile, np.column_stack([time, population]), fmt="%.8f", header=header)
    if args.bootstrap:
        cfile = os.path.join(cdir, f"{fname}_{args.representation}_lbound.dat")
        np.savetxt(cfile, np.column_stack([time, lbound]), fmt="%.8f", header=header)
        cfile = os.path.join(cdir, f"{fname}_{args.representation}_ubound.dat")
        np.savetxt(cfile, np.column_stack([time, ubound]), fmt="%.8f", header=header)
