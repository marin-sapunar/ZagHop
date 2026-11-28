"""Calculate electronic state populations from ZagHop trajectories and write to file."""
import argparse
import json
import os
import sys
import numpy as np
import scipy.stats
from zagpy.rng import add_seed_argument, make_rng


# Files of the character and quantum estimators for each representation (ZagHop print options
# 31, 32 and 41-43).
ESTIMATOR_FILES = {
    "character": {"adiabatic": "cstate_adiab", "diabatic": "cstate_diab"},
    "quantum": {"adiabatic": "qpop_adiab", "diabatic": "qpop_diab", "mch": "qpop_mch"},
}


def add_subparser(subparsers):
    parser = subparsers.add_parser(
        "population",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Calculate populations from an ensemble of ZagHop trajectories.
          Estimators: 'active', the fraction of trajectories in each state, from the active state
          in Results/energy.dat, in the basis of the dynamics (dynamics_basis in
          Results/states.json); 'character', the average squared character of the active state
          in the states of --representation (Results/cstate_<basis>, print options 31, 32);
          'quantum', the average populations |c|^2 from the coefficients of the electronic wave
          function in the states of --representation (Results/qpop_<basis>, print options
          41-43). Output: <output>_<estimator>_<basis>.dat.
          The adiabatic, diabatic and mch bases have a state for each spin component (the
          populations of the components of the mch basis depend on the choice of the
          components; their sums, --sum-mult, do not). Dynamics with one state for each
          multiplet (LZSH) runs in the mch_sf basis, and its 'character' populations are those
          of the diabatic multiplets (diabatic_sf).
          Trajectories which ended early are continued until the end of the longest trajectory
          based on their Results/status file. With the 'active' estimator: in the ground state
          after an intersection with the ground state (4) and in their final state after reaching
          the target state (2, 3). With 'character' and 'quantum': with their last row (2, 3,
          4); after an intersection with the ground state this is the state at the
          intersection, not the ground state. All estimators: in an additional 'unknown' state
          after an error or unfinished run (<= 0).""")
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
        "--estimator",
        type=str,
        default="active",
        choices=["active", "character", "quantum"],
        help="Population estimator.")
    parser.add_argument(
        "--representation",
        type=str,
        default=None,
        choices=["adiabatic", "diabatic", "mch"],
        help="""Basis of the populations. Required for the 'character' (adiabatic, diabatic) and
          'quantum' (adiabatic, diabatic, mch) estimators. Not accepted for 'active', whose
          basis is that of the dynamics.""")
    parser.add_argument(
        "--nstate",
        type=int,
        nargs=3,
        default=None,
        metavar=("NSINGLET", "NDOUBLET", "NTRIPLET"),
        help="""Number of singlet, doublet, and triplet states. Needed for --sum-mult option, and
          for the number of states of the 'active' estimator without Results/states.json.""")
    parser.add_argument(
        "--sum-mult",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="""Sum the populations of the spin components of each state. Not available in the
          mch_sf and diabatic_sf bases, which have one state for each multiplet.""")
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
    add_seed_argument(parser)
    parser.set_defaults(func=run)


def check_options(args):
    """Check the combination of --estimator, --representation and --sum-mult."""
    if args.estimator == "active":
        if args.representation is not None:
            raise ValueError("--representation needs --estimator character or quantum. The "
                             "'active' estimator uses the basis of the dynamics.")
        return
    if args.representation is None:
        raise ValueError(f"--estimator {args.estimator} needs --representation.")
    if args.representation not in ESTIMATOR_FILES[args.estimator]:
        choices = ", ".join(ESTIMATOR_FILES[args.estimator])
        raise ValueError(f"--estimator {args.estimator} is not available in the "
                         f"{args.representation} basis (available: {choices}).")


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


def read_layout(cdir):
    """Basis of the dynamics and its number of states from Results/states.json of a trajectory.

    Up to format version 2 of the file, the basis with one state for each multiplet was called
    'mch'; it is returned under its present name, 'mch_sf'. Returns (None, None) if the file
    does not exist (runs before states.json was written)."""
    lfile = os.path.join(cdir, "Results", "states.json")
    if not os.path.isfile(lfile):
        return None, None
    with open(lfile, "r") as f:
        layout = json.load(f)
    basis = layout["dynamics_basis"]
    nstate = layout["n_state"][basis]
    if layout.get("format_version", 1) <= 2 and basis == "mch":
        basis = "mch_sf"
    return basis, nstate


def spin_free_dynamics(dirs):
    """Whether the dynamics of the trajectories ran in a basis with one state for each multiplet.

    Trajectories without Results/states.json are not counted. Stops if only some of the others
    did."""
    bases = {read_layout(cdir)[0] for cdir in dirs} - {None}
    spin_free = {basis.endswith("_sf") for basis in bases}
    if len(spin_free) > 1:
        raise ValueError(f"The trajectories have different bases of the dynamics: {bases}.")
    return spin_free == {True}


def dynamics_basis(dirs):
    """Common basis of the dynamics of the trajectories and its number of states.

    Stops if the trajectories have different bases. Without Results/states.json in a trajectory,
    a warning is printed, and the basis is that of the other trajectories, or 'dynamics' (with
    unknown number of states) if no trajectory has the file."""
    layouts = {cdir: read_layout(cdir) for cdir in dirs}
    missing = [cdir for cdir, (basis, _) in layouts.items() if basis is None]
    known = {layout for layout in layouts.values() if layout[0] is not None}
    if len(known) > 1:
        raise ValueError(f"The trajectories have different bases of the dynamics: {known}.")
    if missing:
        print(f"Warning: Results/states.json not found in {len(missing)} trajectories (first: "
              f"{missing[0]}). Their active states are assumed to be in the basis of the "
              "dynamics of the other trajectories.", file=sys.stderr)
    if known:
        return known.pop()
    return "dynamics", None


def read_estimator(cdir, estimator, representation, n_time, complete):
    """Populations of one trajectory (n_row, n_state) from the file of the estimator.

    The rows follow those of Results/energy.dat (n_time rows). The numbers of rows must match
    if the trajectory is complete (status > 0); otherwise the files can end at different rows."""
    fname = ESTIMATOR_FILES[estimator][representation]
    pop = np.loadtxt(os.path.join(cdir, "Results", fname), comments="#", ndmin=2)
    if estimator == "character":
        pop = pop**2
    if complete and len(pop) != n_time:
        raise ValueError(f"{cdir}/Results/{fname} has {len(pop)} rows, Results/energy.dat has "
                         f"{n_time}.")
    return pop


def run(args):
    check_options(args)
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

    if args.estimator == "active":
        basis, nstate = dynamics_basis(dirs)
    else:
        basis = args.representation
        # The character of a state of a multiplet basis is given in the diabatic multiplets.
        if args.estimator == "character" and spin_free_dynamics(dirs):
            basis += "_sf"
    if args.sum_mult:
        if basis.endswith("_sf"):
            raise ValueError(f"--sum-mult is not available in the {basis} basis, which has one "
                             "state for each multiplet.")
        if args.nstate is None:
            raise ValueError("--sum-mult needs --nstate.")

    times = []
    pops = []
    final_states = []
    for cdir, cstatus in zip(dirs, status):
        efile = os.path.join(cdir, "Results", "energy.dat")
        edata = np.loadtxt(efile, comments="#", usecols=(0, 1), ndmin=2)
        times.append(edata[:, 0])
        final_states.append(int(edata[-1, 1]))
        if args.estimator == "active":
            pops.append(edata[:, 1].astype(int))
        else:
            pops.append(read_estimator(cdir, args.estimator, args.representation, len(edata),
                                       cstatus > 0))

    time = max(times, key=len)
    ntraj = len(dirs)
    ntime = time.shape[0]

    if args.estimator == "active":
        if nstate is None and args.nstate is not None:
            nstate = sum(args.nstate * np.array([1, 2, 3]))
        elif nstate is None:
            nstate = max(np.max(cstate) for cstate in pops)
        if max(np.max(cstate) for cstate in pops) > nstate:
            raise ValueError(f"Active state larger than the number of states ({nstate}).")
        pops = [cstate[:, None] == np.arange(1, nstate + 1) for cstate in pops]
    else:
        nstate = pops[0].shape[1]
        if any(pop.shape[1] != nstate for pop in pops):
            raise ValueError("The trajectories have different numbers of states.")
    if args.sum_mult and nstate != sum(args.nstate * np.array([1, 2, 3])):
        raise ValueError(f"--nstate {' '.join(map(str, args.nstate))} gives "
                         f"{sum(args.nstate * np.array([1, 2, 3]))} states with their spin "
                         f"components, the populations have {nstate}.")

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
        if cstatus in (2, 3, 4) and args.estimator != "active":
            cpop[i, n:] = pop[n - 1]
        elif cstatus in (2, 3):
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
        conf = scipy.stats.bootstrap([cpop], np.nanmean, vectorized=True, batch=100,
                                     n_resamples=50, rng=make_rng(args.seed))
        lbound = conf.confidence_interval.low
        ubound = conf.confidence_interval.high
    else:
        lbound = None
        ubound = None

    cdir = os.path.dirname(args.output)
    fname = f"{os.path.basename(args.output)}_{args.estimator}_{basis}"
    # Write population
    cfile = os.path.join(cdir, f"{fname}.dat")
    np.savetxt(cfile, np.column_stack([time, population]), fmt="%.8f", header=header)
    if args.bootstrap:
        cfile = os.path.join(cdir, f"{fname}_lbound.dat")
        np.savetxt(cfile, np.column_stack([time, lbound]), fmt="%.8f", header=header)
        cfile = os.path.join(cdir, f"{fname}_ubound.dat")
        np.savetxt(cfile, np.column_stack([time, ubound]), fmt="%.8f", header=header)
