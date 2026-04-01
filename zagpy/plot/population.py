
""" Plot adiabatic population from an ensemble of ZagHop trajectories. """
import argparse
import os
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt


def add_subparser(subparsers):
    """ Register the 'population' subcommand under the 'plot' subparser. """
    parser = subparsers.add_parser(
        "population",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Plot adiabatic population from an ensemble of ZagHop trajectories.")
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
        help="Population representation to plot.")
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.03,
        metavar="P",
        help="Skip plotting states whose population never exceeds this value.")
    parser.add_argument(
        "--save",
        type=str,
        default=None,
        metavar="FILE",
        help="Save the figure to FILE instead of displaying it.")
    parser.add_argument(
        "--nstate-mult",
        type=int,
        nargs=3,
        default=None,
        metavar=("NSINGLET", "NDOUBLET", "NTRIPLET"),
        help="""Number of singlet, doublet, and triplet states.
             Populations of degenerate components are summed.""")
    parser.set_defaults(func=run)


def apop_from_bool(state, axis=0):
    """Calculate adiabatic population from boolean array of state occupation."""
    norm = state.shape[axis]
    apop = np.count_nonzero(state, axis=axis) / norm
    return apop


def get_apop(cstate, bootstrap=True):
    """Calculate adiabatic population from array of state indices."""
    nstate = int(np.max(cstate)) + 1
    populated = np.array([cstate == i for i in range(nstate)])
    populated = np.moveaxis(populated, 0, -1)
    apop = apop_from_bool(populated, axis=0)
    if not bootstrap:
        return apop, None
    conf = scipy.stats.bootstrap([populated], apop_from_bool, batch=100, method='basic', axis=0)
    return apop, conf


def get_cstate_proj(cstate, data):
    """Get projection of data onto the current state."""
    cs = np.array(cstate, ndmin=1)
    cs = np.expand_dims(cs, axis=(-2, -1))
    proj = np.take_along_axis(data, cs, axis=-1).squeeze()
    return proj


def dpop_from_proj(data, axis=None):
    """Calculate diabatic population from projected data."""
    if axis is None:
        return np.sum(data**2) / len(data)
    dpop = np.sum(data**2, axis=axis) / data.shape[axis]
    return dpop


def get_dpop(cstate, adt, bootstrap=True):
    """Calculate diabatic population from array of state indices and transformation matrices."""
    proj = get_cstate_proj(cstate, adt)
    dpop = dpop_from_proj(proj, axis=0)
    if not bootstrap:
        return dpop, None
    conf = scipy.stats.bootstrap([proj], dpop_from_proj, vectorized=True, batch=100, n_resamples=50)
    return dpop, conf


def plot_population(time, pop, ax=None, conf_interval=None,
                    label=None, colors=None, threshold=0.03, **kwargs):
    if ax is None:
        ax = plt.subplots()[1]
    if label is None:
        label = [f"S$_{{{i}}}$" for i in range(pop.shape[1])]
    if colors is None:
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    line_styles = ['-', '--', '-.', ':']
    for i, state_pop in enumerate(pop.T):
        if all(state_pop < threshold):
            continue
        color = colors[i % len(colors)]
        ls = line_styles[(i // len(colors)) % len(line_styles)]
        ax.plot(time, state_pop, ls=ls, c=color, label=label[i], **kwargs)
        if conf_interval is not None:
            low = conf_interval.confidence_interval.low.T[i]
            high = conf_interval.confidence_interval.high.T[i]
            ax.fill_between(time, low, high, color=color, alpha=0.2)
    ax.legend(frameon=False)
    ax.set_xlabel("Time / fs")
    ax.set_ylabel("Population")
    ax.set_ylim([0, 1])
    return ax


def reindex_mult(state, nstate_mult):
    """Convert from raw state index to multiplicity-grouped state index."""
    cs = state
    i0 = 1
    for mult, i in zip([1, 2, 3], nstate_mult):
        if i * mult > cs:
            return i0 + (cs % i) - 1
        cs -= i * mult
        i0 += i
    raise ValueError("State index out of range for given nstate_mult.")


def make_mult_labels(n_singlets, n_doublets, n_triplets):
    """Build state labels for the grouped multiplicity states."""
    labels = [f"S$_{{{i}}}$" for i in range(n_singlets)]
    labels += [f"D$_{{{i + 1}}}$" for i in range(n_doublets)]
    labels += [f"T$_{{{i + 1}}}$" for i in range(n_triplets)]
    return labels


def run(args):
    efile = os.path.join(args.dirs[0], "Results", "energy.dat")
    time = np.loadtxt(efile, comments="#", usecols=0)
    ntime = time.shape[0]

    states = []
    for cdir in args.dirs:
        efile = os.path.join(cdir, "Results", "energy.dat")
        states.append(np.loadtxt(efile, comments="#", usecols=1, dtype=int))
    cstate = np.array(states) - 1  # shape: (ntraj, ntime)

    if args.representation == "diabatic":
        adt = []
        for cdir in args.dirs:
            adt_file = os.path.join(cdir, "Results", "adt")
            data = np.loadtxt(adt_file, comments="t")
            nstate = data.shape[1]
            adt.append(data.reshape(ntime, nstate, nstate))
        adt = np.array(adt)  # (ntraj, ntime, nstate, nstate)

    labels = None
    if args.nstate_mult is not None:
        cstate = np.vectorize(reindex_mult, excluded=[1])(cstate, args.nstate_mult)
        labels = make_mult_labels(*args.nstate_mult)
        nstate_mult = sum(args.nstate_mult)
        if args.representation == "diabatic":

            adt_mult = np.zeros((adt.shape[0], adt.shape[1], nstate_mult, nstate_mult))
            for i in range(nstate_mult):
                for j in range(nstate_mult):
                    mi = reindex_mult(i, args.nstate_mult)
                    mj = reindex_mult(j, args.nstate_mult)
                    adt_mult[..., i, j] += adt[..., mi, mj]
            adt = adt_mult


    if args.representation == "adiabatic":
        population, conf = get_apop(cstate, bootstrap=args.bootstrap)
    else:
        population, conf = get_dpop(cstate, adt, bootstrap=args.bootstrap)

    _, ax = plt.subplots()
    plot_population(time, population, ax=ax, conf_interval=conf,
                    label=labels, threshold=args.threshold)

    if args.save:
        plt.savefig(args.save)
    else:
        plt.show()
