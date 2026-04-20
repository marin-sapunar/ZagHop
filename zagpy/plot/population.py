
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
        "--nstate",
        type=int,
        nargs=3,
        default=None,
        metavar=("NSINGLET", "NDOUBLET", "NTRIPLET"),
        help="""Number of singlet, doublet, and triplet states.""")
    parser.add_argument(
        "--sum-mult",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Sum populations of states with the same multiplicity.")
    parser.set_defaults(func=run)


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
        nstate = np.max(cstate)
        cpop = np.zeros((ntraj, ntime, nstate))
        for i in range(nstate):
            cpop[..., i] = cstate == i + 1
    elif args.representation == "diabatic":
        cpop = []
        for cfile in args.dirs:
            cstate_file = os.path.join(cfile, "Results", "cstate_diab")
            cpop.append(np.loadtxt(cstate_file, comments="#"))
        cpop = np.array(cpop)**2

    if args.sum_mult:
        i0 = args.nstate[0]
        for mult, mult_ns in zip([2, 3], args.nstate[1:]):
            for i in range(1, mult):
                cpop[..., i0:i0+mult_ns] += cpop[..., i0+mult_ns*i:i0+mult_ns*(i+1)]
            i0 += mult_ns * mult
        cpop = cpop[..., :sum(args.nstate)]
    if args.nstate is not None:
        if args.sum_mult:
            labels = [f"S$_{{{i}}}$" for i in range(args.nstate[0])]
            labels += [f"D$_{{{i + 1}}}$" for i in range(args.nstate[1])]
            labels += [f"T$_{{{i + 1}}}$" for i in range(args.nstate[2])]
        else:
            labels = [f"S$_{{{i},0}}$" for i in range(args.nstate[0])]
            labels += [f"D$_{{{i + 1},1/2}}$" for i in range(args.nstate[1])]
            labels += [f"D$_{{{i + 1},-1/2}}$" for i in range(args.nstate[1])]
            labels += [f"T$_{{{i + 1},0}}$" for i in range(args.nstate[2])]
            labels += [f"T$_{{{i + 1},1}}$" for i in range(args.nstate[2])]
            labels += [f"T$_{{{i + 1},-1}}$" for i in range(args.nstate[2])]

    else:
        labels = [f"S$_{{{i}}}$" for i in range(1, cpop.shape[-1] + 1)]

    population = cpop.mean(axis=0)
    if args.bootstrap:
        conf = scipy.stats.bootstrap([cpop], np.mean, vectorized=True, batch=100, n_resamples=50)
    else:
        conf = None

    _, ax = plt.subplots()
    plot_population(time, population, ax=ax, conf_interval=conf,
                    label=labels, threshold=args.threshold)

    if args.save:
        plt.savefig(args.save)
    else:
        plt.show()
