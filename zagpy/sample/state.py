""" Select points/states for starting FSSH dynamics. """
import os
import shutil
import numpy as np
from zagpy.units import eV


EN_FILE = 'qm_en.dat'
ADT_FILE = 'qm_adt.dat'


def add_subparser(subparsers):
    """ Register the 'state' subcommand under the 'sample' subparser. """
    parser = subparsers.add_parser(
        "state",
        formatter_class=__import__('argparse').ArgumentDefaultsHelpFormatter,
        description="Select points/states for starting FSSH dynamics.")
    parser.add_argument("--energy-range", type=float, nargs=2, default=None,
                        help="Energy window (in eV) for sampling.")
    parser.add_argument("--diabatic-state", type=int, default=None,
                        help="Diabatic state to select.")
    parser.add_argument("--adiabatic-state", type=int, default=None,
                        help="Adiabatic state to select.")
    parser.add_argument("--target-dir", type=str, default='start_trajs',
                        help="Directory to save selected points/states.")
    parser.add_argument("dirs", nargs="+", type=str,
                        help="Path to directories with single point calculations.")
    parser.set_defaults(func=run)


def run(args):
    """ Execute state selection with the given parsed arguments. """
    ex_en = []
    for cdir in args.dirs:
        if not os.path.isdir(cdir):
            print(f"Error: {cdir} is not a valid directory.")
            return
        dir_en = np.loadtxt(os.path.join(cdir, EN_FILE), usecols=1)
        ex_en.append(dir_en - dir_en[0])
    ex_en = np.array(ex_en) * eV

    cmask = np.ones_like(ex_en, dtype=bool)

    if args.energy_range is not None:
        cmask = np.logical_and(cmask, ex_en >= args.energy_range[0])
        cmask = np.logical_and(cmask, ex_en <= args.energy_range[1])

    weights = np.ones_like(ex_en)
    if args.diabatic_state is not None:
        adt = []
        for cdir in args.dirs:
            adt.append(np.loadtxt(os.path.join(cdir, ADT_FILE)))
        adt = np.array(adt)
        weights = weights * adt[:, args.diabatic_state - 1, :]**2

    if args.adiabatic_state is not None:
        amask = np.zeros_like(ex_en, dtype=bool)
        amask[:, args.adiabatic_state - 1] = True
        cmask = np.logical_and(cmask, amask)

    rng = np.random.default_rng()
    cmask = np.logical_and(cmask, rng.random(ex_en.shape) < weights)

    selected = []
    for i in range(ex_en.shape[0]):
        for j in range(ex_en.shape[1]):
            if cmask[i, j]:
                selected.append((i, j+1))

    try:
        os.mkdir(args.target_dir)
    except FileExistsError:
        print(f"Error: {args.target_dir} already exists.")
        return

    for i, j in selected:
        tdir = os.path.join(args.target_dir, f"{args.dirs[i]}_{j:02d}")
        shutil.copytree(args.dirs[i], tdir)
        with open(os.path.join(tdir, 'istate'), 'w') as wfile:
            wfile.write(f'{j}\n')
