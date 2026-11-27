""" Select points/states for starting FSSH dynamics. """
import json
import os
import shutil
import numpy as np
from zagpy.units import eV


OUT_FILE = 'qm_out.json'


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
    parser.add_argument("--weight-oscill", action="store_true", default=False,
                        help="Weight selection by oscillator strength.")
    parser.add_argument("--weight-norm", type=str, choices=["none", "max"],
                        default="max",
                        help="Normalization of weights: 'max' divides by the maximum weight.")
    parser.add_argument("--target-dir", type=str, default='start_trajs',
                        help="Directory to save selected points/states.")
    parser.add_argument("dirs", nargs="+", type=str,
                        help="Path to directories with single point calculations (each with "
                        f"the {OUT_FILE} response of eval_point).")
    parser.set_defaults(func=run)


def read_point(cdir):
    """ Read the response of eval_point for a single point.

    Returns the energies of the singlet states, the oscillator strengths from the first state
    (None if not in the response) and the diabatic transformation (None if not in the
    response), with element [d, a] the coefficient of diabatic state d in adiabatic state a."""
    with open(os.path.join(cdir, OUT_FILE), encoding="utf-8") as infile:
        data = json.load(infile)
    energies = None
    for block in data["states"]:
        if block["multiplicity"] == 1:
            energies = np.array(block["energies"])
    if energies is None:
        raise ValueError(f"{cdir}/{OUT_FILE} has no singlet energies.")
    props = data.get("properties", {})
    oscill = props.get("oscillator_strengths")
    adt = props.get("diabatic_transformation")
    return (energies,
            None if oscill is None else np.array(oscill),
            None if adt is None else np.array(adt))


def run(args):
    """ Execute state selection with the given parsed arguments. """
    ex_en = []
    oscill = []
    adt = []
    for cdir in args.dirs:
        if not os.path.isdir(cdir):
            print(f"Error: {cdir} is not a valid directory.")
            return
        en_data, osc_data, adt_data = read_point(cdir)
        ex_en.append(en_data - en_data[0])
        if args.weight_oscill:
            if osc_data is None:
                print(f"Error: {cdir}/{OUT_FILE} has no oscillator strengths "
                      "(property oscillator_strengths of the request).")
                return
            oscill.append(osc_data)
        if args.diabatic_state is not None:
            if adt_data is None:
                print(f"Error: {cdir}/{OUT_FILE} has no diabatic transformation "
                      "(property diabatic_transformation of the request).")
                return
            adt.append(adt_data)
    ex_en = np.array(ex_en) * eV

    cmask = np.ones_like(ex_en, dtype=bool)

    if args.energy_range is not None:
        cmask = np.logical_and(cmask, ex_en >= args.energy_range[0])
        cmask = np.logical_and(cmask, ex_en <= args.energy_range[1])

    weights = np.ones_like(ex_en)
    if args.weight_oscill:
        weights = weights * np.array(oscill)

    if args.diabatic_state is not None:
        # Weight of each adiabatic state: its squared overlap with the diabatic state.
        adt = np.array(adt)
        weights = weights * adt[:, args.diabatic_state - 1, :]**2

    if args.adiabatic_state is not None:
        amask = np.zeros_like(ex_en, dtype=bool)
        amask[:, args.adiabatic_state - 1] = True
        cmask = np.logical_and(cmask, amask)

    if args.weight_norm == "max":
        wmax = weights[cmask].max() if np.any(cmask) else 1.0
        if wmax > 0:
            weights = weights / wmax

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
