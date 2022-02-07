#!/usr/bin/env python3
""" Run a set of QM calculations. """
import os
import itertools
from pathlib import Path
from argparse import ArgumentParser
import numpy as np
import file_utils
import turbomole_nstate as turbomole
import bagel

FILES = {
    "geom": "qm_geom",
    "state": "qm_state",
    "energy": "qm_energy",
    "gradient": "qm_grad",
    "oscill": "qm_oscill",
}

INTERFACES = {
        "tm_adc2": turbomole.ADC2,
        }


def main():
    """ Run a set of QM calculations using a specified interface. """
    parser = ArgumentParser(
        description="Run an electronic structure calculation.")
    parser.add_argument("min_nstate", type=int, help="Minimum number of states.")
    parser.add_argument("max_nstate", type=int, help="Maximum number of states.")
    parser.add_argument("interface", type=str, help="Interface to use.")
    parser.add_argument("work_dir", type=str, help="Working directory.")
    args = parser.parse_args()
    try:
        run(args)
    except:
        cleanup()
        raise


def run(args):
    """ Read input files and run the interface. """
    in_data = dict()
    in_data["geom"] = np.loadtxt(FILES["geom"])
    in_data["state"] = np.loadtxt(FILES["state"], dtype=int).flat[0]
    in_data["natom"] = len(in_data["geom"])
    in_data["min_nstate"] = args.min_nstate
    in_data["max_nstate"] = args.max_nstate
    qm_prog = INTERFACES[args.interface](in_data)
    # Go to work directory and run calculation.
    cwd = Path(os.getcwd())
    os.chdir(args.work_dir)
    qm_prog.update_input()
    qm_prog.run()
    os.chdir(cwd)
    # Write results.
    calc_write(args.max_nstate, qm_prog.results)


def calc_write(max_nstate, results):
    """ Write output files. """
    num_format = " {:20.14f}"
    xyz_format = 3 * " {:20.14f}"
    pad = np.zeros(max_nstate - len(results["energy"]))
    with open(FILES["energy"], "w") as out_f:
        for energy in itertools.chain(results["energy"], pad):
            out_f.write(num_format.format(energy))
        out_f.write("\n")
    with open(FILES["gradient"], "w") as out_f:
        for grad in results["gradient"]:
            out_f.write(xyz_format.format(*grad) + "\n")
    if "oscill" in results:
        with open(FILES["oscill"], "w") as out_f:
            for oscill in itertools.chain(results["oscill"], pad):
                out_f.write(num_format.format(oscill))
            out_f.write("\n")


def cleanup():
    for path in FILES.values():
        file_utils.remove(path)


if __name__ == "__main__":
    main()
