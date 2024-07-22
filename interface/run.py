#!/usr/bin/env python3
""" Run a set of QM calculations. """
import os
from pathlib import Path
from argparse import ArgumentParser
import file_utils
import numpy as np
from turbomole import Turbomole
from orca import Orca

FILES = {
    "geom": "qm_geom",
    "state": "qm_state",
    "energy": "qm_energy",
    "gradient": "qm_grad",
    "oscill": "qm_oscill",
    "mm_geom": "mm_geom",
}

INTERFACES = {
        "turbomole": Turbomole,
        "orca" : Orca
        }


def main():
    """ Run a set of QM calculations using a specified interface. """
    parser = ArgumentParser(
        description="Run an electronic structure calculation.")
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
    in_data["iroot"] = np.loadtxt(FILES["state"], dtype=int)[0]
    in_data["nstate"] = np.loadtxt(FILES["state"], dtype=int)[1]
    try:
        in_data["mm_geom"] = np.loadtxt(FILES["mm_geom"])
        in_data["qmmm"] = True
    except OSError:
        in_data["qmmm"] = False
    # Go to work directory and run calculation.
    cwd = Path(os.getcwd())
    os.chdir(args.work_dir)
    interface = INTERFACES[args.interface]
    qm_prog = interface()
    request = {}
    request["gradient"] =  True
    if in_data["nstate"] > 1:
        request["oscillator_strength"] = True
    else:
        request["oscillator_strength"] = False
    qm_prog.update_input(in_data, request)
    qm_prog.run()
    qm_prog.read()
    os.chdir(cwd)
    # Write results.
    calc_write(qm_prog.data)


def calc_write(results):
    """ Write output files. """
    num_format = " {:20.14f}"
    xyz_format = 3 * " {:20.14f}"
    with open(FILES["energy"], "w") as out_f:
        for energy in results["energy"]:
            out_f.write(num_format.format(energy))
        out_f.write("\n")
    with open(FILES["gradient"], "w") as out_f:
        for grad in results["gradient"]:
            out_f.write(xyz_format.format(*grad) + "\n")
    if "oscillator_strength" in results:
        with open(FILES["oscill"], "w") as out_f:
            for oscill in results["oscillator_strength"]:
                out_f.write(num_format.format(oscill))
            out_f.write("\n")


def cleanup():
    for path in FILES.values():
        file_utils.remove(path)


if __name__ == "__main__":
    main()
