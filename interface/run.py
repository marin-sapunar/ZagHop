#!/usr/bin/env python3
""" Run a set of QM calculations. """
import os
import sys
from pathlib import Path
from argparse import ArgumentParser
import json
import numpy as np
import file_utils
from turbomole import Turbomole
from orca import Orca


INTERFACES = {
        "turbomole": Turbomole,
        "orca" : Orca
        }


def cli():
    parser = ArgumentParser(
        description="Run an electronic structure calculation.")
    parser.add_argument(
        "-init",
        "--init-inputs",
        action="store_true",
        help="Initialize input files for the interface.")
    parser.add_argument("interface", type=str, help="Interface to use.")
    parser.add_argument("work_dir", type=str, help="Working directory.")
    return parser


def run(args):
    """ Read input files and run the interface. """
    with open("qm.json", "r") as infile:
        in_data = json.load(infile)
    # @todo This is a messy patch to update new qm.json format
    #       with qm.yaml format from previous version. Fix.
    states = in_data.pop("states")[0]
    sys_index = states.pop("system")
    system = in_data.pop("system")[sys_index]
    system["states"] = {"nstate" : states.pop("nstate")}
    request = states
    if "geom" not in system:
        print("Error, geom not found in qm.json file.")
        sys.exit(1)
    system["geom"] = np.array(system["geom"], dtype=float)

    # Go to work directory and run calculation.
    cwd = Path(os.getcwd())
    os.chdir(args.work_dir)
    qm_prog = INTERFACES[args.interface]()
    if args.init_inputs:
        options = in_data.pop("options")
        qm_prog.generate_inputs(system, options)
    qm_prog.check_template()
    qm_prog.update(system, request)
    qm_prog.run()
    qm_prog.read()
    os.chdir(cwd)
    return qm_prog.results


def calc_write(results):
    """ Write output files. """
    for key in results:
        # Convert numpy arrays to lists for yaml.dump
        if isinstance(results[key], np.ndarray):
            results[key] = results[key].tolist()
    with open("qm_out.yaml", "w") as outfile:
        yaml.dump(results, outfile)


def main():
    """ Run a set of QM calculations using a specified interface. """
    parser = cli()
    args = parser.parse_args()
    try:
        results = run(args)
        calc_write(results)
    except:
        file_utils.remove("qm_out.yaml")
        raise


if __name__ == "__main__":
    main()
