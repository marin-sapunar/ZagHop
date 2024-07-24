#!/usr/bin/env python3
""" Run a set of QM calculations. """
import os
import yaml
from pathlib import Path
from argparse import ArgumentParser
import file_utils
import numpy as np
from turbomole import Turbomole
from orca import Orca


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
    with open("qm_sys.yaml", "r") as infile:
        in_data = yaml.safe_load(infile)
    if "geom" not in in_data:
        print("Error, geom not found in qm_sys.yaml file.")
        sys.exit(1)
    in_data["geom"] = np.array(in_data["geom"], dtype=float)
    request = in_data.pop("request")

    # Go to work directory and run calculation.
    cwd = Path(os.getcwd())
    os.chdir(args.work_dir)
    interface = INTERFACES[args.interface]
    qm_prog = interface()
    qm_prog.update_input(in_data, request)
    qm_prog.run()
    qm_prog.read()
    os.chdir(cwd)
    # Write results.
    calc_write(qm_prog.results)


def calc_write(results):
    """ Write output files. """
    for key in results:
        # Convert numpy arrays to lists for yaml.dump
        if isinstance(results[key], np.ndarray):
            results[key] = results[key].tolist()
    with open("qm_out.yaml", "w") as outfile:
        yaml.dump(results, outfile)


def cleanup():
    file_utils.remove("qm_out.yaml")


if __name__ == "__main__":
    main()
