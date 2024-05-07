#!/usr/bin/env python3
import numpy as np
import argparse
import glob
import re
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser( 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=
        "Calculate the UV apsorption spectrum using the Nuclear Ensamble Approach.")
    parser.add_argument(
        "-in",
        "--input-file",
        type=str,
        metavar=("input_file"),
        default="qm.out",
        help="Name of the files that contain oscillator strengths and excitation energies.")
    parser.add_argument(
        "-b",
        "--bandwidth",
        default=0.1,
        type=float,
        metavar=("bandwidth"),
        help="Give the value of the bandwith of Gaussian functions.")
    parser.add_argument(
        "-eu",
        "--energy-unit",
        choices=["eV", "nm"],
        type=str,
        metavar=("energy_unit"),
        default="eV",
        help="Choose the units to express the energy on the x-axis of the spectrum. Available options are eV and nm.")
    parser.add_argument(
        "-n",
        "--normalization",
        action='store_true',
        help="Create the spectrum with intensities relative to the maximum peak which is set to have the value of 1.")
    parser.add_argument(
        "-sv",
        "--start-value",
        type=float,
        metavar=("start_value"),
        default=0,
        help="Set the start value of energy on the x-axis.")
    parser.add_argument(
        "-ev",
        "--end-value",
        type=float,
        metavar=("end_value"),
        default=10,
        help="Set the end value of energy on the x-axis.")
    parser.add_argument(
        "-s",
        "--step",
        type=float,
        metavar=("step"),
        default=0.01,
        help="Set the value of the step (distance between two points) of the energy axis.")

    args = parser.parse_args()

    input_file = args.input_file
    bandwidth = args.bandwidth

    all_os = []
    all_excited = []

    files_out = sorted(glob.glob('**/' + input_file, recursive=True))

    for file in files_out:
        os = read_os_from_ricc2out(file)
        all_os.append(os)
        
    for file in files_out:
        excited = read_en_from_ricc2out(file)
        all_excited.append(excited)
        
    def nea_spect(en):
        intensity = 0
        for i in range(len(all_os)):
            for j in range(len(all_os[0])):
                if args.energy_unit == "nm":
                    intensity +=  1/bandwidth * float(all_os[i][j]) * 1239.8/(float(all_excited[i][j])) * np.exp(-1/2 *((1239.8/en - 1239.8/(float(all_excited[i][j])))/bandwidth) **2)
                if args.emergy_unit == "eV":
                    intensity +=  1/bandwidth * float(all_os[i][j]) * float(all_excited[i][j]) * np.exp(-1/2 *((en - float(all_excited[i][j]))/bandwidth) **2)
        return float(intensity)

    f = np.vectorize(nea_spect)

    en = np.arange(args.start_value, args.end_value, args.step)

    plt.plot(en, f(en))

    if args.normalization:
        g= f(en)/f(en).max()



    #en = np.arange(3, 8, 0.1)

    plt.show()


def read_en_from_ricc2out(fname):
    with open(fname, "r") as ifile:
        list_lines = []
        lines = ifile.readlines()
    for line in lines:
        if re.search(r"frequency :", line):
            line_split = line.split()[5]
            list_lines.append(line_split)
    return list_lines

def read_os_from_ricc2out(fname):
    with open(fname, "r") as ifile:
        lines = ifile.readlines()
        list_os = []
    for line in lines:
        if re.search(r"oscillator strength", line):
            line_split = line.split()[5]
            list_os.append(line_split)
    return list_os



if __name__ == "__main__":
    main()
