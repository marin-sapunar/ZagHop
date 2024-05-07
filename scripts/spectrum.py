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
        help="Name of the files that contain oscillator strengths and excitation energies.")
    parser.add_argument(
        "-b",
        "--bandwidth",
        nargs='?',
        type=float,
        metavar=("bandwidth"),
        help="Give the value of the bandwith of Gaussian functions.")
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
               #intensity += 1.22 * 10**(-44)* 1/bandwidth * 1/100 * float(all_os[i][j]) * float(all_excited[i][j]) * np.exp(-1/2 *((en - float(all_excited[i][j]))/bandwidth) **2)
                intensity +=  1/bandwidth * float(all_os[i][j]) * float(all_excited[i][j]) * np.exp(-1/2 *((en - float(all_excited[i][j]))/bandwidth) **2)
        return float(intensity)

    f = np.vectorize(nea_spect)

    en = np.arange(3, 8, 0.1)

    plt.plot(en, f(en))

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
