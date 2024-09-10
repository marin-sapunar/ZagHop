#!/usr/bin/env python3
import argparse
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import statistics
from recombinator.iid_bootstrap import \
    iid_balanced_bootstrap, \
    iid_bootstrap, \
    iid_bootstrap_vectorized, \
    iid_bootstrap_via_choice, \
    iid_bootstrap_via_loop, \
    iid_bootstrap_with_antithetic_resampling
import itertools

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
        help="Insert the name of a file that contains oscillator strengths and excitation energies.")
    parser.add_argument(
        "-b",
        "--bandwidth",
        default=0.15,
        type=float,
        metavar=("bandwidth"),
        help="Give the value of the bandwith for Gaussian functions.")
    parser.add_argument(
        "-start",
        "--start-value",
        type=float,
        metavar=("start_value"),
        default=1.0,
        help="Set the start value of energy on the x-axis to be displayed.")
    parser.add_argument(
        "-end",
        "--end-value",
        type=float,
        metavar=("end_value"),
        default=10.0,
        help="Set the end value of energy on the x-axis to be displayed.")
    parser.add_argument(
        "-dots",
        "--dots",
        type=float,
        metavar=("dots"),
        default=100,
        help="Set the number of dots to create the spectrum.")
    parser.add_argument(
        "-std",
        "--std",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Calculate the standard deviation of the spectrum using Bootstrap method and show it on a spectrum as a shaded area around the mean value.")
    parser.add_argument(
        "-r",
        "--resamples",
        type=int,
        metavar=("resamples"),
        default=100,
        help="If calculating the standard deviation, insert the number of resamples for the Bootstrap method.")

    args = parser.parse_args()

    input_file = args.input_file

    bw = args.bandwidth

    start = args.start_value

    end = args.end_value

    dots = args.dots

    R = args.resamples

    data = np.loadtxt(input_file)

    energies = data[:, 0]
    os_strengths = data[:, 1]

    exc_en = energies.tolist()

    osc_str = os_strengths.tolist()

    en = np.linspace(start, end, dots)

    def nea_spect_initial(en):
        kernel = stats.gaussian_kde(exc_en, bw_method=bw, weights=osc_str)
        return np.array(kernel(en))

    def nea_spect_bootstrap(en):
        results = []
        for i in range(R):
            kernel = stats.gaussian_kde(res_en[i], bw_method=bw, weights=res_os[i])
            results.append(kernel(en))
        return np.array(results)

    def plot_initial():
        results = nea_spect_initial(en)
        plt.plot(en, results)
        plt.xlabel("Energy in eV")
        plt.ylabel("Intensity")
        plt.title("UV spectrum calculated using NEA")
        plt.show()

    def plot_bootstrap():
        results = nea_spect_bootstrap(en)

        mean_kde = np.mean(results, axis=0)
        std_kde = np.std(results, axis=0)

        plt.plot(en, mean_kde, label='Mean KDE')

        plt.fill_between(en, mean_kde - std_kde, mean_kde + std_kde, color='gray', alpha=0.3, label='±1 Std Dev')

        plt.xlabel("Energy in eV")
        plt.ylabel("Intensity")
        plt.title("Bootstrap Resampled Spectra with Standard Deviation")
        plt.legend()
        plt.show()

    if args.std:
        joined = [list(pair) for pair in zip(exc_en, osc_str)]
        merged = np.array(joined)
        resampled = iid_bootstrap(merged, replications=R, replace=True)
        res_en = []
        res_os = []
        for inner_list in resampled:
            first, second = zip(*inner_list)
            res_en.append(list(first))
            res_os.append(list(second))
        plot_bootstrap()
    else:
        plot_initial()

if __name__ == "__main__":
    main()

