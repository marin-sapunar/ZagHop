""" Interface for running Turbomole ADC(2) calculations. """
import subprocess
import re
import numpy as np
import file_utils

STDOUT = "qm.log"
STDERR = "qm.err"


class CASPT2():
    """ Interface. """

    def __init__(self, data):
        self.data = dict()
        self.data["input file"] = "template.json"
        self.data.update(data)
        self.results = dict()

    def update_input(self):
        """ Update input files with values from self.data. """
        update_coord(self.data["input file"], self.data["geom"])
        update_state(self.data["input file"], self.data["state"])

    def run(self):
        """ Run the calculation, check success and read results. """
        with open(STDOUT, "w") as out, open(STDERR, "w") as err:
            subprocess.run(['BAGEL', self.data["input file"]], stdout=out, stderr=err)
           #check()
        self.read()

    def read(self):
        self.results["energy"] = read_energy()
        self.results["gradient"] = read_gradient(self.data["state"])
        return self.results


def update_coord(infile, geom):
    split_re = re.compile(r"([\[\],])")
    file_utils.replace_cols_inplace(infile, geom, r"geometry", [4, 6, 8], split=split_re.split)


def update_state(infile, state):
    state_string = "{},".format(state)
    search_re = r'(\"target\"\ +\:\ +).*,'
    file_utils.replace_inplace(infile, search_re, r"\1 "+state_string)


def read_energy():
    return np.loadtxt("ENERGY.out")


def read_gradient(state):
    gfile = "FORCE_" + str(state) + ".out"
    grad = np.loadtxt(gfile, skiprows=1, usecols=[1, 2, 3])
    return grad
#
#
#def read_oscill():
#    return oscill
#
#
#def read_gradient():
#    """ Read gradient from STDOUT file. """
#    gradient = file_utils.search_file(STDOUT, 'o Atom', after=3, col=1)
#    natom = len(gradient) / 3
#    gradient = np.array(gradient, dtype=np.float64).reshape(natom, 3)
#    return gradient
#
#
#def read_energy():
#    """ Read energies from STDOUT file. """
#    energy = file_utils.search_file(STDOUT, "MS-CASPT2 energy")
#    energy = file_utils.split_columns(energy, col=6)
#    energy = np.array(energy, dtype=np.float64)
#    return energy
