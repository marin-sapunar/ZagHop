""" Interface for running CFOUR calculations. """
import subprocess
import re
import numpy as np
import file_utils

STDOUT = "qm.log"
STDERR = "qm.err"


class CC():
    """ Interface. """

    def __init__(self, data):
        self.data = data
        self.results = dict()

    def update_input(self):
        """ Update input files with values from self.data. """
        file_utils.replace_cols_inplace("ZMAT", self.data["geom"], r"", cols=[1,2,3])

    def run(self):
        """ Run the calculation, check success and read results. """
        with open(STDOUT, "w") as out, open(STDERR, "w") as err:
            subprocess.run("xcfour", stdout=out, stderr=err)
        self.read()

    def read(self):
        self.results["energy"] = read_energy()
        self.results["gradient"] = read_gradient(self.data["natom"])


def read_energy():
    key = 'Total .*CC.* energy'
    energy = file_utils.search_file(STDOUT, key)
    file_utils.split_columns(energy, col=3, convert=np.float64)
    return energy


def read_gradient():
    key = r"reordered gradient in QCOM coords for ZMAT order"
    grad = file_utils.search_file(STDOUT, key, after=natom)
    grad = file_utils.split_columns(grad, convert=np.float64)
    return np.array(grad)


