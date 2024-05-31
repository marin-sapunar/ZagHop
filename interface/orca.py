""" Interface for running Orca calculations. """
import subprocess
import sys
import numpy as np
import file_utils
import interface


class Orca(interface.QMInterface):
    """Interface for QM calculations using Orca.

    This class provides methods for setting up and running calculations, as well as
    reading and retrieving results from the Orca output files.

    See the `QMInterface` class for more details on the methods and attributes.
    """

    def __init__(self, template="qm.inp", log_file="qm.log", err_file="qm.err", **kwargs):
        """ Initialize the Orca interface object. """
        super().__init__(template=template, log_file=log_file, err_file=err_file, **kwargs)
        xyzfile = file_utils.search_file(self.template,
                               r"(?i)xyzfile\s*\d*\s*\d*\s*(\S*)",
                               matching_only=True)
        if self.template.endswith(".inp"):
            self.opts["base_name"] = self.template.replace(".inp", "")
        else:
            self.opts["base_name"] = self.template
        self.opts["xyz_file"] = xyzfile[0][0]

    def run(self):
        """ Run the calculation and check success. """
        with open(self.log_file, "w") as outf:
            with open(self.err_file, "w") as errf:
                try:
                    subprocess.run(["orca",
                                    self.template],
                                    stdout=outf,
                                    stderr=errf,
                                    check=True)
                except subprocess.CalledProcessError:
                    print("Error in Orca interface. QM calculation failed.")
                    print("  Check " + self.err_file + "/" + self.log_file + ".")
                    sys.exit(1)
                try:
                    subprocess.run(["orca_2mkl", self.opts["base_name"], "-molden"],
                                    stdout=outf,
                                    stderr=errf,
                                    check=True)
                except subprocess.CalledProcessError:
                    print("Error in Orca interface. orca_2mkl call failed.")
                    print("  Check " + self.err_file + "/" + self.log_file + ".")
                    sys.exit(1)
        self.calc_done = True

    def read_energy(self):
        """Read energies from the log file into the data dictionary."""
        gs_en = file_utils.search_file(self.log_file,
                                       r"Total Energy\s*\:\s*(\S*)\s*Eh",
                                       matching_only=True)[0][0]
        gs_en = np.array(gs_en, ndmin=1, dtype=float)
        if self.data["nstate"] == 1:
            self.data["energy"] = gs_en
            return
        ex_en = file_utils.search_file(self.log_file,
                                       r"STATE\s*\d*\:\s*E=\s*(\S*)\s*au",
                                       matching_only=True)
        ex_en = np.array([[0.0]] + ex_en, dtype=float).flatten()
        self.data["energy"] = ex_en + gs_en

    def read_gradient(self):
        """ Read gradient from the log file into the data dictionary. """
        grad = file_utils.search_file(self.log_file,
                                      r"CARTESIAN GRADIENT",
                                      after=2+self.data["natom"])[2:]
        grad = file_utils.split_columns(grad, col=[3, 4, 5])
        self.data["gradient"] = np.array(grad, dtype=float)

    def read_oscill(self):
        """ Read oscillator strengths from the log file into the data dictionary. """
        n_ex_state = self.data["nstate"] - 1
        oscill = file_utils.search_file(self.log_file,
                                        r"VIA TRANSITION ELECTRIC DIPOLE MOMENTS",
                                        after=4+n_ex_state)[4:]
        oscill = file_utils.split_columns(oscill, col=3)
        self.data["oscillator_strength"] = np.array(oscill, dtype=float)

    def update_geom(self):
        """ Update geometry in QM input files. """
        self.data["natom"] = len(self.data["geom"])
        file_utils.replace_cols_inplace(self.opts["xyz_file"],
                                        self.data["geom"] * 0.52917721067, #@todo units
                                        r".", # Match first line
                                        skip=1, # Skip comment line
                                        cols = [1, 2, 3])

    def update_nstate(self):
        """ Update number of states requested from the QM code. """
        if self.data["nstate"] > 1:
            file_utils.replace_inplace(self.template,
                                       r"(?i)NROOTS\s*\d*",
                                       "NROOTS " + str(self.data["nstate"] - 1))

    def update_iroot(self):
        """ Request gradient for specified state. """
        if not file_utils.check_in_file(self.template, r"(?i)!.*engrad"):
            file_utils.replace_inplace(self.template,
                                       r"(?i)^!(.*)engrad",
                                       "!\1engrad ")
        file_utils.replace_inplace(self.template,
                                   r"(?i)IROOT\s*\d*",
                                   "IROOT " + str(self.data["iroot"] - 1))
