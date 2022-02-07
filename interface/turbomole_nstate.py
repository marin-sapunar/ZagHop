""" Interface for running Turbomole ADC(2) calculations. """
import subprocess
import re
import numpy as np
import file_utils

STDOUT = "qm.log"
STDERR = "qm.err"


class ADC2():
    """ Interface. """

    def __init__(self, data):
        self.data = data
        self.results = dict()

    def update_input(self):
        """ Update input files with values from self.data. """
        # Update coord file.
        file_utils.replace_cols_inplace("coord", self.data["geom"], r"\$coord")
        # Update requested state in control file.
        if self.data["state"] == 1:
            state_string = "(x)"
        else:
            state_string = "(a {})".format(self.data["state"]-1)
        file_utils.replace_inplace("control",
                                   r"(geoopt +model=adc\(2\) +state=).*",
                                   r"\1" + state_string)
        n_state = min(self.data["max_nstate"]-1, self.data["state"])
        n_state = max(self.data["min_nstate"]-1, n_state)
        file_utils.replace_inplace("control",
                                   r"(\s*irrep.*)nexc\s*=\s*\d+(.*)",
                                   r"\1 nexc="+str(n_state)+r" \2")

    def run(self):
        """ Run the calculation, check success and read results. """
        with open(STDOUT, "w") as out, open(STDERR, "w") as err:
            subprocess.run("dscf", stdout=out, stderr=err)
            actual_check()
            subprocess.run("ricc2", stdout=out, stderr=err)
            actual_check()
        self.read()

    def read(self):
        self.results["energy"] = read_energy()
        self.results["gradient"] = read_gradient()[self.data["state"]]
        self.results["oscill"] = read_oscill()


def read_energy():
    gs_energy = file_utils.search_file(STDOUT, "Final MP2 energy")
    file_utils.split_columns(gs_energy, col=5, convert=np.float64)
    ex_energy = file_utils.search_file(STDOUT, "Energy:")
    ex_energy = file_utils.split_columns(ex_energy, 1, convert=np.float64)
    energy = np.repeat(gs_energy, len(ex_energy) + 1)
    energy[1:] = energy[1:] + ex_energy
    return energy


def read_oscill():
    """ Read oscillator strengths from STDOUT file. """
    oscill = file_utils.search_file(STDOUT,
                                    r"oscillator strength \(length gauge\)")
    oscill = file_utils.split_columns(oscill, col=5)
    return np.array(oscill, dtype=np.float64)


def read_gradient():
    grads = dict()
    # Try to get ground state gradient.
    try:
        cfile = file_utils.go_to_keyword(STDOUT, "GROUND STATE FIRST-ORDER PROPERTIES")[0]
        grads[1] = get_grad_from_stdout(cfile)
    except:
        pass
    # Try to get excited state gradients.
    try:
        cfile = file_utils.go_to_keyword(STDOUT, "EXCITED STATE PROPERTIES")[0]
    except:
        return grads
    while True:
        try:
            line = file_utils.search_file(cfile, 
                    "Excited state reached by transition:", 
                    max_res=1, 
                    close=False, 
                    after=3)
            cstate = int(line[0].split()[4]) + 1
        except:
            cfile.close()
            break
        try:
            file_utils.search_file(cfile, 
                    "cartesian gradient of the energy",
                    max_res=1,
                    stop_at=r"\+={73}\+",
                    close=False)
            grads[cstate] = get_grad_from_stdout(cfile)
        except:
            pass
    return grads


def get_grad_from_stdout(cfile):
    grad = file_utils.search_file(cfile, r"ATOM", after=3, stop_at=r"resulting FORCE", close=False)
    grad = [line[5:] for line in grad]
    grad = [' '.join(grad[0::3]), ' '.join(grad[1::3]), ' '.join(grad[2::3])]
    grad = [line.split() for line in grad]
    grad = [list(map(file_utils.fortran_double, vals)) for vals in grad]
    return np.array(grad).T


def actual_check():
    """ Check that the Turbomole calculation finished without error. """
    check = subprocess.run("actual", stdout=subprocess.PIPE)
    if check.stdout.startswith(b'fine, there is no data group "$actual step"'):
        return
    raise RuntimeError("Turbomole calculation failed.")
