""" Interface for running Turbomole ADC(2) calculations. """
import subprocess
import re
import os
import numpy as np
import file_utils

STDOUT = "qm.log"
STDERR = "qm.err"

class Turbomole():
    """ Interface for turbomole calculations. """
    def __init__(self, data):
        self.data = data
        self.results = dict()

    def update_coord(self):
        """ Update coord file with self.data["geom"]. """
        file_utils.replace_cols_inplace("coord", self.data["geom"], r"\$coord")
        if self.data["qmmm"]:
            fn = file_utils.search_file("control", r"\$point_charges")[0]
            fn = fn.split("=")[1].split()[0]
            file_utils.replace_cols_inplace(fn, self.data["mm_geom"], r"\$point_charges")

    def update_state(self):
        """ Update control file to request gradient of self.data["state"]. """
        raise NotImplementedError("Need to call specific interface.")

    def update_input(self):
        """ Update all input files with values from self.data.  """
        self.update_coord()
        self.update_state()

    def run(self):
        """ Run the calculation, check success and read results. """
        raise NotImplementedError("Need to call specific interface.")

    def read(self):
        """ Read calculation results. """
        raise NotImplementedError("Need to call specific interface.")


class ricc2(Turbomole):
    """ Interface for turbomole ricc2 calculations. """

    def __init__(self, model, data):
        self.data = data
        self.model = model
        self.gs_model = model
        if model == "adc(2)":
            self.gs_model = "mp2"
        self.results = dict()

    def update_state(self):
        """ Update control file to request gradient of self.data["state"]. """
        if self.data["state"] == 1:
            state = r"(x)"
        else:
            state = r"(a {})".format(self.data["state"]-1)
        search_string = r"(geoopt +model={} +state=).*".format(re.escape(self.model))
        n_sub = file_utils.replace_inplace("control", search_string, r"\1" + state)
        if n_sub == 0:
            repl = r"$ricc2\n  geoopt model={} state={}".format(self.model, state)
            n_sub = file_utils.replace_inplace("control", r"\$ricc2", repl)

    def run(self):
        """ Run the calculation, check success and read results. """
        with open(STDOUT, "w") as out, open(STDERR, "w") as err:
            subprocess.run("dscf", stdout=out, stderr=err)
            actual_check()
            subprocess.run("ricc2", stdout=out, stderr=err)
            actual_check()

    def read(self):
        self.results["energy"] = ricc2_energy(STDOUT, self.gs_model)
        self.results["gradient"] = ricc2_gradient()[self.data["state"]]
        try:
            self.results["oscill"] = ricc2_oscill(STDOUT)
        except:
            pass


class mp2(ricc2):
    """ Interface for MP2 ground state calculations. """

    def read(self):
        self.results["energy"] = ricc2_gs_energy(self.gs_model)
        self.results["gradient"] = ricc2_gradient()[1]


class egrad(Turbomole):
    def __init__(self, data):
        self.data = data
        try:
           soes = file_utils.search_file("control",  r"\$soes", after=1)
           soes = file_utils.split_columns(soes, col=1, convert=int)[0]
           self.data["n_ex_state"] = soes
        except:
            self.data["n_ex_state"] = 0
        try:
            _ = file_utils.search_file("control", r"\$rij")
            self.data["ri"] = True
        except:
            self.data["ri"] = False
        self.results = dict()


    def update_state(self):
        ex_state = self.data["state"] - 1
        if self.data["n_ex_state"] < ex_state:
            raise ValueError("Not enough states selected in QM calculation.")
        if ex_state == 0:
            return
        n_sub = file_utils.replace_inplace("control",
                r"\$exopt.*", 
                r"$exopt {}".format(ex_state))
        if n_sub == 0:
            file_utils.replace_inplace("control",
                r"\$end", 
                r"$exopt {}\n$end".format(ex_state))


    def run(self):
        # Remove existing gradient file to avoid the file becoming huge.
        try:
            os.remove("gradient")
        except:
            pass
        # Run new calculation.
        with open(STDOUT, "w") as out, open(STDERR, "w") as err:
            if self.data["ri"]:
                subprocess.run("ridft", stdout=out, stderr=err)
            else:
                subprocess.run("dscf", stdout=out, stderr=err)
            actual_check()
            if self.data["state"] == 1:
                if self.data["ri"]:
                    subprocess.run("rdgrad", stdout=out, stderr=err)
                else:
                    subprocess.run("grad", stdout=out, stderr=err)
                if self.data["n_ex_state"] > 0:
                    subprocess.run("escf", stdout=out, stderr=err)
            else:
                subprocess.run("egrad", stdout=out, stderr=err)
            actual_check()


    def read(self):
        self.results["energy"] = tddft_energy(STDOUT)
        self.results["gradient"] = get_grad_from_gradient(STDOUT, self.data["natom"])
        if self.data["n_ex_state"] > 0:
            self.results["oscill"] = tddft_oscill(STDOUT)


# TDDFT Calculation keywords
DSCF_EN = re.escape(r"|  total energy      =")
ESCF_EN = r"Total energy:"
GRAD_GRAD = r"SCF ENERGY GRADIENT with respect to NUCLEAR COORDINATES"
EGRAD_GRAD = r"Excited state no.*chosen for optimization"
RDGRAD_GRAD = r"RDGRAD - INFORMATION"

def tddft_energy(fname):
    try:
        energy = file_utils.search_file(fname, ESCF_EN)
        col = 2
    except ValueError:
        energy = file_utils.search_file(fname, DSCF_EN)
        col = 4
    file_utils.split_columns(energy, col=col, convert=np.float64)
    return np.array(energy)


def tddft_gradient(fname, target, ri):
    if target == 1:
        if ri:
            cfile = file_utils.go_to_keyword(fname, RDGRAD_GRAD)[0]
        else:
            cfile = file_utils.go_to_keyword(fname, GRAD_GRAD)[0]
    else:
        cfile = file_utils.open_if_needed(fname)
        while True:
            cfile, cstate = file_utils.go_to_keyword(cfile, EGRAD_GRAD)
            cstate = int(cstate.split()[3]) + 1
            if cstate == target:
                break
    grad = get_grad_from_stdout(cfile)
    return grad


def tddft_oscill(fname):
    oscill = file_utils.search_file(fname, "mixed representation:")
    file_utils.split_columns(oscill, col=2, convert=np.float64)
    return np.array(oscill)


def ricc2_energy(fname, model):
    search_string = "Final " + re.escape(model.upper()) + " energy"
    gs_energy = file_utils.search_file(fname, search_string)
    file_utils.split_columns(gs_energy, col=5, convert=np.float64)
    ex_energy = file_utils.search_file(fname, "Energy:")
    ex_energy = file_utils.split_columns(ex_energy, 1, convert=np.float64)
    energy = np.repeat(gs_energy, len(ex_energy) + 1)
    energy[1:] = energy[1:] + ex_energy
    return energy


def ricc2_gs_energy(model):
    energy = file_utils.search_file(STDOUT, "Total Energy  ")
    energy = np.float(energy[0].split()[3])
    return np.array([energy])


def ricc2_oscill(fname):
    """ Read oscillator strengths from STDOUT file. """
    oscill = file_utils.search_file(fname,
                                    r"oscillator strength \(length gauge\)")
    oscill = file_utils.split_columns(oscill, col=5)
    return np.array(oscill, dtype=np.float64)


def ricc2_gradient():
    grads = dict()
    # Try to get ground state gradient.
    try:
        cfile = file_utils.go_to_keyword(STDOUT, "GROUND STATE FIRST-ORDER PROPERTIES")[0]
        grads[1] = get_grad_from_stdout(cfile)
    except ValueError:
        pass
    # Try to get excited state gradients.
    try:
        cfile = file_utils.go_to_keyword(STDOUT, "EXCITED STATE PROPERTIES")[0]
    except ValueError:
        return grads
    while True:
        try:
            line = file_utils.search_file(cfile, 
                    "Excited state reached by transition:", 
                    max_res=1, 
                    close=False, 
                    after=3)
            cstate = int(line[0].split()[4]) + 1
        except ValueError:
            cfile.close()
            break
        try:
            file_utils.search_file(cfile, 
                    "cartesian gradient of the energy",
                    max_res=1,
                    stop_at=r"\+={73}\+",
                    close=False)
            grads[cstate] = get_grad_from_stdout(cfile)
        except ValueError:
            pass
    return grads


def get_grad_from_gradient(natom):
    grad = file_utils.search_file("gradient", r"cycle", after=2*natom)[-natom:]
    grad = file_utils.split_columns(grad, col=[0, 1, 2],  convert=file_utils.fortran_double)
    return grad


def get_grad_from_stdout(cfile):
    grad = file_utils.search_file(cfile, r"^  ATOM", after=3, stop_at=r"resulting FORCE", close=False)
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
