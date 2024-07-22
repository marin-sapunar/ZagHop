""" Interface for running Turbomole RICC2 or TDDFT calculations. """
import subprocess
import re
import os
import sys
import numpy as np
import file_utils
from interface import QMInterface


class Turbomole(QMInterface):
    """ Interface for turbomole calculations. """
    def __init__(self, template="control", log_file="qm.log", err_file="qm.err", **kwargs):
        super().__init__(template=template, log_file=log_file, err_file=err_file, **kwargs)
        check, lines = get_group("ricc2")
        if check != 0:
            self.opts["ricc2"] = False
        else:
            self.opts["ricc2"] = True
            for line in lines:
                if line.strip() in ["adc(2)", "mp2"]:
                    self.opts["ricc2_model"] = line.strip()
                    break
            if "ricc2_model" not in self.opts.keys():
                print("Error in Turbomole interface.")
                print("Unrecognized wave function model in $ricc2 section.")
                sys.exit(1)
            # MP2 is the ground state for ADC(2) calculations
            if self.opts["ricc2_model"] == "adc(2)":
                self.opts["ricc2_gs_model"] = "mp2"
            else:
                self.opts["ricc2_gs_model"] = self.opts["ricc2_model"]
        check, lines = get_group("soes")
        if check != 0:
            self.opts["egrad"] = False
        else:
            self.opts["egrad"] = True
            # @todo egrad options.
        if file_utils.check_in_file(self.template, r"\$rij"):
            self.opts["ri"] = True
        else:
            self.opts["ri"] = False
        self.opts["qmmm"] = False

    def run(self):
        """ Run Turbomole calculation. """
        # Remove existing gradient file to avoid the file becoming huge.
        try:
            os.remove("gradient")
        except FileNotFoundError:
            pass
        # All subprocess calls use check=False since Turbomole programs don't
        # return a non-zero exit code on error anyways. Instead, to check for errors
        # `actual_check` calls "actual" and parses the output.
        with open(self.log_file, "w") as out, open(self.err_file, "w") as err:
            # Run the SCF calculation.
            if self.opts["ri"]:
                subprocess.run("ridft", stdout=out, stderr=err, check=False)
            else:
                subprocess.run("dscf", stdout=out, stderr=err, check=False)
            actual_check()
            if file_utils.check_in_file(self.err_file, r"ended abnormally"):
                print("Turbomole calculation failed. Check output in " + os.getcwd() + ".")
                sys.exit(1)
            # Run the excited state + gradient calculation(s).
            if self.opts["ricc2"]:
                subprocess.run("ricc2", stdout=out, stderr=err, check=False)
            else:
                if self.data["iroot"] == 1:
                    # Ground state gradient calculation.
                    if self.opts["ri"]:
                        subprocess.run("rdgrad", stdout=out, stderr=err, check=False)
                    else:
                        subprocess.run("grad", stdout=out, stderr=err, check=False)
                    actual_check()
                    # Separate excited state calculation.
                    if self.data["nstate"] > 1:
                        subprocess.run("escf", stdout=out, stderr=err, check=False)
                else:
                    subprocess.run("egrad", stdout=out, stderr=err, check=False)
            actual_check()

    def read_energy(self):
        """ Read energy from log file. """
        if self.opts["ricc2"]:
            self.data["energy"] = get_ricc2_energy(self.log_file, self.opts["ricc2_gs_model"])
        else:
            self.data["energy"] = get_dft_energy(self.log_file)

    def read_gradient(self):
        """ Read gradient from log file. """
        if self.opts["ricc2"]:
            self.data["gradient"] = get_ricc2_gradient(self.log_file, self.data["iroot"])
        else:
            self.data["gradient"] = get_dft_gradient(self.log_file, self.data["iroot"],
                                                     self.opts["ri"])

    def read_oscill(self):
        """ Read oscillator strengths from log file. """
        if self.opts["ricc2"]:
            self.data["oscillator_strength"] = get_ricc2_oscill(self.log_file)
        else:
            self.data["oscillator_strength"] = get_tddft_oscill(self.log_file)

    def update_geom(self):
        """ Update coord file with self.data["geom"]. """
        file_utils.replace_cols_inplace("coord", self.data["geom"], r"\$coord")
        if self.opts["qmmm"]:
            fn = file_utils.search_file("control", r"\$point_charges")[0]
            fn = fn.split("=")[1].split()[0]
            file_utils.replace_cols_inplace(fn, self.data["mm_geom"], r"\$point_charges")

    def update_nstate(self):
        """ Update control file to request self.data["nstate"] states. """
        if self.opts["ricc2"]:
            n_ex_state = self.data["nstate"] - 1
            file_utils.replace_inplace("control",
                                       r"(\s*irrep.*)nexc\s*=\s*\d+(.*)",
                                       r"\1nexc="+str(n_ex_state)+r"\2")
        if self.opts["egrad"]:
            n_ex_state = self.data["nstate"] - 1
            if n_ex_state == 0:
                return
            n_sub = file_utils.replace_inplace("control",
                    r"^\s*a\s+\d+\s*$",
                    rf" a  {n_ex_state}\n")
            if n_sub != 1:
                raise ValueError("Failed to update number of states.")

    def update_iroot(self):
        ex_state = self.data["iroot"]-1
        if self.opts["ricc2"]:
            if ex_state == 0:
                state = r"(x)"
            else:
                state = rf"(a {ex_state})"
            re_grad = re.escape(self.opts["ricc2_model"])
            re_grad = rf"(geoopt +model={re_grad} +state=).*"
            n_sub = file_utils.replace_inplace(self.template, re_grad, r"\1" + state)
            if n_sub == 0:
                repl = r"$ricc2\n  geoopt model={} state={}"
                repl = repl.format(self.opts["ricc2_model"], state)
                n_sub = file_utils.replace_inplace(self.template, r"\$ricc2", repl)
        elif self.opts["egrad"]:
            if ex_state == 0:
                return
            n_sub = file_utils.replace_inplace(self.template,
                    r"\$exopt.*", 
                    rf"$exopt {ex_state}")
            if n_sub == 0:
                file_utils.replace_inplace(self.template,
                        r"\$end", 
                        rf"$exopt {ex_state}\n$end")


def get_group(group):
    """ Get the data group from the actual output file. """
    sdg_run = subprocess.run(["sdg", group], capture_output=True, check=False)
    return sdg_run.returncode, sdg_run.stdout.decode().splitlines()


def actual_check():
    """ Check that the Turbomole calculation finished without error. """
    check = subprocess.run("actual", stdout=subprocess.PIPE, check=True)
    if check.stdout.startswith(b'fine, there is no data group "$actual step"'):
        return
    print("Turbomole calculation failed. Check output in " + os.getcwd() + ".")
    sys.exit(1)


def get_ricc2_energy(log_file, gs_model):
    """ Read energy from ricc2 output file."""
    try:
        re_gs = "Final " + re.escape(gs_model.upper()) + " energy"
        gs_energy = file_utils.search_file(log_file, re_gs)
        gs_energy = file_utils.split_columns(gs_energy, col=5, convert=np.float64)
    except ValueError:
        re_gs = rf"Method          :  {gs_model.upper()}"
        gs_energy = file_utils.search_file(log_file, re_gs, after=1)[0]
        energy = np.array(gs_energy.split()[3], dtype=float, ndmin=1)
        return energy
    ex_energy = file_utils.search_file(log_file, "Energy:")
    ex_energy = file_utils.split_columns(ex_energy, 1, convert=np.float64)
    energy = np.repeat(gs_energy, len(ex_energy) + 1)
    energy[1:] = energy[1:] + ex_energy
    return energy


def get_dft_energy(log_file):
    """ Read energy from dscf/ridft/escf/egrad output file."""
    try:
        # Look for ESCF output
        energy = file_utils.search_file(log_file, r"Total energy:")
        col = 2
    except ValueError:
        # Look for DSCF output
        re_en = re.escape(r"|  total energy      =")
        energy = file_utils.search_file(log_file, re_en)
        col = 4
    energy = file_utils.split_columns(energy, col=col, convert=np.float64)
    return np.array(energy)


def get_ricc2_gradient(log_file, iroot):
    """ Read gradient from ricc2 output file."""
    if iroot == 1:
        re_gs = r"GROUND STATE FIRST-ORDER PROPERTIES"
        cfile = file_utils.go_to_keyword(log_file, re_gs)[0]
    else:
        # For excited state gradients first go to properties section in ricc2 output.
        re_grad = r"EXCITED STATE PROPERTIES"
        cfile = file_utils.go_to_keyword(log_file, re_grad)[0]
        # Then look for the gradient of the correct state.
        re_grad = r"Excited state reached by transition:"
        while True:
            line = file_utils.search_file(cfile,
                    re_grad,
                    max_res=1,
                    close=False,
                    after=3)
            cstate = int(line[0].split()[4]) + 1
            if cstate == iroot:
                break
    grad = get_grad_from_stdout(cfile)
    cfile.close()
    return grad


def get_dft_gradient(log_file, iroot, ri):
    """ Read gradient from grad/rdgrad/egrad output files. """
    if iroot == 1:
        if ri:
            re_grad = r"RDGRAD - INFORMATION"
        else:
            re_grad = r"SCF ENERGY GRADIENT with respect to NUCLEAR COORDINATES"
    else:
        re_grad = f"Excited state no.{iroot:5d} chosen for optimization"
    cfile = file_utils.go_to_keyword(log_file, re_grad)[0]
    grad = get_grad_from_stdout(cfile)
    cfile.close()
    return grad


def get_grad_from_stdout(cfile):
    """ Read gradient from any Turbomole output file. """
    grad = file_utils.search_file(cfile, r"^  ATOM", after=3, 
                                  stop_at=r"resulting FORCE", close=False)
    grad = [line[5:] for line in grad]
    grad = [' '.join(grad[0::3]), ' '.join(grad[1::3]), ' '.join(grad[2::3])]
    grad = [line.split() for line in grad]
    grad = [list(map(file_utils.fortran_double, vals)) for vals in grad]
    return np.array(grad).T


def get_tddft_oscill(fname):
    """ Read oscillator strengths from escf/egrad output file. """
    oscill = file_utils.search_file(fname, "mixed representation:")
    file_utils.split_columns(oscill, col=2, convert=np.float64)
    return np.array(oscill)


def get_ricc2_oscill(fname):
    """ Read oscillator strengths from ricc2 output file. """
    oscill = file_utils.search_file(fname,
                                    r"oscillator strength \(length gauge\)")
    oscill = file_utils.split_columns(oscill, col=5)
    return np.array(oscill, dtype=np.float64)


def get_grad_from_gradient(natom):
    """ Get gradient from gradient file."""
    grad = file_utils.search_file("gradient", r"cycle", after=2*natom)[-natom:]
    grad = file_utils.split_columns(grad, col=[0, 1, 2],  convert=file_utils.fortran_double)
    return grad
