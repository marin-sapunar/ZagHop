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
    def __init__(self,
                 template="control",
                 log_file="qm.log",
                 err_file="qm.err",
                 **kwargs):
        """ Initialize the interface. """
        super().__init__(template=template,
                         log_file=log_file,
                         err_file=err_file,
                         **kwargs)
        self.read_funcs = {
            "energy": self.read_energy,
            "gradient": self.read_gradient,
            "oscillator_strength": self.read_oscill
        }
        check, lines = get_group("ricc2")
        if check != 0:
            self.options["ricc2"] = False
        else:
            self.options["ricc2"] = True
            for line in lines:
                if line.strip() in ["adc(2)", "mp2"]:
                    self.options["ricc2_model"] = line.strip()
                    break
            if "ricc2_model" not in self.options:
                print("Error in Turbomole interface.")
                print("Unrecognized wave function model in $ricc2 section.")
                sys.exit(1)
            # MP2 is the ground state for ADC(2) calculations
            if self.options["ricc2_model"] == "adc(2)":
                self.options["ricc2_gs_model"] = "mp2"
            else:
                self.options["ricc2_gs_model"] = self.options["ricc2_model"]
        check, lines = get_group("soes")
        if check != 0:
            self.options["egrad"] = False
        else:
            self.options["egrad"] = True
        if file_utils.check_in_file(template, r"\$rij"):
            self.options["ri"] = True
        else:
            self.options["ri"] = False
        self.options["qmmm"] = False

    def run(self):
        """ Run Turbomole calculation. """
        # Remove existing gradient file to avoid the file becoming huge.
        file_utils.remove("gradient")
        file_utils.remove(self.options["log_file"])
        file_utils.remove(self.options["err_file"])

        # SCF calculation.
        if self.options["ri"]:
            self.run_prog("ridft")
        else:
            self.run_prog("dscf")

        # Excited state + gradient calculation(s).
        if self.options["ricc2"]:
            # ricc2 does everything with a single call.
            self.run_prog("ricc2")
        else:
            # rdgrad, grad and egrad calculate the gradient.
            # egrad also calculates the excited states.
            # escf is called if excited state energies are requested without
            #   an excited state gradient.
            if "gradient" in self.request:
                if self.system["iroot"] == 1:
                    # Ground state gradient calculation.
                    if self.options["ri"]:
                        self.run_prog("rdgrad")
                    else:
                        self.run_prog("grad")
                    # Separate excited state calculation.
                    if self.system["nstate"] > 1:
                        self.run_prog("escf")
                elif self.system["iroot"] > 1:
                    self.run_prog("egrad")
            elif self.system["nstate"] > 1:
                self.run_prog("escf")

    def run_prog(self, prog):
        """ Run a Turbomole program and check it finished successfully."""
        # All subprocess calls use check=False since Turbomole programs don't
        # return a non-zero exit code on error anyways. Instead, to check for
        # errors `actual_check` calls "actual" and parses the output.
        with open(self.options["log_file"], "a") as out:
            with open(self.options["err_file"], "a") as err:
                subprocess.run(prog, stdout=out, stderr=err, check=False)
        actual_check(self.options["err_file"])

    def read_energy(self):
        """ Read energy from log file. """
        if self.options["ricc2"]:
            self.results["energy"] = get_ricc2_energy(
                self.options["log_file"], self.options["ricc2_gs_model"])
        else:
            self.results["energy"] = get_dft_energy(self.options["log_file"])

    def read_gradient(self):
        """ Read gradient from log file. """
        if self.options["ricc2"]:
            self.results["gradient"] = get_ricc2_gradient(
                self.options["log_file"], self.system["iroot"])
        else:
            self.results["gradient"] = get_dft_gradient(
                self.options["log_file"], self.system["iroot"],
                self.options["ri"])

    def read_oscill(self):
        """ Read oscillator strengths from log file. """
        if self.options["ricc2"]:
            self.results["oscillator_strength"] = get_ricc2_oscill(
                self.options["log_file"])
        else:
            self.results["oscillator_strength"] = get_tddft_oscill(
                self.options["log_file"])

    def update_geom(self):
        """ Update coord file with self.system["geom"]. """
        file_utils.replace_cols_inplace("coord", self.system["geom"],
                                        r"\$coord")
        if self.options["qmmm"]:
            pc_file = file_utils.search_file("control", r"\$point_charges")[0]
            pc_file = pc_file.split("=")[1].split()[0]
            file_utils.replace_cols_inplace(pc_file, self.system["mm_geom"],
                                            r"\$point_charges")

    def update_nstate(self):
        """ Update control file to request self.system["nstate"] states. """
        if self.options["ricc2"]:
            n_ex_state = self.system["nstate"] - 1
            file_utils.replace_inplace("control",
                                       r"(\s*irrep.*)nexc\s*=\s*\d+(.*)",
                                       r"\1nexc=" + str(n_ex_state) + r"\2")
        if self.options["egrad"]:
            n_ex_state = self.system["nstate"] - 1
            if n_ex_state == 0:
                return
            n_sub = file_utils.replace_inplace("control", r"^\s*a\s+\d+\s*$",
                                               rf" a  {n_ex_state}\n")
            if n_sub != 1:
                raise ValueError("Failed to update number of states.")

    def update_iroot(self):
        ex_state = self.system["iroot"] - 1
        if self.options["ricc2"]:
            if ex_state == 0:
                state = r"(x)"
            else:
                state = rf"(a {ex_state})"
            re_grad = re.escape(self.options["ricc2_model"])
            re_grad = rf"(geoopt +model={re_grad} +state=).*"
            n_sub = file_utils.replace_inplace(self.options["template"],
                                               re_grad, r"\1" + state)
            if n_sub == 0:
                repl = r"$ricc2\n  geoopt model={} state={}"
                repl = repl.format(self.options["ricc2_model"], state)
                n_sub = file_utils.replace_inplace(self.options["template"],
                                                   r"\$ricc2", repl)
        elif self.options["egrad"]:
            if ex_state == 0:
                return
            n_sub = file_utils.replace_inplace(self.options["template"],
                                               r"\$exopt.*",
                                               rf"$exopt {ex_state}")
            if n_sub == 0:
                file_utils.replace_inplace(self.options["template"], r"\$end",
                                           rf"$exopt {ex_state}\n$end")


def get_group(group):
    """ Get the data group from the actual output file. """
    sdg_run = subprocess.run(["sdg", group], capture_output=True, check=False)
    return sdg_run.returncode, sdg_run.stdout.decode().splitlines()


def actual_check(err_file):
    """ Check that the Turbomole calculation finished without error. """
    check_str = b'fine, there is no data group "$actual step"'
    check1 = file_utils.check_in_file(err_file, r"ended abnormally")
    check2 = subprocess.run("actual", stdout=subprocess.PIPE, check=True)
    check2 = not check2.stdout.startswith(check_str)
    if check1 or check2:
        errmsg = ("Turbomole calculation failed.\n"
                  "    Check output in " + os.getcwd() + ".")
        print(errmsg)
        sys.exit(1)


def get_ricc2_energy(log_file, gs_model):
    """ Read energy from ricc2 output file."""
    try:
        re_gs = "Final " + re.escape(gs_model.upper()) + " energy"
        gs_energy = file_utils.search_file(log_file, re_gs)
        gs_energy = file_utils.split_columns(gs_energy,
                                             col=5,
                                             convert=np.float64)
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


def get_dft_gradient(log_file, iroot, rdgrad):
    """ Read gradient from grad/rdgrad/egrad output files. """
    if iroot == 1:
        if rdgrad:
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
    grad = file_utils.search_file(cfile,
                                  r"^  ATOM",
                                  after=3,
                                  stop_at=r"resulting FORCE",
                                  close=False)
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
    grad = file_utils.search_file("gradient", r"cycle",
                                  after=2 * natom)[-natom:]
    grad = file_utils.split_columns(grad,
                                    col=[0, 1, 2],
                                    convert=file_utils.fortran_double)
    return grad
