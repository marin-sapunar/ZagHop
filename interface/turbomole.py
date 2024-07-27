""" Interface for running Turbomole RICC2 or TDDFT calculations. """
import subprocess
import re
import os
import sys
import numpy as np
from states import ElStates
import file_utils
from interface import QMInterface


class Turbomole(QMInterface):
    """ Interface for turbomole calculations. """
    def __init__(self,
                 template="control",
                 log_file="qm.log",
                 err_file="qm.err",
                 **kwargs):
        """ Init of base class, only defaults changed. """
        super().__init__(template=template,
                         log_file=log_file,
                         err_file=err_file,
                         **kwargs)

    @classmethod
    def generate_inputs(cls, system, options):
        """ Run define to generate inputs for Turbomole. """
        inst = cls(**options)
        states = ElStates(system["states"])
        write_coord("coord", system["geom"], system["atoms"])
        define_inp = generate_define_input(inst.options, states)
        define_run = subprocess.run("define",
                                    input=define_inp.encode(),
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    check=False)
        # @todo Add checks to see if everything was set correctly
        return inst

    def check_template(self):
        """ Read options from template. """
        check, lines = get_group("ricc2")
        self.options["ricc2"] = check
        if self.ricc2:
            for line in lines:
                if line.strip() in ["adc(2)", "mp2"]:
                    self.options["model"] = line.strip()
                    break
            if "model" not in self.options:
                print("Error in Turbomole interface.")
                print("Unrecognized wave function model in $ricc2 section.")
                sys.exit(1)
            # MP2 is the ground state for ADC(2) calculations
            if self.model == "adc(2)":
                self.options["gs_model"] = "mp2"
            else:
                self.options["gs_model"] = self.model
        check, lines = get_group("soes")
        self.options["egrad"] = check
        check = get_group("rij")[0]
        self.options["ri"] = check
        if self.ricc2:
            self.read_funcs = {
                "energy": self.read_energy_ricc2,
                "gradient": self.read_gradient_ricc2,
                "oscillator_strength": self.read_oscill_ricc2
            }
        else:
            self.read_funcs = {
                "energy": self.read_energy_dft,
                "gradient": self.read_gradient_dft,
                "oscillator_strength": self.read_oscill_dft
            }

    def update(self, system, request):
        self.system = system
        self.system["states"] = ElStates(self.system["states"])
        self.request = request
        # Update coord file.
        file_utils.replace_cols_inplace("coord", self.geom, r"\$coord")
        # Update number of states.
        self.update_nstate()
        # Update iroot
        if "gradient" in self.request:
            self.options["iroot"] = self.request["gradient"]
            self.update_iroot()
        else:
            self.options["iroot"] = None

    def run(self):
        """ Run Turbomole calculation. """
        # Remove existing gradient file to avoid the file becoming huge.
        file_utils.remove("gradient")
        file_utils.remove(self.log_file)
        file_utils.remove(self.err_file)

        # SCF calculation.
        if self.ri:
            self.run_prog("ridft")
        else:
            self.run_prog("dscf")

        # Excited state + gradient calculation(s).
        if self.ricc2:
            # ricc2 does everything with a single call.
            self.run_prog("ricc2")
        else:
            # rdgrad, grad and egrad calculate the gradient.
            # egrad also calculates the excited states.
            # escf is called if excited state energies are requested without
            #   an excited state gradient.
            if "gradient" in self.request:
                if self.request["gradient"] == 1:
                    # Ground state gradient calculation.
                    if self.ri:
                        self.run_prog("rdgrad")
                    else:
                        self.run_prog("grad")
                    # Separate excited state calculation.
                    if self.states.nstate > 1:
                        self.run_prog("escf")
                elif self.request["gradient"] > 1:
                    self.run_prog("egrad")
            elif self.states.nstate > 1:
                self.run_prog("escf")

    def run_prog(self, prog):
        """ Run a Turbomole program and check it finished successfully."""
        # All subprocess calls use check=False since Turbomole programs don't
        # return a non-zero exit code on error anyways. Instead, to check for
        # errors `actual_check` calls "actual" and parses the output.
        with open(self.log_file, "a") as out:
            with open(self.err_file, "a") as err:
                subprocess.run(prog, stdout=out, stderr=err, check=False)
        actual_check(self.err_file)

    def read_energy_ricc2(self):
        """ Read energy from log file. """
        self.results["energy"] = get_ricc2_energy(self.log_file, self.gs_model)

    def read_energy_dft(self):
        """ Read energy from log file. """
        self.results["energy"] = get_dft_energy(self.log_file)

    def read_gradient_ricc2(self):
        """ Read gradient from log file. """
        self.results["gradient"] = get_ricc2_gradient(self.log_file,
                                                      self.iroot)

    def read_gradient_dft(self):
        """ Read gradient from log file. """
        self.results["gradient"] = get_dft_gradient(self.log_file, self.iroot,
                                                    self.ri)

    def read_oscill_ricc2(self):
        """ Read oscillator strengths from log file. """
        self.results["oscillator_strength"] = get_ricc2_oscill(self.log_file)

    def read_oscill_dft(self):
        """ Read oscillator strengths from log file. """
        self.results["oscillator_strength"] = get_tddft_oscill(self.log_file)

    def update_nstate(self):
        """ Update control file to request correct number of states."""
        n_exci = self.states.nstate - 1
        if n_exci == 0:
            return
        if self.ricc2:
            n_sub = file_utils.replace_inplace(
                "control", r"(\s*irrep.*)nexc\s*=\s*\d+(.*)",
                r"\1nexc=" + str(n_exci) + r"\2")
        if self.egrad:
            n_sub = file_utils.replace_inplace("control", r"^\s*a\s+\d+\s*$",
                                               rf" a  {n_exci}\n")
        if n_sub != 1:
            raise ValueError("Failed to update number of states.")

    def update_iroot(self):
        """ Update control file to request gradient of specific state. """
        eiroot = self.iroot - 1
        if self.ricc2:
            if eiroot == 0:
                state = r"(x)"
            else:
                state = rf"(a {eiroot})"
            re_grad = re.escape(self.model)
            re_grad = rf"(geoopt +model={re_grad} +state=).*"
            n_sub = file_utils.replace_inplace(self.template, re_grad,
                                               r"\1" + state)
            if n_sub == 0:
                repl = rf"$ricc2\n  geoopt model={self.model} state={state}"
                n_sub = file_utils.replace_inplace(self.template, r"\$ricc2",
                                                   repl)
        elif self.egrad:
            if eiroot == 0:
                return
            n_sub = file_utils.replace_inplace(self.template, r"\$exopt.*",
                                               rf"$exopt {eiroot}")
            if n_sub == 0:
                file_utils.replace_inplace(self.template, r"\$end",
                                           rf"$exopt {eiroot}\n$end")


def get_group(group):
    """ Get the data group from the actual output file. """
    sdg_run = subprocess.run(["sdg", group], capture_output=True, check=False)
    check = sdg_run.returncode == 0
    return check, sdg_run.stdout.decode().splitlines()


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


def write_coord(fname, geom, atoms):
    """ Writes a coord file for Turbomole input. Takes output file name,
        numpy array containing geometry and array containing atoms and
        creates a coord file. """
    natom = len(atoms)
    gformat = 3 * '{:20.13f} ' + '  {}\n'
    with open(fname, 'w') as ofile:
        ofile.write('$coord\n')
        for atom, xyz in zip(atoms, geom):
            ofile.write(gformat.format(*xyz, atom.lower()))
        ofile.write('$end\n')


def generate_define_input(options, states):
    """ String for input to `define`. """
    nstate = states.nstate
    define_str = "\n\na coord\n*\nno\n"  # Geometry definition
    define_str += f'b all {options["basis"]}\n*\n'  # Basis set definition
    define_str += "eht\n\n\n\n"  # MO guess
    define_str += "scf\nconv\n7\n\n"
    if options["method"] == "ricc2":
        define_str += "cc\n"  # Open ricc2 submenu
        define_str += "freeze\n*\ncbas\n*\n"  # Default freeze + cbas
        define_str += f'ricc2\n{options["model"]}\n*\n'
        if nstate > 1:
            define_str += "exci\n"
            define_str += f"irrep=a nexc={nstate}\n"
            define_str += "*\n"
        define_str += "*\n"
    else:
        define_str += f'dft\non\nfunc {options["model"]}\n\n'
        if nstate > 1:
            define_str += f"ex\nrpas\n*\na {nstate}\n*\nrpaconv 6\n*\n\n"
    define_str += "*\n"
    return define_str
