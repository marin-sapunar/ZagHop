""" Base class for interfacing with electronic structure codes. """
import os.path
import sys


class QMInterface():
    """ 
    Generic interface for electronic structure calculations.

    The main methods of this class are `update_input`, `run`, and `read`. The `update_input`
    method is used to set up the calculation, `run` is used to execute the calculation, and
    `read` is used to extract the results from the output files.

    The derived classes should implement the following methods:
    - `run`: Execute the QM calculation.
    - `read_energy`: Read the energies from the QM output.
    - `read_gradient`: Read the gradient from the QM output.
    - `read_oscill`: Read the oscillator strengths from the QM output.
    - `update_geom`: Update the geometry in the QM input files.
    - `update_nstate`: Update the number of states requested from the QM code.
    - `update_iroot`: Request the gradient for a specified state.

    Methods starting with `update_` are called by the `update_input` method and should update
    the input files for the QM calculation based on the provided geometry and options.
    
    Methods starting with `read_` are called by the `read` method and should store the
    calculated values in the `data` dictionary.
        
    Attributes:
        template (str): The name of the template file for the QM calculation.
        log_file (str): The name of the log file to store the QM output.
        err_file (str): The name of the error file to store any error messages.
        opts (dict): A dictionary of additional options for the QM calculation.
        data (dict): A dictionary to store the calculated data from the QM output.
        need (dict): A dictionary specifying the values to be calculated/read.
        calc_done (bool): A flag indicating whether the calculation has been completed.
    """
    def __init__(self, template="qm.inp", log_file="qm.log", err_file="qm.err", **kwargs):
        """
        Initialize the QM interface object.
        
        Parameters:
            - template (str): The path to the template file. Default is "qm.inp".
            - log_file (str): The path to the log file. Default is "qm.log".
            - err_file (str): The path to the error file. Default is "qm.err".
            - **kwargs: Additional keyword arguments.
        
        Raises:
            - SystemExit: If the template file is not found.
        """
        self.template = template
        self.log_file = log_file
        self.err_file = err_file
        self.calc_done = False
        self.data = {}
        self.opts = {}
        self.need = {}
        self.opts.update(kwargs)
        if not os.path.isfile(self.template):
            print("Error in QM interface. Template file not found")
            sys.exit(1)

    def update_input(self, in_data, request):
        """ 
        Update input files to prepare a calculation. 
        
        Parameters:
            - in_data: The input data for the calculation.
            - request: A dictionary specifying values to be calculated/read.
        """
        self.data.update(in_data)
        self.need.update(request)
        self.update_geom()
        self.update_nstate()
        if "gradient" in self.need:
            self.update_iroot()
        self.calc_done = False

    def run(self):
        """ Run the calculation and check success. """
        raise NotImplementedError("Need to call specific interface.")

    def read(self):
        """ 
        Read all recognized values from "request" dictionary.
        
        Returns:
            The data dictionary containing the calculated values.
        """
        self.read_energy()
        if self.need["gradient"]:
            self.read_gradient()
        if self.need["oscillator_strength"]:
            self.read_oscill()
        return self.data

    def read_energy(self):
        """ 
        Read energies from a completed calculation.
        
        The energies should be stored in the `energy` field of the `data` dictionary.
        """
        raise NotImplementedError("Need to call specific interface.")

    def read_gradient(self):
        """ 
        Read gradient from a completed calculation.
        
        The gradient should be stored in the `gradient` field of the `data` dictionary.
        """
        raise NotImplementedError("Need to call specific interface.")

    def read_oscill(self):
        """ 
        Read oscillator strengths from a completed calculation.
        
        The oscillator strengths should be stored in the `oscillator_strength` field
        of the `data` dictionary.
        """
        raise NotImplementedError("Need to call specific interface.")

    def update_geom(self):
        """ Update geometry in the QM input files. """
        raise NotImplementedError("Need to call specific interface.")

    def update_nstate(self):
        """ Update number of states requested from the QM code. """
        raise NotImplementedError("Need to call specific interface.")

    def update_iroot(self):
        """ Request gradient for specified state. """
        raise NotImplementedError("Need to call specific interface.")
