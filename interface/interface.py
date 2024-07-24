""" Base class for interfacing with electronic structure codes. """
import os.path
import sys


class QMInterface():
    """ 
    Generic interface for electronic structure calculations.

    The main methods of this class are `update_input`, `run`, and `read`. The
    `update_input` method is used to set up the calculation, `run` is used to
    execute the calculation, and `read` is used to extract the requested
    results from the output files.

    Derived classes should implement the following methods:
    - `__init__`: Set up everything needed for the interface to run the other
                  methods. Specifically, it should always define the read_funcs
                  dictionary.
    - `run`: Execute the QM calculation.
    - `update_geom`: Update the geometry in the QM input files.
    - `update_nstate`: Update the number of states requested from the QM code.
    - `update_iroot`: Update the target state (for gradient calculation).
    - `read_*` : Methods to be called to read values from the QM output.

    Methods starting with `update_` are called by the `update_input` method and
    should update the input files for the QM calculation based on the provided 
    geometry and options.
    
    Methods starting with `read_` are called by the `read` method and should
    store the calculated values in the `results` dictionary.
        
    Attributes:
        read_funcs (dict) : Key, value pairs of values which can be read by
                            the interface and the corresponding read functions.
        options (dict): Options and paths relevant for the QM calculation.
        system (dict): Information on the system.
        request (list): Values which should be calculated/read.
        results (dict): Results read from the QM output.
        calc_done (bool): True after a calculation has been completed.
    """
    def __init__(self,
                 template="qm.inp",
                 log_file="qm.log",
                 err_file="qm.err",
                 **kwargs):
        """
        Initialize the QM interface object.
        
        Parameters:
            - template (str): Path to the template file. Default is "qm.inp".
            - log_file (str): Path to the log file. Default is "qm.log".
            - err_file (str): Path to the error file. Default is "qm.err".
            - **kwargs: Additional keyword arguments.
        
        Raises:
            - SystemExit: If the template file is not found.
        """
        self.calc_done = False
        self.read_funcs = {}
        self.system = {}
        self.results = {}
        self.request = []
        self.options = {
            "template": template,
            "log_file": log_file,
            "err_file": err_file
        }
        self.options.update(kwargs)
        if not os.path.isfile(self.options["template"]):
            print("Error in QM interface. Template file not found.")
            sys.exit(1)

    def update_input(self, in_data, request):
        """ 
        Update input files to prepare a calculation. 
        
        Parameters:
            - in_data: The input data for the calculation.
            - request: List of values to be calculated/read.
        """
        self.system.update(in_data)
        self.request = request
        self.update_geom()
        self.update_nstate()
        self.update_iroot()
        self.calc_done = False

    def run(self):
        """ Run the calculation and check success. """
        raise NotImplementedError("Need to call specific interface.")

    def read(self):
        """ 
        Read all requested values.
        
        Returns:
            The data dictionary containing the calculated values.
        """
        for req in self.request:
            if req not in self.read_funcs:
                print(f"Error in QM interface. {req} not available.")
                sys.exit(1)
            self.read_funcs[req]()
        return self.results

    def update_geom(self):
        """ Update geometry in the QM input files. """
        raise NotImplementedError("Need to call specific interface.")

    def update_nstate(self):
        """ Update number of states requested from the QM code. """
        raise NotImplementedError("Need to call specific interface.")

    def update_iroot(self):
        """ Request gradient for specified state. """
        raise NotImplementedError("Need to call specific interface.")
